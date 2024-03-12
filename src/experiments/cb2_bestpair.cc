#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <ilcplex/cplex.h>

#include <inputs.hh>
#include <live.hh>

// parameters

// seems faster on bigger problems
constexpr bool PARAM_DYNAMIC_SEARCH = true;
// greatly improves the LP
constexpr bool PARAM_LOCAL_L = true;
// makes nodes evaluation about 3 to 5 times slower, does not seem to improve the bound
constexpr bool PARAM_LOCAL_M = true;
// -1 to disable
constexpr double PARAM_TIME_LIMIT = 45;

constexpr bool PARAM_LIVE_SOL = false;

#define _c(what) if (int _error = what) { \
	cout << "CPX error: " #what << endl; cout << "CPX error: " << _error << endl; abort(); }

using namespace std;

int B_order[MAX_N][MAX_N];
int B_jk_order[MAX_N][MAX_N][MAX_N];

int i_x(int i) {
	assert(0 <= i && i < N);
	return i;
}

int i_w(int i) {
	assert(0 <= i && i < N);
	return N + i;
}

void preorder_B() {
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) B_order[i][j] = j;
		sort(B_order[i], B_order[i]+N, [i](int ja, int jb) {
			return B[i][ja] < B[i][jb];
		});
	}
	for (int j=0; j<N; j++)
	for (int k=0; k<N; k++) {
		for (int i=0; i<N; i++) B_jk_order[j][k][i] = i;
		sort(B_jk_order[j][k], B_jk_order[j][k]+N, [j,k](int ia, int ib) {
			return B[ia][j]+B[ia][k] < B[ib][j]+B[ib][k];
		});
	}
}

void add_variables_x(CPXENVptr env, CPXLPptr lp) {
	vector<double> ub(N, 1.0);
	vector<char> type(N, 'I');
	vector<string> name(N);
	vector<char*> c_name(N);
	for (int i=0; i<N; i++) {
		name[i] = "x_" + to_string(i);
		c_name[i] = const_cast<char*>(name[i].c_str());
	}
	_c(CPXnewcols(env, lp, N, NULL, NULL, ub.data(), type.data(), c_name.data()));
}

void add_variables_w(CPXENVptr env, CPXLPptr lp) {
	vector<double> lb(N, 0.0);
	vector<double> obj(N, 1.0);
	vector<string> name(N);
	vector<char*> c_name(N);
	for (int i=0; i<N; i++) {
		name[i] = "w_" + to_string(i);
		c_name[i] = const_cast<char*>(name[i].c_str());
	}
	_c(CPXnewcols(env, lp, N, obj.data(), NULL, NULL, NULL, c_name.data()));
}

// select M locations constraint
// sum_i ( x_i ) = M
void add_cons_M(CPXENVptr env, CPXLPptr lp) {
    double rhs[1] = {(double) M};
    int rmatbeg[1] = {0};
	vector<int> rmatind(N);
	vector<double> rmatval(N, 1.0);
    for (int i = 0; i < N; i++) {
        rmatind[i] = i_x(i);
    }
    _c(CPXaddrows(env, lp, 0, 1, N, rhs, "E", rmatbeg, rmatind.data(), rmatval.data(), NULL, NULL));
}

double calc_Mi(int i) {
	int take = M;
	double tot = 0;
	for (int k = N-1; take && k >= 0; k--) {
		int j = B_order[i][k];
		if (j == i) continue; // fixed x_i = 0
		tot += B[i][j];
		take--;
	}
	return tot;
}

// w_i := { 0 if x_i=0; sum_j( b_ij x_j ) if x_j !=0 }
// => w_i >= sum_j( b_ij x_j ) - M_i * (1 - x_i)
// => sum_{j!=i}( b_ij x_j ) + (M_i + b_ii) x_i - w_i <= M_i
void add_cons_xw(CPXENVptr env, CPXLPptr lp) {
	vector<double> rhs(N, 0);
	vector<char> sense(N, 'L');
	vector<int> rmatbeg(N, 0);
	vector<int> rmatind(N*(N+1), 0);
	vector<double> rmatval(N*(N+1), 0);
    for (int i = 0; i < N; i++) {
		double Mi = calc_Mi(i);
		int beg = i * (N+1);
		rhs[i] = Mi;
		rmatbeg[i] = beg;
		for (int j = 0; j < N; j++) {
			rmatind[beg + j] = i_x(j);
			rmatval[beg + j] = B[i][j] + (i==j ? Mi : 0);
		}
		rmatind[beg + N] = i_w(i);
		rmatval[beg + N] = -1;
	}
    _c(CPXaddrows(env, lp, 0, N, N*(N+1), rhs.data(), sense.data(), rmatbeg.data(), rmatind.data(), rmatval.data(), NULL, NULL));
}

double calc_Li(int i) {
	int take = M-1;
	double tot = B[i][i]; // fixed x_i = 1
	for (int k = 0; take && k < N; k++) {
		int j = B_order[i][k];
		if (j == i) continue; // already counted
		tot += B[i][j];
		take--;
	}
	return tot;
}

// w_i >= L_i x_i  <=>  w_i - L_i x_i >= 0 
void add_cons_L(CPXENVptr env, CPXLPptr lp) {
	vector<double> rhs(N, 0);
	vector<char> sense(N, 'G');
	vector<int> rmatbeg(N, 0);
	vector<int> rmatind(N*2, 0);
	vector<double> rmatval(N*2, 0);
    for (int i = 0; i < N; i++) {
		double Li = calc_Li(i);
		int beg = i * 2;
		rmatbeg[i] = beg;
		rmatind[beg + 0] = i_w(i);
		rmatval[beg + 0] = 1;
		rmatind[beg + 1] = i_x(i);
		rmatval[beg + 1] = -Li;
	}
    _c(CPXaddrows(env, lp, 0, N, N*2, rhs.data(), sense.data(), rmatbeg.data(), rmatind.data(), rmatval.data(), NULL, NULL));
}

static double L_rhs[MAX_N];
static char L_sense[MAX_N];
static int L_rmatbeg[MAX_N];
static int L_rmatind[MAX_N*2];
static int L_purgeable[MAX_N];
static int L_local[MAX_N];

void init_local_L_data() {
	fill_n(L_rhs, N, 0);
	fill_n(L_sense, N, 'G');
	fill_n(L_rmatbeg, N, 0);
	fill_n(L_purgeable, N, CPX_USECUT_FORCE);
	fill_n(L_local, N, 1);
    for (int i = 0; i < N; i++) {
		int beg = i * 2;
		L_rmatbeg[i] = beg;
		L_rmatind[beg + 0] = i_w(i);
		L_rmatind[beg + 1] = i_x(i);
	}
}

void improve_local_L(CPXCALLBACKCONTEXTptr context, bool* fixed0, bool* fixed1) {
	// calculate all local L_i
	thread_local static double ll[MAX_N];
	for (int i=0; i<N; i++) {
		int take = M;
		double tot = 0;
		// always take i
		take--;
		tot += B[i][i];
		// take all the fixed variables
		for (int j = 0; take && j < N; j++) {
			if (j == i) continue; // already taken
			if (fixed1[j]) {
				tot += B[i][j];
				take--;
			}
		}
		// take the smallest remaining values
		for (int k = 0; take && k < N; k++) {
			int j = B_order[i][k];
			if (j == i) continue; // already taken
			if (fixed1[j]) continue; // already counted
			if (fixed0[j]) continue; // cannot be used
			tot += B[i][j];
			take--;
		}
		assert(take == 0);
		ll[i] = tot;
	}

	// generate new contraints
	thread_local static double L_rmatval[MAX_N*2];
    for (int i = 0; i < N; i++) {
		int beg = i * 2;
		L_rmatval[beg + 0] = 1;
		L_rmatval[beg + 1] = -ll[i];
	}
    _c(CPXcallbackaddusercuts(context, N, N*2, L_rhs, L_sense, L_rmatbeg, L_rmatind, L_rmatval, L_purgeable, L_local));

	// calculate all local L_jk
	thread_local static double lljk[MAX_N][MAX_N];
	int up = 0, down = 0;
	for (int j=0; j<N; j++) if (fixed1[j]) for (int k=0; k<N; k++) {
		if (fixed1[k] || fixed0[k]) continue;
		int take = M;
		double tot = 0;
		// take all the fixed variables
		for (int i = 0; take && i < N; i++) {
			if (fixed1[i]) {
				tot += B[i][j] + B[i][k];
				take--;
			}
		}
		// take the smallest remaining values
		for (int ord = 0; take && ord < N; ord++) {
			int i = B_jk_order[j][k][ord];
			if (fixed1[i]) continue; // already counted
			if (fixed0[i]) continue; // cannot be taken
			tot += B[i][j] + B[i][k];
			take--;
		}
		assert(take == 0);
		lljk[j][k] = tot;
		if (tot > ll[j]+ll[k]) up++;
		else down++;
	}

	// generate new Ljk contraints
	int nn_max = N;
	vector<double> rhs(nn_max, 0);
	vector<char> sense(nn_max, 'G');
	vector<int> rmatbeg(nn_max, 0);
	vector<int> rmatind(nn_max*3, 0);
	vector<double> rmatval(nn_max*3, 0);
	vector<int> purgeable(nn_max, CPX_USECUT_FORCE);
	vector<int> local(nn_max, 1);
	int nn = 0, beg = 0;
	for (int k=0; k<N; k++) {
		if (fixed0[k] || fixed1[k]) continue;
		int bj = -1;
		double best_l = -1;
		for (int j=0; j<N; j++) if (fixed1[j]) {
			if (lljk[j][k] > best_l) {
				bj = j;
				best_l = lljk[j][k];
			}
		}
		if (bj == -1) continue;
		double Ljk = lljk[bj][k];
		rmatbeg[nn] = beg;
		rmatind[beg  ] = i_w(bj);
		rmatval[beg++] = 1;
		rmatind[beg  ] = i_w(k);
		rmatval[beg++] = 1;
		rmatind[beg  ] = i_x(k);
		rmatval[beg++] = -Ljk;
		nn++;
	}
	assert(nn == nn);
	assert(beg == nn*3);
    _c(CPXcallbackaddusercuts(context, nn, nn*3, rhs.data(), sense.data(), rmatbeg.data(), rmatind.data(), rmatval.data(), purgeable.data(), local.data()));
}

static char M_sense[MAX_N];
static int M_rmatbeg[MAX_N];
static int M_rmatind[MAX_N*(MAX_N+1)];
static int M_purgeable[MAX_N];
static int M_local[MAX_N];

void init_local_M_data() {
	fill_n(M_sense, N, 'L');
	fill_n(M_rmatbeg, N, 0);
	fill_n(M_rmatind, N*(N+1), 0);
	fill_n(M_purgeable, N, CPX_USECUT_FORCE);
	fill_n(M_local, N, 1);
    for (int i = 0; i < N; i++) {
		int beg = i * (N+1);
		M_rmatbeg[i] = beg;
		for (int j = 0; j < N; j++) {
			M_rmatind[beg + j] = i_x(j);
		}
		M_rmatind[beg + N] = i_w(i);
	}
}

void improve_local_M(CPXCALLBACKCONTEXTptr context, bool* fixed0, bool* fixed1) {
	// calculate all local M_i
	thread_local static double lm[MAX_N];
	for (int i=0; i<N; i++) {
		int take = M;
		double tot = 0;
		// take the biggest values fixed to 1
		for (int k = N-1; take && k >= 0; k--) {
			int j = B_order[i][k];
			if (fixed1[j]) {
				tot += B[i][j];
				take--;
			}
		}
		// take the biggest remaining values
		for (int k = N-1; take && k >= 0; k--) {
			int j = B_order[i][k];
			if (j == i) continue; // never take x_i
			if (fixed1[j]) continue; // already taken
			if (fixed0[j]) continue; // cannot be taken
			tot += B[i][j];
			take--;
		}
		assert(take == 0);
		lm[i] = tot;
	}

	// generate new constraints
	thread_local static double M_rhs[MAX_N];
	thread_local static double M_rmatval[MAX_N*(MAX_N+1)];
    for (int i = 0; i < N; i++) {
		int beg = i * (N+1);
		M_rhs[i] = lm[i];
		for (int j = 0; j < N; j++) {
			M_rmatval[beg + j] = B[i][j] + (i==j ? lm[i] : 0);
		}
		M_rmatval[beg + N] = -1;
	}
    _c(CPXcallbackaddusercuts(context, N, N*(N+1), M_rhs, M_sense, M_rmatbeg, M_rmatind, M_rmatval, M_purgeable, M_local));
}

int cuts_generator(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle) {
	if (PARAM_LIVE_SOL) live_display(context);

	// producing cuts more then once for the same node is wasteful
	thread_local static long long prev_id = -1;
	long long node_id;
	_c(CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEUID, &node_id));
	if (node_id == prev_id) return 0;
	else prev_id = node_id;

	// lower/upperbound of the x variables at the local node
	thread_local static double lb[MAX_N];
	_c(CPXcallbackgetlocallb((CPXCALLBACKCONTEXTptr)context, lb, i_x(0), i_x(N-1)));
	thread_local static double ub[MAX_N];
	_c(CPXcallbackgetlocalub((CPXCALLBACKCONTEXTptr)context, ub, i_x(0), i_x(N-1)));
	// keeps track of which variables are fixed to 0/1
	thread_local static bool fixed1[MAX_N];
	for (int i=0; i<N; i++) fixed1[i] = lb[i] > .5;
	int num_fixed1 = count(fixed1, fixed1+N, true);
	thread_local static bool fixed0[MAX_N];
	for (int i=0; i<N; i++) fixed0[i] = ub[i] < .5;
	int num_fixed0 = count(fixed0, fixed0+N, true);
	assert(num_fixed0 <= N-M);

	if (PARAM_LOCAL_L)
		improve_local_L(context, fixed0, fixed1);
	if (PARAM_LOCAL_M)
		improve_local_M(context, fixed0, fixed1);
	
	return 0;
}

int main() {
	read_inputs();

	preorder_B();
	init_local_L_data();
	init_local_M_data();

	int status;
    CPXENVptr env = CPXopenCPLEX(&status);
	if (status) abort();
    CPXLPptr lp = CPXcreateprob(env, &status, "qap");
	if (status) abort();

	// enable solver logging
	_c(CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON));
	_c(CPXsetintparam(env, CPXPARAM_Threads, 1));

    _c(CPXchgobjsen(env, lp, CPX_MIN));

	add_variables_x(env, lp);
	add_variables_w(env, lp);

	add_cons_M(env, lp);
	add_cons_xw(env, lp);
	if (!PARAM_LOCAL_L) // not needed as they are added locally
		add_cons_L(env, lp);

	// write problem
	//_c(CPXwriteprob(env, lp, "p.lp", NULL));

	// add a callback to generate cuts
	_c(CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_RELAXATION, cuts_generator, NULL));

    // solve as mip
	if (PARAM_TIME_LIMIT != -1)
		_c(CPXsetdblparam(env, CPXPARAM_TimeLimit, PARAM_TIME_LIMIT));
	if (!PARAM_DYNAMIC_SEARCH)
		_c(CPXsetintparam(env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL));
	double t_start, t_end;
	_c(CPXgettime(env, &t_start));
    _c(CPXmipopt(env, lp));
	_c(CPXgettime(env, &t_end));

    // print solution
    double objval, best_boud;
	int solstat;
	vector<double> sol(N + N);
	_c(CPXsolution(env, lp, &solstat, &objval, sol.data(), NULL, NULL, NULL));
	_c(CPXgetbestobjval(env, lp, &best_boud));

	cout << "status = " << solstat << endl;
	cout << "solution = " << objval << endl;
	cerr << N << endl;
	for (int i=0; i<N; i++) {
		cerr << int(sol[i_x(i)]) << endl;
	}
	cerr << "# obj = " << objval << endl;
	if (solstat != CPXMIP_OPTIMAL)
		cerr << "# bound = " << best_boud << endl;
	cerr << "# time = " << (t_end - t_start) << endl;
	cerr << "# cb2 best pair"
		<< ", dynamic search = " << PARAM_DYNAMIC_SEARCH
		<< ", local L = " << PARAM_LOCAL_L
		<< ", local M = " << PARAM_LOCAL_M
		<< ", tl = " << PARAM_TIME_LIMIT << endl;
	cerr << "# n1 = " << M << endl;

    // clean up
    _c(CPXfreeprob(env, &lp));
    _c(CPXcloseCPLEX(&env));

    return 0;
} 
