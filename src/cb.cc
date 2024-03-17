#include <algorithm>
#include <atomic>
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
// makes nodes enumeration slower, slightly improves the LP
constexpr bool PARAM_LOCAL_L_PAIRS = false;
// makes nodes enumeration slower, slightly improves the LP
constexpr bool PARAM_LOCAL_L_ALL = false;
// makes nodes enumeration slighlty slower, improves the LP
constexpr bool PARAM_LOCAL_M = true;
// whether cuts are calculated only once per node, even after a new relaxation
constexpr bool PARAM_CUT_ONCE = false;
// display nodes as they are explored
constexpr bool PARAM_LIVE_SOL = false;
// don't bother adding cuts if less then a minimum are found
constexpr bool PARAM_CUTS_MIN = 4;
// -1 to disable
constexpr double PARAM_TIME_LIMIT = -1;
// single vs multi thread
constexpr bool PARAM_SINGLE_THREAD = false;

constexpr double EPS = 1; // add cuts that are violated by this much
constexpr int CUT_TYPE = CPX_USECUT_FORCE;

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
}

void preorder_B_pairs() {
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
static int L_purgeable[MAX_N];
static int L_local[MAX_N];

void init_L_data() {
	fill_n(L_rhs, N, 0);
	fill_n(L_sense, N, 'G');
	fill_n(L_rmatbeg, N, 0);
	fill_n(L_purgeable, N, CUT_TYPE);
	fill_n(L_local, N, 1);
    for (int i = 0; i < N; i++) {
		int beg = i * 2;
		L_rmatbeg[i] = beg;
	}
}

double calc_local_L(int i, bool* fixed0, bool* fixed1) {
	int take = M;
	double tot = 0;
	// take all the fixed variables
	for (int j = 0; j < N; j++) {
		if (j==i || fixed1[j]) {
			tot += B[i][j];
			take--;
		}
	}
	assert(take >= 0);
	// take the smallest remaining values
	for (int k = 0; take && k < N; k++) {
		int j = B_order[i][k];
		if (j==i || fixed1[j]) continue; // already counted
		if (fixed0[j]) continue; // cannot be used
		tot += B[i][j];
		take--;
	}
	assert(take == 0);
	return tot;
}

void improve_local_L(CPXCALLBACKCONTEXTptr context, bool* fixed0, bool* fixed1, double* x, double* w) {
	thread_local static double L_rmatval[MAX_N*2];
	thread_local static int L_rmatind[MAX_N*2];
	int nn = 0;
    for (int i = 0; i < N; i++) {
		if (fixed0[i] || fixed1[i]) continue;
		double ll = calc_local_L(i, fixed0, fixed1);
		double delta = w[i] - ll * x[i];
		if (delta >= -EPS) continue;
		int beg = nn * 2;
		L_rmatind[beg + 0] = i_w(i);
		L_rmatval[beg + 0] = 1;
		L_rmatind[beg + 1] = i_x(i);
		L_rmatval[beg + 1] = -ll;
		nn++;
	}
	if (nn >= PARAM_CUTS_MIN)
    	_c(CPXcallbackaddusercuts(context, nn, nn*2, L_rhs, L_sense, L_rmatbeg, L_rmatind, L_rmatval, L_purgeable, L_local));
}

static char Lp_sense[MAX_N*MAX_N];
static int Lp_rmatbeg[MAX_N*MAX_N];
static int Lp_purgeable[MAX_N*MAX_N];
static int Lp_local[MAX_N*MAX_N];

void init_Lp_data() {
	fill_n(Lp_sense, N*N, 'G');
	fill_n(Lp_rmatbeg, N*N, 0);
	fill_n(Lp_purgeable, N*N, CUT_TYPE);
	fill_n(Lp_local, N*N, 1);
	for (int i=0; i<N*N; i++) Lp_rmatbeg[i] = i * 3;
}

double calc_local_Ljk(int j, int k, bool* fixed0, bool* fixed1) {
	int take = M;
	double tot = 0;
	// take all the fixed variables
	for (int i = 0; i < N; i++) {
		if (fixed1[i] || i==k) {
			tot += B[i][j] + B[i][k];
			take--;
		}
	}
	assert(take >= 0);
	// take the smallest remaining values
	for (int ord = 0; take && ord < N; ord++) {
		int i = B_jk_order[j][k][ord];
		if (fixed1[i] || i==k) continue; // already counted
		if (fixed0[i]) continue; // cannot be taken
		tot += B[i][j] + B[i][k];
		take--;
	}
	assert(take == 0);
	return tot;
}

static atomic_int tot_lp_cuts = 0;

void improve_local_L_pairs(CPXCALLBACKCONTEXTptr context, bool* fixed0, bool* fixed1, double* x, double* w) {
	thread_local static double Lp_rhs[MAX_N*MAX_N];
	thread_local static int Lp_rmatind[MAX_N*MAX_N*3];
	thread_local static double Lp_rmatval[MAX_N*MAX_N*3];
	int nn = 0, beg = 0;
	for (int j=0; j<N; j++) if (fixed1[j]) for (int k=0; k<N; k++) {
		if (fixed0[k] || fixed1[k]) continue;
		double Ljk = calc_local_Ljk(j, k, fixed0, fixed1);
		double Lj = calc_local_L(j, fixed0, fixed1);
		if (w[j] + w[k] + (-Ljk+Lj) * x[k] - Lj > -EPS) continue;
		Lp_rhs[nn] = Lj;
		Lp_rmatind[beg  ] = i_w(j);
		Lp_rmatval[beg++] = 1;
		Lp_rmatind[beg  ] = i_w(k);
		Lp_rmatval[beg++] = 1;
		Lp_rmatind[beg  ] = i_x(k);
		Lp_rmatval[beg++] = -Ljk + Lj;
		nn++;
	}
	if (nn >= PARAM_CUTS_MIN) {
		_c(CPXcallbackaddusercuts(context, nn, nn*3, Lp_rhs, Lp_sense, Lp_rmatbeg, Lp_rmatind, Lp_rmatval, Lp_purgeable, Lp_local));
		tot_lp_cuts += nn;
	}
}

static char La_sense[MAX_N];
static int La_purgeable[MAX_N];
static int La_local[MAX_N];

void init_La_data() {
	fill_n(La_sense, N, 'G');
	fill_n(La_purgeable, N, CUT_TYPE);
	fill_n(La_local, N, 1);
}

static atomic_int tot_la_cuts = 0;

double calc_local_F0(bool* fixed0, bool* fixed1) {
	// calculate all C_i for this k
	thread_local static double C[MAX_N];
	fill_n(C, N, 0);
	for (int i=0; i<N; i++) {
		C[i] = 0;
		for (int j=0; j<N; j++) {
			if (fixed1[j]) C[i] += B[i][j];
		}
	}
	thread_local static int C_ord[MAX_N];
	for (int i=0; i<N; i++) C_ord[i] = i;
	sort(C_ord, C_ord+N, [&](int ia, int ib){ return C[ia] < C[ib]; });
	// calculate F_k
	int take = M;
	double tot = 0;
	// take all the fixed variables
	for (int i = 0; i < N; i++) {
		if (fixed1[i]) {
			tot += C[i];
			take--;
		}
	}
	assert(take >= 0);
	// take the smallest remaining values
	for (int ord = 0; take && ord < N; ord++) {
		int i = C_ord[ord];
		if (fixed1[i]) continue; // already counted
		// we could do slightly better here be calculating a F0 for each k, skipping here when i==k
		if (fixed0[i]) continue; // cannot be used
		tot += C[i];
		take--;
	}
	assert(take == 0);
	return tot;
}

double calc_local_F1(int k, bool* fixed0, bool* fixed1) {
	// calculate all B_i for this k
	thread_local static double C[MAX_N];
	fill_n(C, N, 0);
	for (int i=0; i<N; i++) {
		C[i] = B[i][k];
		for (int j=0; j<N; j++) {
			if (fixed1[j]) C[i] += B[i][j];
		}
	}
	thread_local static int C_ord[MAX_N];
	for (int i=0; i<N; i++) C_ord[i] = i;
	sort(C_ord, C_ord+N, [&](int ia, int ib){ return C[ia] < C[ib]; });
	// calculate F_k
	int take = M;
	double tot = 0;
	// take all the fixed variables
	for (int i = 0; i < N; i++) {
		if (i == k || fixed1[i]) {
			tot += C[i];
			take--;
		}
	}
	assert(take >= 0);
	// take the smallest remaining values
	for (int ord = 0; take && ord < N; ord++) {
		int i = C_ord[ord];
		if (i == k || fixed1[i]) continue; // already counted
		if (fixed0[i]) continue; // cannot be used
		tot += C[i];
		take--;
	}
	assert(take == 0);
	return tot;
}

void improve_local_L_all(CPXCALLBACKCONTEXTptr context, bool* fixed0, bool* fixed1, double* x, double* w) {
	// use the same value for all `k`s as they would not differ by much anyways
	double f0 = calc_local_F0(fixed0, fixed1);

	thread_local static double La_rhs[MAX_N];
	thread_local static int La_rmatbeg[MAX_N];
	thread_local static int La_rmatind[MAX_N*MAX_N];
	thread_local static double La_rmatval[MAX_N*MAX_N];
	int nn = 0, beg = 0;
	for (int k=0; k<N; k++) {
		if (fixed0[k] || fixed1[k]) continue;
		// check if this cut is violated
		double f1k = calc_local_F1(k, fixed0, fixed1);
		double delta = w[k] + (-f1k+f0) * x[k] - f0;
		for (int j=0; j<N; j++) if (fixed1[j]) {
			delta += w[j];
		}
		if (delta > -EPS) continue;

		La_rhs[nn] = f0;
		La_rmatbeg[nn] = beg;
		for (int j=0; j<N; j++) if (fixed1[j]) {
			La_rmatind[beg  ] = i_w(j);
			La_rmatval[beg++] = 1;
		}
		La_rmatind[beg  ] = i_w(k);
		La_rmatval[beg++] = 1;
		La_rmatind[beg  ] = i_x(k);
		La_rmatval[beg++] = -f1k + f0;
		nn++;
	}
	if (nn >= PARAM_CUTS_MIN) {
		_c(CPXcallbackaddusercuts(context, nn, beg, La_rhs, La_sense, La_rmatbeg, La_rmatind, La_rmatval, La_purgeable, La_local));
		tot_la_cuts += nn;
	}
}

static char M_sense[MAX_N];
static int M_rmatbeg[MAX_N];
static int M_purgeable[MAX_N];
static int M_local[MAX_N];

void init_M_data() {
	fill_n(M_sense, N, 'L');
	fill_n(M_rmatbeg, N, 0);
	fill_n(M_purgeable, N, CUT_TYPE);
	fill_n(M_local, N, 1);
    for (int i = 0; i < N; i++) {
		M_rmatbeg[i] = i * (N+1);
	}
}

static atomic_int tot_m_cuts = 0;

double calc_local_M(int i, bool* fixed0, bool* fixed1) {
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
	return tot;
}

void improve_local_M(CPXCALLBACKCONTEXTptr context, bool* fixed0, bool* fixed1, double* x, double* w) {
	// calculate all local M_i
	thread_local static double lm[MAX_N];
	for (int i=0; i<N; i++) {
		if (fixed0[i] || fixed1[i]) continue;
		lm[i] = calc_local_M(i, fixed0, fixed1);
	}

	// generate new constraints
	thread_local static double M_rhs[MAX_N];
	thread_local static int M_rmatind[MAX_N*(MAX_N+1)];
	thread_local static double M_rmatval[MAX_N*(MAX_N+1)];
	int nn = 0;
    for (int i = 0; i < N; i++) {
		if (fixed0[i] || fixed1[i]) continue;
		double delta = -w[i] - lm[i];
		for (int j = 0; j < N; j++)
			delta += x[j] * (B[i][j] + (j==i ? lm[i] : 0));
		if (delta < EPS) continue;

		int beg = nn * (N+1);
		M_rhs[nn] = lm[i];
		for (int j = 0; j < N; j++) {
			M_rmatind[beg + j] = i_x(j);
			M_rmatval[beg + j] = B[i][j] + (i==j ? lm[i] : 0);
		}
		M_rmatind[beg + N] = i_w(i);
		M_rmatval[beg + N] = -1;
		nn++;
	}
	if (nn >= PARAM_CUTS_MIN) {
		_c(CPXcallbackaddusercuts(context, nn, nn*(N+1), M_rhs, M_sense, M_rmatbeg, M_rmatind, M_rmatval, M_purgeable, M_local));
		tot_m_cuts += nn;
	}
}

void find_fixings(CPXCALLBACKCONTEXTptr context, bool* fixed0, bool* fixed1) {
	// lower/upperbound of the x variables at the local node
	thread_local static double lb[MAX_N], ub[MAX_N];
	_c(CPXcallbackgetlocallb((CPXCALLBACKCONTEXTptr)context, lb, i_x(0), i_x(N-1)));
	_c(CPXcallbackgetlocalub((CPXCALLBACKCONTEXTptr)context, ub, i_x(0), i_x(N-1)));
	for (int i=0; i<N; i++) fixed1[i] = lb[i] > .5;
	for (int i=0; i<N; i++) fixed0[i] = ub[i] < .5;
}

void cuts_generator(CPXCALLBACKCONTEXTptr context) {
	if (PARAM_LIVE_SOL) live_display(context);

	// determine if we have already seen this node (best effort: only thread local)
	thread_local static long long prev_id = -1;
	long long node_id;
	_c(CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEUID, &node_id));
	bool old_node = node_id == prev_id;
	prev_id = node_id;

	//double obj;
	//_c(CPXcallbackgetrelaxationpoint(context, NULL, 0, -1, &obj));
	//cout << node_id << " \t-> " << tot_m_cuts << "/" << tot_lp_cuts << "/" << tot_la_cuts << " \t -> " << obj << "  \t/ " << old_node << very_old_node << endl;

	// producing cuts more then once for the same node could be wasteful
	if (PARAM_CUT_ONCE && old_node) return;

	// keeps track of which variables are fixed to 0/1
	thread_local static bool fixed0[MAX_N], fixed1[MAX_N];
	find_fixings(context, fixed0, fixed1);

	// solution at the current node
	thread_local static double x[MAX_N], w[MAX_N];
	_c(CPXcallbackgetrelaxationpoint(context, x, i_x(0), i_x(N-1), NULL));
	_c(CPXcallbackgetrelaxationpoint(context, w, i_w(0), i_w(N-1), NULL));

	if (PARAM_LOCAL_L) // these cuts are always the same
		improve_local_L(context, fixed0, fixed1, x, w);
	if (PARAM_LOCAL_M)
		improve_local_M(context, fixed0, fixed1, x, w);
	if (PARAM_LOCAL_L_PAIRS)
		improve_local_L_pairs(context, fixed0, fixed1, x, w);
	if (PARAM_LOCAL_L_ALL)
		improve_local_L_all(context, fixed0, fixed1, x, w);
}

double tot_callback_time = 0;

int cplex_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle) {
	CPXENVptr env = (CPXENVptr) userhandle;
	double ts, te;
	CPXgettime(env, &ts);
	cuts_generator(context);
	CPXgettime(env, &te);
	tot_callback_time += te - ts;
	return 0;
}

int main() {
	read_inputs();

	preorder_B();
	if (PARAM_LOCAL_L_PAIRS) preorder_B_pairs();

	if (PARAM_LOCAL_L) init_L_data();
	if (PARAM_LOCAL_L_PAIRS) init_Lp_data();
	if (PARAM_LOCAL_L_ALL) init_La_data();
	if (PARAM_LOCAL_M) init_M_data();

	int status;
    CPXENVptr env = CPXopenCPLEX(&status);
	if (status) abort();
    CPXLPptr lp = CPXcreateprob(env, &status, "qap");
	if (status) abort();

	_c(CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON)); // enable solver logging
	if (PARAM_SINGLE_THREAD)
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
	_c(CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_RELAXATION, cplex_callback, env));

    // solve as mip
	_c(CPXsetintparam(env, CPXPARAM_Preprocessing_Symmetry, 5));
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
	int solstat, nodecnt;
	vector<double> sol(N + N);
	_c(CPXsolution(env, lp, &solstat, &objval, sol.data(), NULL, NULL, NULL));
	_c(CPXgetbestobjval(env, lp, &best_boud));
	nodecnt = CPXgetnodecnt(env, lp);

	cout << "time in user callback = " << tot_callback_time << endl;
	cout << "status = " << solstat << endl;
	cout << "solution = " << objval << endl;
	cout << "applied lp / la / m cuts = " << tot_lp_cuts << " / " << tot_la_cuts << " / " << tot_m_cuts << endl;
	cerr << N << endl;
	for (int i=0; i<N; i++) {
		cerr << int(round(sol[i_x(i)])) << endl;
	}
	cerr << "# obj = " << objval << endl;
	if (solstat != CPXMIP_OPTIMAL)
		cerr << "# bound = " << best_boud << endl;
	cerr << "# time = " << (t_end - t_start) << ", in cb = " << tot_callback_time << endl;
	cerr << "# nodes = " << nodecnt << endl;
	cerr << "# cb"
		<< ", dynamic search = " << PARAM_DYNAMIC_SEARCH
		<< ", local L = " << PARAM_LOCAL_L
		<< ", local L pairs = " << PARAM_LOCAL_L_PAIRS
		<< ", local L all = " << PARAM_LOCAL_L_ALL
		<< ", local M = " << PARAM_LOCAL_M
		<< ", cut once = " << PARAM_CUT_ONCE
		<< ", single thread = " << PARAM_SINGLE_THREAD
		<< ", tl = " << PARAM_TIME_LIMIT << endl;
	cerr << "# n1 = " << M << endl;

    // clean up
    _c(CPXfreeprob(env, &lp));
    _c(CPXcloseCPLEX(&env));

    return 0;
} 
