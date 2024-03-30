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
#include <params.hh>

// static paramteres
constexpr double EPS = 1; // add cuts that are violated by this much
constexpr int CUT_TYPE = CPX_USECUT_FORCE;

#define _c(what) if (int _error = what) { \
	cout << "CPX error: " #what << endl; cout << "CPX error: " << _error << endl; abort(); }

using namespace std;

constexpr int MAX_RC = MAX_N + MAX_N + MAX_N + MAX_N;
constexpr int MAX_NZC = MAX_N*2 + MAX_N*(MAX_N+1) + MAX_N*3 + MAX_N*(MAX_N+1);

struct relpoint_t {
	double x[MAX_N], w[MAX_N];
	double lb[MAX_N], ub[MAX_N];
	bool fixed0[MAX_N], fixed1[MAX_N];
};

struct cutbuf_t {
	int rc = 0, nzc = 0;
	int n_l = 0, n_m = 0, n_p = 0, n_a = 0;
	double rhs[MAX_RC];
	char sense[MAX_RC];
	int rmatbeg[MAX_RC];
	int rmatind[MAX_NZC];
	double rmatval[MAX_NZC];
};

// proerdered values
int B_order[MAX_N][MAX_N];
int B_jk_order[MAX_N][MAX_N][MAX_N];
// global buffers
int gb_purgeable[MAX_RC];
int gb_local[MAX_RC];
// total cuts statistics
atomic_int tot_l_cuts = 0;
atomic_int tot_p_cuts = 0;
atomic_int tot_a_cuts = 0;
atomic_int tot_m_cuts = 0;
double tot_callback_time = 0;

int i_x(int i) {
	assert(0 <= i && i < N);
	return i;
}

int i_w(int i) {
	assert(0 <= i && i < N);
	return N + i;
}

void cut_begin(cutbuf_t& buf, double rhs, char sense) {
	buf.rhs[buf.rc] = rhs;
	buf.sense[buf.rc] = sense;
	buf.rmatbeg[buf.rc] = buf.nzc;
	buf.rc++;
}

void cut_add(cutbuf_t& buf, int rmatind, double rmatval) {
	buf.rmatind[buf.nzc] = rmatind;
	buf.rmatval[buf.nzc] = rmatval;
	buf.nzc++;
}

void init_global_buffers() {
	fill_n(gb_purgeable, MAX_RC, CUT_TYPE);
	fill_n(gb_local, MAX_RC, 1);
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

double calc_local_L(int i, const bool* fixed0, const bool* fixed1) {
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

void calc_L_cuts(const relpoint_t& rp, cutbuf_t& buf) {
	thread_local static double L_rmatval[MAX_N*2];
	thread_local static int L_rmatind[MAX_N*2];
	for (int i = 0; i < N; i++) {
		if (rp.fixed0[i] || rp.fixed1[i]) continue;
		double ll = calc_local_L(i, rp.fixed0, rp.fixed1);
		double delta = rp.w[i] - ll * rp.x[i];
		if (delta > -EPS) continue;

		cut_begin(buf, 0, 'G');
		cut_add(buf, i_w(i), 1);
		cut_add(buf, i_x(i), -ll);
		buf.n_l++;
	}
}

double calc_local_M(int i, const bool* fixed0, const bool* fixed1) {
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

void calc_M_cuts(const relpoint_t& rp, cutbuf_t& buf) {
	for (int i = 0; i < N; i++) {
		if (rp.fixed0[i] || rp.fixed1[i]) continue;
		double lm = calc_local_M(i, rp.fixed0, rp.fixed1);
		double delta = -rp.w[i] - lm;
		for (int j = 0; j < N; j++)
			delta += rp.x[j] * (B[i][j] + (j==i ? lm : 0));
		if (delta < EPS) continue;

		cut_begin(buf, lm, 'L');
		for (int j = 0; j < N; j++) {
			cut_add(buf, i_x(j), B[i][j] + (i==j ? lm : 0));
		}
		cut_add(buf, i_w(i), -1);
		buf.n_m++;
	}
}

double calc_local_Ljk(int j, int k, const bool* fixed0, const bool* fixed1) {
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

void calc_P_cuts(const relpoint_t& rp, cutbuf_t& buf) {
	double Ljks[MAX_N], Ljs[MAX_N];
	for (int j=0; j<N; j++) if (rp.fixed1[j]) {
		// we could do slightly better with different Ljs for each k where we fix k to 0
		Ljs[j] = calc_local_L(j, rp.fixed0, rp.fixed1);
	}
	for (int k=0; k<N; k++) {
		if (rp.fixed0[k] || rp.fixed1[k]) continue;
		int best_j = -1;
		double best_delta = -EPS;
		for (int j=0; j<N; j++) if (rp.fixed1[j]) {
			double ljk = Ljks[j] = calc_local_Ljk(j, k, rp.fixed0, rp.fixed1);
			double lj = Ljs[j];
			double delta = rp.w[j] + rp.w[k] + (-ljk+lj) * rp.x[k] - lj;
			if (delta < best_delta) {
				best_delta = delta;
				best_j = j;
			}
		}
		if (best_j != -1) {
			int j = best_j;
			double ljk = Ljks[j], lj = Ljs[j];
			cut_begin(buf, lj, 'G');
			cut_add(buf, i_w(j), 1);
			cut_add(buf, i_w(k), 1);
			cut_add(buf, i_x(k), -ljk + lj);
			buf.n_p++;
		}
	}
}

double calc_local_F0(const bool* fixed0, const bool* fixed1) {
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

double calc_local_F1(int k, const bool* fixed0, const bool* fixed1) {
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

void calc_A_cuts(const relpoint_t& rp, cutbuf_t& buf) {
	// use the same value for all `k`s as they would not differ by much anyways
	double f0 = calc_local_F0(rp.fixed0, rp.fixed1);

	for (int k=0; k<N; k++) {
		if (rp.fixed0[k] || rp.fixed1[k]) continue;
		// check if this cut is violated
		double f1k = calc_local_F1(k, rp.fixed0, rp.fixed1);
		double delta = rp.w[k] + (-f1k+f0) * rp.x[k] - f0;
		for (int j=0; j<N; j++) if (rp.fixed1[j]) {
			delta += rp.w[j];
		}
		if (delta > -EPS) continue;

		cut_begin(buf, f0, 'G');
		for (int j=0; j<N; j++) if (rp.fixed1[j]) {
			cut_add(buf, i_w(j), 1);
		}
		cut_add(buf, i_w(k), 1);
		cut_add(buf, i_x(k), -f1k + f0);
		buf.n_a++;
	}
}

void find_fixings(CPXCALLBACKCONTEXTptr context, relpoint_t& rp) {
	// lower/upperbound of the x variables at the local node
	_c(CPXcallbackgetlocallb(context, rp.lb, i_x(0), i_x(N-1)));
	_c(CPXcallbackgetlocalub(context, rp.ub, i_x(0), i_x(N-1)));
	for (int i=0; i<N; i++) rp.fixed0[i] = rp.ub[i] < .5;
	for (int i=0; i<N; i++) rp.fixed1[i] = rp.lb[i] > .5;
}


void cuts_generator(CPXCALLBACKCONTEXTptr context) {
	if (PARAM_LIVE_SOL) live_display(context);

	// determine if we have already seen this node (best effort: only thread local)
	thread_local static long long prev_id = -1;
	long long node_id;
	_c(CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEUID, &node_id));
	bool old_node = node_id == prev_id;
	prev_id = node_id;

	// producing cuts more then once for the same node could be wasteful
	if (PARAM_CUT_ONCE == 2 && old_node) return;

	// display ram progress
	if (!old_node) {
		long long node_c;
		_c(CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODECOUNT, &node_c));
		if (node_c % 50000 == 0) {
				cout << "applied l / m / lp / la cuts = " << tot_l_cuts << " / " << tot_m_cuts << " / " << tot_p_cuts << " / " << tot_a_cuts << endl;
				system("free -h");
		}
	}

	// local buffers
	relpoint_t rp;
	cutbuf_t buf;
	buf.rc = buf.nzc = buf.n_l = buf.n_m = buf.n_p = buf.n_a = 0;

	// keeps track of which variables are fixed to 0/1
	find_fixings(context, rp);

	// solution at the current node
	_c(CPXcallbackgetrelaxationpoint(context, rp.x, i_x(0), i_x(N-1), NULL));
	_c(CPXcallbackgetrelaxationpoint(context, rp.w, i_w(0), i_w(N-1), NULL));

	// compute enabled cuts
	if (PARAM_LOCAL_L)
		calc_L_cuts(rp, buf);
	if (PARAM_LOCAL_M)
		calc_M_cuts(rp, buf);
	if (PARAM_LOCAL_L_PAIRS && (!old_node || PARAM_CUT_ONCE == 0))
		calc_P_cuts(rp, buf);
	if (PARAM_LOCAL_L_ALL && (!old_node || PARAM_CUT_ONCE == 0))
		calc_A_cuts(rp, buf);

	assert(buf.rc <= MAX_RC);
	assert(buf.nzc <= MAX_NZC);
	if (buf.rc >= PARAM_CUTS_MIN) {
		_c(CPXcallbackaddusercuts(context, buf.rc, buf.nzc, buf.rhs, buf.sense, buf.rmatbeg, buf.rmatind, buf.rmatval, gb_purgeable, gb_local));
		tot_l_cuts += buf.n_l;
		tot_m_cuts += buf.n_m;
		tot_p_cuts += buf.n_p;
		tot_a_cuts += buf.n_a;
	}
}

int cplex_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle) {
	CPXENVptr env = (CPXENVptr) userhandle;
	double ts, te;
	CPXgettime(env, &ts);
	cuts_generator(context);
	CPXgettime(env, &te);
	tot_callback_time += te - ts;
	return 0;
}

int main(int argc, char** argv) {
	read_parameters(argc, argv);
	read_inputs();

	preorder_B();
	if (PARAM_LOCAL_L_PAIRS) preorder_B_pairs();

	init_global_buffers();

	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	if (status) abort();
	CPXLPptr lp = CPXcreateprob(env, &status, "qap");
	if (status) abort();

	_c(CPXchgobjsen(env, lp, CPX_MIN));

	add_variables_x(env, lp);
	add_variables_w(env, lp);

	add_cons_M(env, lp);
	add_cons_xw(env, lp);
	if (!PARAM_LOCAL_L)
		add_cons_L(env, lp);

	// add a callback to generate cuts
	_c(CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_RELAXATION, cplex_callback, env));

	// solve as mip
	apply_parameters(env, lp);
	double t_start, t_end;
	_c(CPXgettime(env, &t_start));
	_c(CPXmipopt(env, lp));
	_c(CPXgettime(env, &t_end));

	// get the solutions
	double objval, best_boud;
	int solstat, nodecnt;
	vector<double> sol(N + N);
	_c(CPXsolution(env, lp, &solstat, &objval, sol.data(), NULL, NULL, NULL));
	_c(CPXgetbestobjval(env, lp, &best_boud));
	nodecnt = CPXgetnodecnt(env, lp);

	// print statistics and solutions
	cout << "time in user callback = " << tot_callback_time << endl;
	cout << "status = " << solstat << endl;
	cout << "solution = " << objval << endl;
	cout << "applied l / m / lp / la cuts = " << tot_l_cuts << " / " << tot_m_cuts << " / " << tot_p_cuts << " / " << tot_a_cuts << endl;
	cerr << N << endl;
	for (int i=0; i<N; i++) {
		cerr << int(round(sol[i_x(i)])) << endl;
	}
	cerr << "# status = " << solstat << endl;
	cerr << "# obj = " << objval << endl;
	if (solstat != CPXMIP_OPTIMAL)
		cerr << "# bound = " << best_boud << endl;
	cerr << "# time = " << (t_end - t_start) << ", in cb = " << tot_callback_time << endl
		<< "# nodes = " << nodecnt << endl
		<< "# l / m / lp / la cuts = " << tot_l_cuts << " / " << tot_m_cuts << " / " << tot_p_cuts << " / " << tot_a_cuts << endl
		<< "# cb"
		<< ", dynamic search = " << PARAM_DYNAMIC_SEARCH
		<< ", local L = " << PARAM_LOCAL_L
		<< ", local L pairs = " << PARAM_LOCAL_L_PAIRS
		<< ", local L all = " << PARAM_LOCAL_L_ALL
		<< ", local M = " << PARAM_LOCAL_M
		<< ", cut once = " << PARAM_CUT_ONCE
		<< ", min cuts = " << PARAM_CUTS_MIN
		<< ", cplex cuts = " << PARAM_CPLEX_CUTS
		<< ", single thread = " << PARAM_SINGLE_THREAD
		<< ", opportunistic = " << PARAM_OPPORTUNISTIC
		<< ", eps = " << EPS
		<< ", cut type = " << CUT_TYPE
		<< ", tl = " << PARAM_TIME_LIMIT
		<< ", ml = " << PARAM_MEMLIMIT
		<< endl;
	cerr << "# n1 = " << M << endl;

	// clean up
	_c(CPXfreeprob(env, &lp));
	_c(CPXcloseCPLEX(&env));

	return 0;
} 
