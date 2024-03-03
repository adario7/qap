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

#define _c(what) if (int _error = what) { \
	cout << "CPX error: " #what << endl; cout << "CPX error: " << _error << endl; abort(); }

using namespace std;

int B_order[MAX_N][MAX_N];

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

static double C_rhs[MAX_N];
static char C_sense[MAX_N];
static int C_rmatbeg[MAX_N];
static int C_rmatind[MAX_N*2];
thread_local static double C_rmatval[MAX_N*2];
static int C_purgeable[MAX_N];
static int C_local[MAX_N];

void init_cuts_data() {
	fill_n(C_rhs, N, 0);
	fill_n(C_sense, N, 'G');
	fill_n(C_rmatbeg, N, 0);
	fill_n(C_purgeable, N, CPX_USECUT_FILTER);
	fill_n(C_local, N, 1);
    for (int i = 0; i < N; i++) {
		int beg = i * 2;
		C_rmatbeg[i] = beg;
		C_rmatind[beg + 0] = i_w(i);
		C_rmatind[beg + 1] = i_x(i);
	}
}

//int cuts_generator(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p) {
//int cuts_generator(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *nodeindex_p, int *useraction_p) {
int cuts_generator(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle) {
	// upperbound of the x variables at the local node
	thread_local static double lb[MAX_N];
	_c(CPXcallbackgetlocallb((CPXCALLBACKCONTEXTptr)context, lb, i_x(0), i_x(N-1)));
	//thread_local static double ub[MAX_N];
	//_c(CPXcallbackgetlocalub((CPXCALLBACKCONTEXTptr)context, ub, i_x(0), i_x(N-1)));

	// keeps track of which variables are fixed to 1
	thread_local static bool fixed[MAX_N];
	for (int i=0; i<N; i++) fixed[i] = lb[i] > .5;
	int num_fixed = count(fixed, fixed+N, true);
	assert(num_fixed <= M); // should never fix more than M varaibles to 1

	// calculate all L_i
	thread_local static double L[MAX_N];
	for (int i=0; i<N; i++) {
		int take = M;
		double tot = 0;
		// always take i
		take--;
		tot += B[i][i];
		// take all the fixed variables
		for (int j = 0; take && j < N; j++) {
			if (j == i) continue; // already counted
			if (fixed[j]) {
				tot += B[i][j];
				take--;
			}
		}
		// take the smallest remaining values
		for (int k = 0; take && k < N; k++) {
			int j = B_order[i][k];
			if (j == i) continue; // already counted
			if (fixed[j]) continue; // already counted
			tot += B[i][j];
			take--;
		}
		L[i] = tot;
	}

	// generate new contraints
    for (int i = 0; i < N; i++) {
		int beg = i * 2;
		C_rmatval[beg + 0] = 1;
		C_rmatval[beg + 1] = -L[i];
	}
    _c(CPXcallbackaddusercuts(context, N, N*2, C_rhs, C_sense, C_rmatbeg, C_rmatind, C_rmatval, C_purgeable, C_local));
	
	return 0;
}

int main() {
	read_inputs();
	preorder_B();
	init_cuts_data();

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
	//add_cons_L(env, lp); // not needed as they are added locally

	// write problem
	//_c(CPXwriteprob(env, lp, "p.lp", NULL));

	// add a callback to generate cuts
	_c(CPXsetintparam(env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF)); // we generate cuts for the original problem, not the presolved one
	//_c(CPXsetusercutcallbackfunc(env, cuts_generator, NULL));
	//_c(CPXsetnodecallbackfunc(env, cuts_generator, NULL));
	_c(CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_RELAXATION, cuts_generator, NULL));


    // solve as mip
	_c(CPXsetintparam(env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL)); // TODO dynamic search vs traditional search
	double t_start, t_end;
	_c(CPXgettime(env, &t_start));
    _c(CPXmipopt(env, lp));
	_c(CPXgettime(env, &t_end));

    // print solution
    double objval;
	int solstat;
	vector<double> sol(N + N);
	_c(CPXsolution(env, lp, &solstat, &objval, sol.data(), NULL, NULL, NULL));

	cout << "solution = " << objval << endl;
	cerr << N << " " << M << endl;
	for (int i=0; i<N; i++) {
		cout << "  x_" << i << " = " << int(sol[i_x(i)]) << endl;
		cerr << i << " " << int(sol[i_x(i)]) << endl;
	}
	cerr << "# obj = " << objval << endl;
	cerr << "# time = " << (t_end - t_start) << endl;
	cerr << "# cb, no dynamic search, static alloc" << endl;

    // clean up
    _c(CPXfreeprob(env, &lp));
    _c(CPXcloseCPLEX(&env));

    return 0;
} 
