#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <ilcplex/cplex.h>

#include <inputs.hh>

#define _c(what) if(what) abort()

// seems faster on bigger problems
constexpr bool PARAM_DYNAMIC_SEARCH = true;
// -1 to disable
constexpr double PARAM_TIME_LIMIT = 45;

using namespace std;

int i_x(int i) {
	assert(0 <= i && i < N);
	return i;
}

int i_y(int i, int j) {
	assert(0 <= i && i < N);
	assert(0 <= j && j < N);
	return N + i * N + j;
}

void add_variables_x(CPXENVptr env, CPXLPptr lp) {
	vector<double> lb(N, 0.0);
	vector<double> ub(N, 1.0);
	vector<char> type(N, 'I');
	vector<string> name(N);
	vector<char*> c_name(N);
	for (int i=0; i<N; i++) {
		name[i] = "x_" + to_string(i);
		c_name[i] = const_cast<char*>(name[i].c_str());
	}
	_c(CPXnewcols(env, lp, N, NULL, lb.data(), ub.data(), type.data(), c_name.data()));
}

void add_variables_y(CPXENVptr env, CPXLPptr lp) {
	vector<double> lb(N*N, 0.0);
	vector<double> ub(N*N, 1.0);
	vector<char> type(N*N, 'I');
	vector<double> ycost(N*N);
	vector<string> name(N*N);
	vector<char*> c_name(N*N);
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++) {
		int k = i_y(i, j) - i_y(0,0);
		ycost[k] = B[i][j];
		name[k] = "y_" + to_string(i) + "_" + to_string(j);
		c_name[k] = const_cast<char*>(name[k].c_str());
	}
	_c(CPXnewcols(env, lp, N*N, ycost.data(), lb.data(), ub.data(), type.data(), c_name.data()));
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

// y_ij <= x_i <=> y_ij - x_i <= 0
void add_cons_and1(CPXENVptr env, CPXLPptr lp) {
	vector<double> rhs(N*N, 0);
	vector<char> sense(N*N, 'L');
	vector<int> rmatbeg(N*N, 0);
	vector<int> rmatind(N*N*2, 0);
	vector<double> rmatval(N*N*2, 0);
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
		int k = i*N + j;
		rmatbeg[k] = 2*k;
		rmatind[2*k+0] = i_y(i,j);
		rmatval[2*k+0] = 1;
		rmatind[2*k+1] = i_x(i);
		rmatval[2*k+1] = -1;
	}
    _c(CPXaddrows(env, lp, 0, N*N, 2*N*N, rhs.data(), sense.data(), rmatbeg.data(), rmatind.data(), rmatval.data(), NULL, NULL));
}

// y_ij <= x_j <=> y_ij - x_j <= 0
void add_cons_and2(CPXENVptr env, CPXLPptr lp) {
	vector<double> rhs(N*N, 0);
	vector<char> sense(N*N, 'L');
	vector<int> rmatbeg(N*N, 0);
	vector<int> rmatind(N*N*2, 0);
	vector<double> rmatval(N*N*2, 0);
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
		int k = i*N + j;
		rmatbeg[k] = 2*k;
		rmatind[2*k+0] = i_y(i,j);
		rmatval[2*k+0] = 1;
		rmatind[2*k+1] = i_x(j);
		rmatval[2*k+1] = -1;
	}
    _c(CPXaddrows(env, lp, 0, N*N, 2*N*N, rhs.data(), sense.data(), rmatbeg.data(), rmatind.data(), rmatval.data(), NULL, NULL));
}

// y_ij >= x_i + x_j - 1 <=> y_ij - x_i - x_j >= -1
void add_cons_and3(CPXENVptr env, CPXLPptr lp) {
	vector<double> rhs(N*N, -1);
	vector<char> sense(N*N, 'G');
	vector<int> rmatbeg(N*N, 0);
	vector<int> rmatind(N*N*3, 0);
	vector<double> rmatval(N*N*3, 0);
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
		int k = i*N + j;
		rmatbeg[k] = 3*k;
		rmatind[3*k+0] = i_y(i,j);
		rmatval[3*k+0] = 1;
		rmatind[3*k+1] = i_x(i);
		rmatval[3*k+1] = -1;
		rmatind[3*k+2] = i_x(j);
		rmatval[3*k+2] = -1;
	}
    _c(CPXaddrows(env, lp, 0, N*N, 3*N*N, rhs.data(), sense.data(), rmatbeg.data(), rmatind.data(), rmatval.data(), NULL, NULL));
}

int main() {
	read_inputs();

	int status;
    CPXENVptr env = CPXopenCPLEX(&status);
	if (status) abort();
    CPXLPptr lp = CPXcreateprob(env, &status, "qap");
	if (status) abort();

	// enable solver logging
	_c(CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON));
	_c(CPXsetintparam(env, CPXPARAM_Threads, 1));

    _c(CPXchgobjsen(env, lp, CPX_MIN));

	add_variables_x(env, lp);
	add_variables_y(env, lp);

	add_cons_M(env, lp);
	add_cons_and1(env, lp);
	add_cons_and2(env, lp);
	add_cons_and3(env, lp);

	// write problem
	//_c(CPXwriteprob(env, lp, "p.lp", NULL));

    // solve as mip
	if (!PARAM_DYNAMIC_SEARCH)
		_c(CPXsetintparam(env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL));
	if (PARAM_TIME_LIMIT != -1)
		_c(CPXsetdblparam(env, CPXPARAM_TimeLimit, PARAM_TIME_LIMIT));
	double t_start, t_end;
	_c(CPXgettime(env, &t_start));
    _c(CPXmipopt(env, lp));
	_c(CPXgettime(env, &t_end));

    // print solution
    double objval;
	int solstat;
	vector<double> sol(N*N + N);
	_c(CPXsolution(env, lp, &solstat, &objval, sol.data(), NULL, NULL, NULL));

	cout << "solution = " << objval << endl;
	cerr << N << " " << M << endl;
	for (int i=0; i<N; i++) {
		cerr << i << " " << int(round(sol[i_x(i)])) << endl;
	}
	cerr << "# obj = " << objval << endl;
	cerr << "# time = " << (t_end - t_start) << endl;
	cerr << "# toff, no dynamic search" << endl;

    // clean up
    _c(CPXfreeprob(env, &lp));
    _c(CPXcloseCPLEX(&env));

    return 0;
} 
