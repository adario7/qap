#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <ilcplex/cplex.h>

#include <inputs.hh>
#include <live.hh>
#include <params.hh>

// static paramteres
constexpr double EPS = 10; // add cuts that are violated by this much
constexpr int CUT_TYPE = CPX_USECUT_FILTER;

#define _c(what) if (int _error = what) { \
	cout << "CPX error " << _error << ": " #what << endl; abort(); }

#define _w(what) if (int _error = what) { \
	cout << "Warn: CPX error " << _error << ": " #what << endl;  }

using namespace std;

int i_x(int i) {
	assert(0 <= i && i < N);
	return i;
}

void print_results(CPXENVptr env, CPXLPptr lp, double time) {
	// get the solutions
	double objval, best_boud;
	int solstat, nodecnt;
	vector<double> sol(N + N);
	_c(CPXsolution(env, lp, &solstat, &objval, sol.data(), NULL, NULL, NULL));
	_c(CPXgetbestobjval(env, lp, &best_boud));
	nodecnt = CPXgetnodecnt(env, lp);

	// print statistics and solutions
	cout << "status = " << solstat << endl;
	cout << "solution = " << objval << endl;
	cerr << N << endl;
	for (int i=0; i<N; i++) {
		cerr << int(round(sol[i_x(i)])) << endl;
	}
	cerr << "# status = " << solstat << endl;
	cerr << "# obj = " << objval << endl;
	if (solstat != CPXMIP_OPTIMAL && solstat != CPXMIP_OPTIMAL_TOL) {
		cerr << "# bound = " << best_boud << endl;
		cerr << "# gap = " << 100*(objval - best_boud)/objval << endl;
	}
	cerr << "# nodes = " << nodecnt << endl
		<< "# cb"
		<< ", seed = " << PARAM_SEED
		<< ", dynamic search = " << PARAM_DYNAMIC_SEARCH
		<< ", local L = " << PARAM_LOCAL_L
		<< ", local L pairs = " << PARAM_LOCAL_L_PAIRS
		<< ", local L all = " << PARAM_LOCAL_L_ALL
		<< ", local M = " << PARAM_LOCAL_M
		<< ", local F = " << PARAM_LOCAL_FREE
		<< ", cut once = " << PARAM_CUT_ONCE
		<< ", rel d = " << PARAM_REL_DELTA
		<< ", exp later = " << PARAM_EXP_LATER
		<< ", min cuts = " << PARAM_CUTS_MIN
		<< ", cplex cuts = " << PARAM_CPLEX_CUTS
		<< ", branching = " << PARAM_BRANCHING
		<< ", single thread = " << PARAM_SINGLE_THREAD
		<< ", opportunistic = " << PARAM_OPPORTUNISTIC
		<< ", eps = " << EPS
		<< ", cut type = " << CUT_TYPE
		<< ", tl = " << PARAM_TIME_LIMIT
		<< ", ml = " << PARAM_MEMLIMIT
		<< endl;
	cerr << "# n1 = " << M << endl;
}

int main(int argc, char** argv) {
	read_parameters(argc, argv);
	read_inputs();

    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    int status;

    env = CPXopenCPLEX(&status);
    if (env == NULL) {
        char errmsg[CPXMESSAGEBUFSIZE];
        CPXgeterrorstring(env, status, errmsg);
        fprintf(stderr, "%s", errmsg);
		abort();
    }

    lp = CPXcreateprob(env, &status, "quadratic_mip");
    if (lp == NULL) {
        fprintf(stderr, "Failed to create LP.\n");
		abort();
    }

    int num_vars = N;
    double obj[num_vars];
    char coltype[num_vars];
    for (int i = 0; i < num_vars; i++) {
        obj[i] = 0.0;
        coltype[i] = 'B'; // Binary variables
    }

    _c(CPXnewcols(env, lp, num_vars, obj, NULL, NULL, coltype, NULL));

    int matbeg[1] = {0};
    int matind[num_vars];
    double matval[num_vars];
    for (int i = 0; i < num_vars; i++) {
        matind[i] = i;
        matval[i] = 1.0;
    }
    double rhs[1] = {(double)M};
    char sense[1] = {'E'};

    _c(CPXaddrows(env, lp, 0, 1, num_vars, rhs, sense, matbeg, matind, matval, NULL, NULL));

    int qmatbeg[num_vars];
    int qmatcnt[num_vars];
    int qmatind[num_vars * num_vars];
    double qmatval[num_vars * num_vars];
    int k = 0;

    for (int i = 0; i < N; i++) {
        qmatbeg[i] = k;
        qmatcnt[i] = N;
        for (int j = 0; j < N; j++) {
            qmatind[k] = j;
            qmatval[k] = 2*B[i][j];
            k++;
        }
    }
    _c(CPXcopyquad(env, lp, qmatbeg, qmatcnt, qmatind, qmatval));

	apply_parameters(env, lp);
    _c(CPXmipopt(env, lp));

	print_results(env, lp, 0);

    if (lp != NULL) CPXfreeprob(env, &lp);
    if (env != NULL) CPXcloseCPLEX(&env);

    return 0;
}
