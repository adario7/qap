#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <ilcplex/cplex.h>
#include <vector>

#define _c(what) status =  what; \
	if (status) abort();

using namespace std;

const int MAX_N = 256;

uint N;
double C;
double V[MAX_N];
double W[MAX_N];

void read_inputs() {
	cin >> N; // T
	cin >> N >> C;
	C += 0.5;
	assert(N <= MAX_N);
	for (int i = 0; i < N; i++)
		cin >> V[i];
	for (int i = 0; i < N; i++)
		cin >> W[i];
}

int main() {
	read_inputs();

	int status;
    CPXENVptr env = CPXopenCPLEX(&status);
	if (status) abort();
    CPXLPptr lp = CPXcreateprob(env, &status, "knapsack");
	if (status) abort();

	// enable solver logging
	_c(CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON));

    _c(CPXchgobjsen(env, lp, CPX_MAX));

    // variables
	vector<double> lb(N, 0.0);
	vector<double> ub(N, 1.0);
	vector<char> xctype(N, 'I');
    _c(CPXnewcols(env, lp, N, V, lb.data(), ub.data(), xctype.data(), NULL));

    // capacity constraint
    double rhs[1] = {C};
    int rmatbeg[1] = {0};
	vector<int> rmatind(N);
    for (int i = 0; i < N; i++) {
        rmatind[i] = i;
    }
    _c(CPXaddrows(env, lp, 0, 1, N, rhs, "L", rmatbeg, rmatind.data(), W, NULL, NULL));

    // solve as mip
    _c(CPXmipopt(env, lp));

    // print solution
    double objval;
	int solstat;
	static double x[MAX_N];
	_c(CPXsolution(env, lp, &solstat, &objval, x, NULL, NULL, NULL));

    printf("Solution: %f\n", objval);
    printf("Taking:");
	for (int i=0; i<N; i++) {
		if (x[i] > 0.5)
			printf(" %i", i);
	}
	printf("\n");

    // clean up
    _c(CPXfreeprob(env, &lp));
    _c(CPXcloseCPLEX(&env));

    return 0;
} 
