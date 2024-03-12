#include <cstdlib>
#include <live.hh>
#include <inputs.hh>
#include <iostream>
#include <fstream>

using namespace std;

#define _c(what) if (int _error = what) { \
	cout << "CPX error: " #what << endl; cout << "CPX error: " << _error << endl; abort(); }

void live_display(CPXCALLBACKCONTEXTptr context) {
	if (rand()%100) return;

	long long node_id;
	_c(CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEUID, &node_id));

	static double x[MAX_N], w[MAX_N];
	double obj;
	_c(CPXcallbackgetrelaxationpoint(context, x, 0, N-1, &obj));
	_c(CPXcallbackgetrelaxationpoint(context, w, N, 2*N-1, NULL));

	static double ub[MAX_N], lb[MAX_N];
	_c(CPXcallbackgetlocallb((CPXCALLBACKCONTEXTptr)context, lb, 0, N));
	_c(CPXcallbackgetlocalub((CPXCALLBACKCONTEXTptr)context, ub, 0, N));

	ofstream of("/tmp/sol");
	of << N << endl;
	for (int i=0; i<N; i++) {
		of << (ub[i]<.5 ? 0 : lb[i]>.5 ? 1 : .5) << endl;
	}
	for (int i=0; i<N; i++) {
		of << x[i] << endl;
	}
	for (int i=0; i<N; i++) {
		of << w[i] << endl;
	}

	cout << "node " << hex << node_id << dec << ", obj = " << obj << endl;
	system("python ../../plot/live.py /tmp/sol");
}
