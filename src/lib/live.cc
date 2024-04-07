#include <atomic>
#include <cstdlib>
#include <live.hh>
#include <inputs.hh>
#include <iostream>
#include <fstream>

using namespace std;

#define _c(what) if (int _error = what) { \
	cout << "CPX error: " #what << endl; cout << "CPX error: " << _error << endl; abort(); }

extern atomic_int tot_l_cuts;
extern atomic_int tot_p_cuts;
extern atomic_int tot_a_cuts;
extern atomic_int tot_m_cuts;
extern atomic_int tot_f_cuts;

void live_display(CPXCALLBACKCONTEXTptr context) {
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

	static long long prev_id = -1;
	static int prev_l = 0, prev_m = 0, prev_p = 0, prev_a = 0, prev_f = 0;

	cout << "node " << hex << node_id << dec << ", obj = " << obj << endl;
	if (node_id == prev_id) {
		cout << " * after " << (tot_l_cuts-prev_l) << " / " << (tot_m_cuts-prev_m) << " / " << (tot_p_cuts-prev_p) << " / " << (tot_a_cuts-prev_a) << " / " << (tot_f_cuts-prev_f) << " cuts" << endl;
	}
	system("python ../../plot/live.py /tmp/sol");

	prev_id = node_id;
	prev_l = tot_l_cuts;
	prev_m = tot_m_cuts;
	prev_p = tot_p_cuts;
	prev_a = tot_a_cuts;
	prev_f = tot_f_cuts;
}
