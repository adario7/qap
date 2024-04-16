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
extern atomic_int tot_m_cuts;
extern atomic_int tot_p_cuts;
extern atomic_int tot_a_cuts;
extern atomic_int tot_f_cuts;

double calc_local_L(int i, const bool* fixed0, const bool* fixed1);
double calc_local_M(int i, const bool* fixed0, const bool* fixed1);

void live_display(CPXCALLBACKCONTEXTptr context) {
	long long node_id, node_cnt;
	double best_sol;
	_c(CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEUID, &node_id));
	_c(CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODECOUNT, &node_cnt));
	_c(CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &best_sol));
	if (node_cnt<100) return;

	static double x[MAX_N], w[MAX_N];
	double obj;
	_c(CPXcallbackgetrelaxationpoint(context, x, 0, N-1, &obj));
	_c(CPXcallbackgetrelaxationpoint(context, w, N, 2*N-1, NULL));

	static double ub[MAX_N], lb[MAX_N];
	static bool fixed0[MAX_N], fixed1[MAX_N];
	_c(CPXcallbackgetlocallb((CPXCALLBACKCONTEXTptr)context, lb, 0, N-1));
	_c(CPXcallbackgetlocalub((CPXCALLBACKCONTEXTptr)context, ub, 0, N-1));
	for (int i=0; i<N; i++) fixed0[i] = ub[i] < .5;
	for (int i=0; i<N; i++) fixed1[i] = lb[i] > .5;


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

	double quadr_obj = 0;
	for (int i=0; i<N; i++) {
		double qw = 0;
		for (int j=0; j<N; j++) {
			qw += B[i][j] * x[i] * x[j];
		}
		of << qw << endl;
		quadr_obj += qw;
	}

	for (int i=0; i<N; i++) {
		double unfeas = 0.5 - abs(x[i] - 0.5);
		of << unfeas*(calc_local_M(i, fixed0, fixed1)-calc_local_L(i, fixed0, fixed1)) << endl;
	}

	for (int i=0; i<N; i++) {
		double unfeas = 0.5 - abs(x[i] - 0.5);
		double gap = -w[i];
		for (int j=0; j<N; j++) gap += x[i] * x[j] * B[i][j];
		of << (unfeas*gap) << endl;
	}

	static long long prev_id = -1;
	static int prev_l = 0, prev_m = 0, prev_p = 0, prev_a = 0, prev_f = 0;

	cout << "node " << node_id << ", \tobj = " << obj << ", \tqadr obj = " << quadr_obj << endl;
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
