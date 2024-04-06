#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <ilcplex/cplex.h>

#include <inputs.hh>

using namespace std;

// parameters
// -1 to disable
constexpr double PARAM_TIME_LIMIT = -1;

#define _c(what) if (int _error = what) { \
	cout << "CPX error: " #what << endl; cout << "CPX error: " << _error << endl; abort(); }

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
	vector<string> name(N);
	vector<char*> c_name(N);
	for (int i=0; i<N; i++) {
		name[i] = "x_" + to_string(i);
		c_name[i] = const_cast<char*>(name[i].c_str());
	}
	_c(CPXnewcols(env, lp, N, NULL, NULL, ub.data(), NULL, c_name.data()));
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

using nbitset = bitset<MAX_N>;

struct node_t {
	nbitset fixed;
	nbitset fvalue;
	double obj;
	int branch;
	double unfeas;
};

double calc_local_L(int i, nbitset fixed0, nbitset fixed1) {
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

double calc_local_M(int i, nbitset fixed0, nbitset fixed1) {
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

int main() {
	read_inputs();

	preorder_B();

	int status;
    CPXENVptr env = CPXopenCPLEX(&status);
	if (status) abort();
    CPXLPptr lp = CPXcreateprob(env, &status, "qap");
	if (status) abort();

	// enable solver logging
	_c(CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_OFF));
	_c(CPXsetintparam(env, CPXPARAM_Threads, 1));

    _c(CPXchgobjsen(env, lp, CPX_MIN));

	add_variables_x(env, lp);
	add_variables_w(env, lp);

	// conss [0 -> N)
	add_cons_xw(env, lp);
	// conss [N -> 2N)
	add_cons_L(env, lp);
	// cons 2N
	add_cons_M(env, lp);

	// write problem
	//_c(CPXwriteprob(env, lp, "p.lp", NULL));

	double t_start, t_end;
	_c(CPXgettime(env, &t_start));

	double incumbent = INFINITY;
	static double inc_sol[MAX_N];

	uint node_n = 0;
	auto cmp_feas = [&](const node_t& a, const node_t& b) {
		return a.unfeas > b.unfeas;
	};
	auto cmp_bound = [&](const node_t& a, const node_t& b) {
		return a.obj > b.obj;
	};
	priority_queue<node_t, vector<node_t>, decltype(cmp_bound)> b_queue(cmp_bound);
	priority_queue<node_t, vector<node_t>, decltype(cmp_feas)> f_queue(cmp_feas);

	map<double, int, less<double>> glob_bound;

	auto curr_bound = [&]() { return glob_bound.empty() ? incumbent : glob_bound.begin()->first; };
	auto print_info = [&]() {
		_c(CPXgettime(env, &t_end));
		double b = curr_bound();
		//double b = b_queue.empty() ? incumbent : b_queue.top().obj;
		double gap = 100 * (incumbent - b) / incumbent;
		cout << int(t_end-t_start) << "s /\t "<< node_n << " /\t " << (b_queue.size()+f_queue.size()) << " /\t " << incumbent << " /\t " << b << " /\t " << gap << "%" << endl;
	};

	auto search_node = [&](nbitset fixed, nbitset fvalue) {
		//cout << "search " << bitset<16>(fixed.to_ulong()) << " / " << bitset<16>(fvalue.to_ulong()) << endl;
		// fix variables bounds
		static int bound_indices[MAX_N];
		static char bound_lu[MAX_N];
		static double bound[MAX_N];
		for (int i=0; i<N; i++) {
			bound_indices[i] = i;
			bound_lu[i] = 'L';
			if (fixed.test(i) && fvalue.test(i)==1) bound[i] = 1;
			else bound[i] = 0;
		}
		_c(CPXchgbds(env, lp, N, bound_indices, bound_lu, bound));
		for (int i=0; i<N; i++) {
			bound_indices[i] = i;
			bound_lu[i] = 'U';
			if (fixed.test(i) && fvalue.test(i)==0) bound[i] = 0;
			else bound[i] = 1;
		}
		_c(CPXchgbds(env, lp, N, bound_indices, bound_lu, bound));
		// local L, M
		nbitset fixed0 = fixed & ~fvalue, fixed1 = fixed & fvalue;
		int rowlist[2 * MAX_N], collist[2 * MAX_N];
		double vallist[2 * MAX_N];
		double localL[MAX_N], localM[MAX_N];
		for (int i=0; i<N; i++) if (!fixed[i]) {
			localL[i] = calc_local_L(i, fixed0, fixed1);
			localM[i] = calc_local_M(i, fixed0, fixed1);
		}
		int nr = 0;
		for (int i=0; i<N; i++) {
			// w_i - L_i x_i >= 0, conss [N -> 2N)
			rowlist[nr] = N + i;
			collist[nr] = i_x(i);
			vallist[nr] = fixed[i] ? 0 : -localL[i];
			nr++;
		}
		for (int i=0; i<N; i++) {
			if (fixed[i]) continue; // an old M on a fixed variable does not hurt
			// sum_{j!=i}( b_ij x_j ) + (M_i + b_ii) x_i - w_i <= M_i, conss [0 -> N)
			rowlist[nr] = i;
			collist[nr] = i_x(i);
			vallist[nr] = localM[i] + B[i][i];
			nr++;
		}
		_c(CPXchgcoeflist(env, lp, nr, rowlist, collist, vallist));
		nr = 0;
		for (int i=0; i<N; i++) {
			if (fixed[i]) continue;
			rowlist[nr] = i;
			vallist[nr] = localM[i];
			nr++;
		}
		_c(CPXchgrhs(env, lp, nr, rowlist, vallist));
		// solve the relaxation
		_c(CPXlpopt(env, lp));
		// get the solution
		double objval;
		int solstat;
		static double sol[2 * MAX_N];
		_c(CPXsolution(env, lp, &solstat, &objval, sol, NULL, NULL, NULL));
		// problem infeasable
		if (solstat != CPX_STAT_OPTIMAL) return;
		// cannot improve incumbent
		if (objval > incumbent) return;
		int sp_every = 5 + node_n/8096;
		bool speculative = incumbent == INFINITY || (node_n/1024)%sp_every==0;
		int branch = -1;
		double most_unfeas = 0;
		double tot_unfeas = 0;
		constexpr double EPS = 1E-6;
		for (int i=0; i<N; i++) {
			double x = sol[i];
			double unfeas = 0.5 - abs(x - 0.5);
			tot_unfeas += unfeas*unfeas;
			if (!speculative && !fixed[i]) unfeas *= localM[i] - localL[i];
			if (unfeas > most_unfeas) {
				most_unfeas = unfeas;
				branch = i;
			}
		}
		if (most_unfeas <= EPS) { // integer solution
			if (objval < incumbent) {
				incumbent = objval;
				copy_n(sol, N, inc_sol);
				print_info();
			}
		} else { // noninteger solution
			node_t next = {
				fixed, fvalue, objval, branch, tot_unfeas
			};
			if (speculative)
				f_queue.push(next);
			else
				b_queue.push(next);
			++glob_bound[objval];
		}

	};

	search_node(0, 0);

	while (b_queue.size() + f_queue.size()) {
		if (PARAM_TIME_LIMIT != -1) {
			_c(CPXgettime(env, &t_end));
			if (t_end - t_start > PARAM_TIME_LIMIT) {
				print_info();
				break;
			}
		}

		node_t node;
		if (!f_queue.empty()) {
			node = f_queue.top();
			f_queue.pop();
		} else {
			node = b_queue.top();
			b_queue.pop();
		}

		auto itr = glob_bound.find(node.obj);
		assert(itr != glob_bound.end());
		if (--itr->second == 0) glob_bound.erase(itr);

		if (node.obj > incumbent) { // we got a better incumbent in the meantime
			continue;
		}

		if (node_n && node_n % 500 == 0) {
			print_info();
		}
		node_n++;

		for (bool direction : {0, 1}) {
			nbitset f = node.fixed;
			f.set(node.branch);
			nbitset fval = node.fvalue;
			fval.set(node.branch, direction);
			search_node(f, fval);
		}
	}

	_c(CPXgettime(env, &t_end));
	double t_took  = t_end - t_start;
	cout << "solution = " << incumbent << endl;
	cout << "took " << t_took << "s" << endl;
	cerr << N << endl;
	for (int i=0; i<N; i++) {
		//cout << "  x_" << i << " = " << int(inc_sol[i_x(i)]) << endl;
		cerr << int(round(inc_sol[i_x(i)])) << endl;
	}
	cerr << "# obj = " << incumbent << endl;
	cerr << "# bound = " << curr_bound() << endl;
	cerr << "# nodes = " << node_n << endl;
	cerr << "# time = " << t_took << endl;
	cerr << "# custom" << endl;

    // clean up
    _c(CPXfreeprob(env, &lp));
    _c(CPXcloseCPLEX(&env));

    return 0;
} 
