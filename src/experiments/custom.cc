#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <mutex>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <ilcplex/cplex.h>

#include <inputs.hh>
#include <params.hh>

using namespace std;

#define _c(what) if (int _error = what) { \
	cout << "CPX error: " #what << endl; cout << "CPX error: " << _error << endl; abort(); }

using nbitset = bitset<MAX_N>;

struct node_t {
	uint id;
	nbitset fixed;
	nbitset fvalue;
	double obj;
	int branch;
	double unfeas;
};

struct cmp_feas { bool operator()(const node_t& a, const node_t& b) const {
	return a.unfeas > b.unfeas;
} };
struct cmp_bound { bool operator()(const node_t& a, const node_t& b) const {
	return a.obj > b.obj;
} };

struct node_queue {
	CPXENVptr env;
	double t_start;
	mutex monitor;
    condition_variable cv;
	priority_queue<node_t, vector<node_t>, cmp_bound> b_queue;
	priority_queue<node_t, vector<node_t>, cmp_feas> f_queue;
	map<double, int, less<double>> glob_bound;
	long tot_depth = 0;
	uint n_count = 0;
	uint n_working = 0;
	double incumbent = INFINITY;
	double inc_sol[MAX_N];
};

int B_order[MAX_N][MAX_N];
double rootL[MAX_N], rootM[MAX_N];

int i_x(int i) {
	assert(0 <= i && i < N);
	return i;
}

int i_w(int i) {
	assert(0 <= i && i < N);
	return N + i;
}

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

void calc_root_LM() {
	nbitset fixed0 = 0, fixed1 = 1;
	for (int i=0; i<N; i++) {
		rootL[i] = calc_local_L(i, fixed0, fixed1);
		rootM[i] = calc_local_M(i, fixed0, fixed1);
	}
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
		double Mi = rootM[i];
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

// w_i >= L_i x_i  <=>  w_i - L_i x_i >= 0 
void add_cons_L(CPXENVptr env, CPXLPptr lp) {
	vector<double> rhs(N, 0);
	vector<char> sense(N, 'G');
	vector<int> rmatbeg(N, 0);
	vector<int> rmatind(N*2, 0);
	vector<double> rmatval(N*2, 0);
    for (int i = 0; i < N; i++) {
		double Li = rootL[i];
		int beg = i * 2;
		rmatbeg[i] = beg;
		rmatind[beg + 0] = i_w(i);
		rmatval[beg + 0] = 1;
		rmatind[beg + 1] = i_x(i);
		rmatval[beg + 1] = -Li;
	}
    _c(CPXaddrows(env, lp, 0, N, N*2, rhs.data(), sense.data(), rmatbeg.data(), rmatind.data(), rmatval.data(), NULL, NULL));
}

void queue_print_info(node_queue& q) {
	double b = q.glob_bound.empty() ? q.incumbent : q.glob_bound.begin()->first;
	double t_end;
	_c(CPXgettime(q.env, &t_end));
	double gap = 100 * (q.incumbent - b) / q.incumbent;
	double avgd = double(q.tot_depth) / q.glob_bound.size();
	cout << int(t_end-q.t_start) << "s /\t "<< q.n_count << " /\t " << (q.b_queue.size()+q.f_queue.size()) << " /\t " << avgd << " /\t " << q.incumbent << " /\t " << b << " /\t " << gap << "%" << endl;
}

void forget_bound(node_queue& q, const node_t& n) {
	q.tot_depth -= n.fixed.count();
	auto itr = q.glob_bound.find(n.obj);
	assert(itr != q.glob_bound.end());
	if (--itr->second == 0) q.glob_bound.erase(itr);
}

bool queue_request(node_queue& q, vector<node_t>& out, double& local_inc) {
	unique_lock<mutex> lock(q.monitor);
	q.cv.wait(lock, [&q] {
		bool empty = q.b_queue.empty() && q.f_queue.empty();
		return !empty || (empty && q.n_working == 0);
	});
	local_inc = q.incumbent;
	out.clear();
	int size = q.b_queue.size() + q.f_queue.size();
	if (size == 0 && q.n_working == 0) return false;
	int take = clamp(size / 4, 4, 64);
	while((q.b_queue.size() || q.f_queue.size()) && out.size() < take) {
		node_t n;
		if (q.f_queue.size()) {
			n = q.f_queue.top();
			q.f_queue.pop();
		} else if (q.b_queue.size()) {
			n = q.b_queue.top();
			q.b_queue.pop();
		} else abort();
		if (n.obj > q.incumbent) {
			forget_bound(q, n);
			continue; // we got a better incumbent in the meantime
		}
		out.push_back(n);
	}
	q.n_working++;
	return true;
}

void queue_submit(node_queue& q, const vector<node_t>& input, vector<node_t>& data, double local_inc, double* local_sol) {
	unique_lock<mutex> lock(q.monitor);
	for (const node_t& n : input)
		forget_bound(q, n);
	q.n_working--;
	bool info = false;
	for (auto& n : data) {
		if (n.obj >= q.incumbent) continue;
		n.id = q.n_count++;
		int sp_every = 5 + q.n_count/8096;
		bool speculative = q.incumbent == INFINITY || (q.n_count/1024)%sp_every==0;
		if (speculative) q.f_queue.push(n);
		else q.b_queue.push(n);
		q.glob_bound[n.obj]++;
		q.tot_depth += n.fixed.count();
		if (n.id && n.id % 5000 == 0) info = true;
	}
	if (local_inc < q.incumbent) {
		info = true;
		q.incumbent = local_inc;
		copy_n(local_sol, N, q.inc_sol);
	}
	if (info)
		queue_print_info(q);
	q.cv.notify_all();
}

CPXLPptr create_lp(CPXENVptr env) {
	int status;
    CPXLPptr lp = CPXcreateprob(env, &status, "qap");
	if (status) abort();
    _c(CPXchgobjsen(env, lp, CPX_MIN));

	add_variables_x(env, lp);
	add_variables_w(env, lp);

	// conss [0 -> N)
	add_cons_xw(env, lp);
	// conss [N -> 2N)
	add_cons_L(env, lp);
	// cons 2N
	add_cons_M(env, lp);

	return lp;
}

void eval_state(CPXENVptr env, CPXLPptr lp, nbitset fixed, nbitset fvalue, vector<node_t>& out_buf, double& local_inc, double* local_sol) {
	// fix variables bounds
	int bound_indices[MAX_N];
	char bound_lu[MAX_N];
	double bound[MAX_N];
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
	for (int i=0; i<N; i++)  {
		localL[i] = fixed[i] ? rootL[i] : calc_local_L(i, fixed0, fixed1);
		localM[i] = fixed[i] ? rootM[i] : calc_local_M(i, fixed0, fixed1);
	}
	int nr = 0;
	for (int i=0; i<N; i++) {
		// w_i - L_i x_i >= 0, conss [N -> 2N)
		rowlist[nr] = N + i;
		collist[nr] = i_x(i);
		vallist[nr] = -localL[i];
		nr++;
	}
	for (int i=0; i<N; i++) {
		// sum_{j!=i}( b_ij x_j ) + (M_i + b_ii) x_i - w_i <= M_i, conss [0 -> N)
		rowlist[nr] = i;
		collist[nr] = i_x(i);
		vallist[nr] = localM[i] + B[i][i];
		nr++;
	}
	_c(CPXchgcoeflist(env, lp, nr, rowlist, collist, vallist));
	nr = 0;
	for (int i=0; i<N; i++) {
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
	double sol[2 * MAX_N];
	_c(CPXsolution(env, lp, &solstat, &objval, sol, NULL, NULL, NULL));
	// problem infeasable
	if (solstat != CPX_STAT_OPTIMAL) return;
	// cannot improve incumbent
	if (objval > local_inc) return;
	int branch = -1;
	double best_score = 0;
	double most_unfeas = 0, tot_unfeas = 0;
	constexpr double EPS = 1E-6;
	for (int i=0; i<N; i++) {
		double x = sol[i_x(i)];
		double unfeas = 0.5 - abs(x - 0.5);
		most_unfeas = max(most_unfeas, unfeas);
		tot_unfeas += unfeas*unfeas;
		double gap;
		if (PARAM_BRANCHING == 2)  {
			gap = -sol[i_w(i)];
			for (int j=0; j<N; j++) gap += B[i][j] * x * sol[i_x(j)];
		} else {
			gap = localM[i] - localL[i];
		}
		double score = unfeas * gap;
		if (score > best_score) {
			best_score = score;
			branch = i;
		}
	}
	if (most_unfeas <= EPS) { // integer solution
		if (objval < local_inc) {
			local_inc = objval;
			copy_n(sol, N, local_sol);
		}
	} else { // noninteger solution
		out_buf.push_back({
			.fixed = fixed, .fvalue = fvalue,
			.obj = objval, .branch = branch, .unfeas = tot_unfeas
		});
	}
}

void expand_node(CPXENVptr env, CPXLPptr lp, const node_t& node, vector<node_t>& out_buf, double& local_inc, double* local_sol) {
	for (bool direction : {0, 1}) {
		nbitset f = node.fixed;
		f.set(node.branch);
		nbitset fval = node.fvalue;
		fval.set(node.branch, direction);
		eval_state(env, lp, f, fval, out_buf, local_inc, local_sol);
	}
}

int main(int argc, char** argv) {
	read_parameters(argc, argv);
	read_inputs();
	preorder_B();
	calc_root_LM();

	int status;
    CPXENVptr env = CPXopenCPLEX(&status);
	if (status) abort();
	_c(CPXsetintparam(env, CPXPARAM_Threads, 1));

	// write problem
	//_c(CPXwriteprob(env, lp, "p.lp", NULL));

	double t_start;
	_c(CPXgettime(env, &t_start));

	node_queue nq = {
		.env = env,
		.t_start = t_start,
	};

	int nthreads = PARAM_SINGLE_THREAD ? 1 : 4;

	nq.n_working = 1;
	vector<thread> threads;
	for (int th_id=0; th_id<nthreads; th_id++) {
		threads.emplace_back([th_id, &nq]() {
			int status;
			CPXENVptr env = CPXopenCPLEX(&status);
			if (status) abort();
			_c(CPXsetintparam(env, CPXPARAM_Threads, 1));
			CPXLPptr lp = create_lp(env);
			double local_inc = INFINITY;
			double local_sol[MAX_N];
			vector<node_t> read_buf, out_buf;
			if (th_id == 0) {
				eval_state(env, lp, 1, 1, out_buf, local_inc, local_sol); // assume x_0=1 due to symmetry w.l.o.g.
				queue_submit(nq, read_buf,out_buf, local_inc, local_sol);
			}
			while (queue_request(nq, read_buf, local_inc)) {
				out_buf.clear();
				for (node_t& n : read_buf)
					expand_node(env, lp, n, out_buf, local_inc, local_sol);
				queue_submit(nq, read_buf, out_buf, local_inc, local_sol);
			}
			_c(CPXfreeprob(env, &lp));
			_c(CPXcloseCPLEX(&env));
		});
	}
	for (thread& th : threads)
		if (th.joinable()) th.join();
	

	double t_end;
	_c(CPXgettime(env, &t_end));
	double t_took  = t_end - t_start;
	cout << "solution = " << nq.incumbent << endl;
	cout << "took " << t_took << "s" << endl;
	cout << "nodes = " << nq.n_count << endl;
	cerr << N << endl;
	for (int i=0; i<N; i++) {
		cerr << int(round(nq.inc_sol[i_x(i)])) << endl;
	}
	cerr << "# obj = " << nq.incumbent << endl;
	//cerr << "# bound = " << curr_bound << endl;
	cerr << "# nodes = " << nq.n_count << endl;
	cerr << "# time = " << t_took << endl;
	cerr << "# custom" << endl;

    // clean up
    _c(CPXcloseCPLEX(&env));

    return 0;
} 
