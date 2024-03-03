#include <algorithm>
#include <cassert>
#include <iostream>

#include <inputs.hh>

using namespace std;

uint N, M;
double B[MAX_N][MAX_N];

void read_inputs() {
	cin >> N;
	assert(N <= MAX_N);
	static double A[MAX_N][MAX_N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			cin >> A[i][j];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			cin >> B[i][j];
	M = count(A[0], A[0]+N, 1.0);
	cout << "N = " << N << ", M = " << M << endl;
	assert(M > 0);
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++) {
		double a = A[i][j];
		bool t1 = i < M && j < M;
		if ((t1 && a != 1) || (!t1 && a != 0)) {
			cerr << "not a type C instance! A[" << i << "][" << j << "] = " << a << endl;
			abort();
		}
	}
}
