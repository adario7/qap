from pyscipopt import Model, quicksum
import sys

N = int(input())
A = [[0.0] * N for _ in range(N)]
for i in range(N):
	A[i] = list(map(float, input().split()))
B = [[0.0] * N for _ in range(N)]
for i in range(N):
	B[i] = list(map(float, input().split()))
M = sum(1 for val in A[0] if val == 1.0)
print("N = {}, M = {}".format(N, M))
assert M > 0
for i in range(N):
	for j in range(N):
		a = A[i][j]
		t1 = i < M and j < M
		if (t1 and a != 1) or (not t1 and a != 0):
			sys.stderr.write("not a type C instance! A[{}][{}] = {}\n".format(i, j, a))
			sys.exit(1)

model = Model()

x = [ model.addVar(f"x_{i}", vtype="B") for i in range(0,N)  ]

model.addCons(quicksum(x) == M)

cost = quicksum(
	B[i][j] * x[i] * x[j]
	for i in range(N)
	for j in range(N)
)

z = model.addVar("z")
model.addCons(z >= cost)
model.setObjective(z, sense="minimize")

model.setRealParam('expr/pow/minzerodistance', .9)  # Adjust the value as needed
model.optimize()

sol = model.getBestSol()
for i in range(N):
	print(f"x_{i} = {round(sol[x[i]])}")
