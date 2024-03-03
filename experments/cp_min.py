from pyscipopt import Model, quicksum
import matplotlib.pyplot as plt
import numpy as np

N1, N2 = 5, 5
N = N1 * N2
D = 0.2
M = round(N * D)
print(f"grid {N1}*{N2} = {N}, d={D}, m={M}")

model = Model()

x = [ model.addVar(f"x_{i}", lb=0, ub=N1-1, vtype="C") for i in range(0,M)  ]
y = [ model.addVar(f"y_{i}", lb=0, ub=N2-1, vtype="C") for i in range(0,M)  ]

def sq(x):
	return x*x

def dist_sq(i, j, v, w):
	return sq(x[i] - x[j] + N1*v) + sq(y[i] - y[j] + N2*w)

for i in range(0, M):
	for j in range(i+1, M):
		model.addCons(dist_sq(i, j, 0, 0) >= 1)

dists = []
for i in range(0, M):
	for j in range(i+1, M):
		neigh = [
			1.0 / dist_sq(i, j, v, w)
			for v in [-1, 0, 1]
			for w in [-1, 0, 1]
		]
		z = model.addVar(f"z_{i}_{j}")
		for k in neigh:
			model.addCons(z >= k)
		dists.append(z)

cost = quicksum(
	z for z in dists
)

z = model.addVar("z")
model.addCons(z >= cost)
model.setObjective(z, sense="minimize")

model.setRealParam('expr/pow/minzerodistance', .9)  # Adjust the value as needed
model.optimize()

sol = model.getBestSol()
sx = [ sol[x[i]] for i in range(0,M) ]
sy = [ sol[y[i]] for i in range(0,M) ]

for i in range(0,M):
	print(f"{sx[i]} {sy[i]}")

#def plot_points_on_grid(N1, N2, sx, sy):
# Create a 2D grid
x = np.linspace(0, N1-1, N1)
y = np.linspace(0, N2-1, N2)
X, Y = np.meshgrid(x, y)

# Plot the grid
plt.figure(figsize=(8, 6))
plt.plot(X, Y, marker='.', color='k', linestyle='none')

# Plot the points
plt.scatter(sx, sy, color='r', label='Points')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('2D Grid with Points')
plt.legend()
plt.grid(True)
plt.show()
