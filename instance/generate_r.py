import numpy as np
import argparse
import random

# (Density of grey)
# Function that returns the force value (used later as a distance) between two given unit locations of the matrix: n1 x n2
def force(n1, n2, r, s, t, u):
	if (r,s) == (t, u): return 0
	return max(
		1/((r-t+n1*v)**2 + (s-u+n2*w)**2)
		for v in {-1, 0, 1}
		for w in {-1, 0, 1}
	)

# Function that generates the flows matrix (A)
def A_generator(n, d):
	m = round(n*(d/100))
	A = [0]*n
	for i in range(n):
		if(i < m):
			A[i] = [1]*m + [0]*(n-m)
		else:
			A[i] = [0]*n
	return A

# Function that generates the distances matrix (B)
def B_generator(n):
	# in our case, we suppose n to be a perfect square
	B = np.zeros((n, n), dtype=int)
	# (note: B is always simmetric)
	scale = 100000
	for i in range(n):
		for j in range(n):
			if(B[i, j] != 0):
				continue
			else:
				B[i, j] = random.randint(0, scale)
				B[j, i] = B[i, j]
	return B


def main():
	# Parse command line arguments
	parser = argparse.ArgumentParser(description='Generate matrices A and B')
	parser.add_argument('n', type=int, help='Integer n')
	parser.add_argument('d', type=int, help='Integer d')
	args = parser.parse_args()

	# Calculate n
	n = args.n

	# Generate matrices A and B
	A = A_generator(n, args.d)
	B = B_generator(n)

	# Write matrices to file
	filename = f"tai{n}c_rand_{args.d}.txt"
	with open(filename, 'w') as file:
		file.write(f"{n}\n")
		for row in A:
			file.write(' '.join(map(str, row)) + '\n')
		for row in B:
			file.write(' '.join(map(str, row)) + '\n')

if __name__ == "__main__":
	main()
