import numpy as np
import argparse

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
def B_generator(n, n1, n2):
    # in our case, we suppose n to be a perfect square
    B = np.zeros((n, n), dtype=int)
    # (note: B is always simmetric)
    scale = 100000
    for r in range(n1):
        for s in range(n2):
            for t in range(n1):
                for u in range(n2):
                    i = n2*(r-1)+s
                    j = n2*(t-1)+u
                    if(B[i, j] != 0):
                        continue
                    else:
                        B[i, j] = round(force(n1, n2, r, s, t, u)*scale)
                        B[j, i] = B[i, j]
    return B

def gen_one(n1, n2, d):
    # Calculate n
    n = n1 * n2

    # Generate matrices A and B
    A = A_generator(n, d)
    B = B_generator(n, n1, n2)

    # Write matrices to file
    filename = f"tai{n}c_{n1}x{n2}_{d}"
    with open(filename, 'w') as file:
        file.write(f"{n}\n")
        for row in A:
            file.write(' '.join(map(str, row)) + '\n')
        for row in B:
            file.write(' '.join(map(str, row)) + '\n')

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Generate matrices A and B')
    parser.add_argument('n1', type=int, help='Integer n1')
    parser.add_argument('n2', type=int, help='Integer n2')
    parser.add_argument('d', type=int, help='Integer d')
    args = parser.parse_args()

    gen_one(args.n1, args.n2, args.d)

if __name__ == "__main__":
    main()
