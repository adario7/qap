import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys, math

if __name__=='__main__':
	assert(len(sys.argv) == 2)
	f = open(sys.argv[1], 'r')

	params = f.readline().split(' ')
	if len(params) == 2:
		n1 = int(params[0])
		n2 = int(params[1])	
		n = n1*n2
	else:
		n = int(params[0])
		n1, n2 = int(math.sqrt(n)), int(math.sqrt(n))

	F = np.zeros((n1, n2), dtype=float)
	S = np.zeros((n1, n2), dtype=float)
	for i in range(n):
		v = float(f.readline())
		F[int(i/n2), i%n2] = v
	for i in range(n):
		v = float(f.readline())
		S[int(i/n2), i%n2] = v

	f.close()


	fig, (ax1, ax2) = plt.subplots(1,2)
	
	sns.heatmap(F, vmin=0, vmax=1, linewidths=.5, linecolor='#000', cmap='Greys', cbar=False, square=True, xticklabels=False, yticklabels=False, ax=ax1)
	#plt.savefig(f'images/tai{n}c_{n1}x{n2}_{d}.eps', format='eps', bbox_inches='tight')
	sns.heatmap(S, vmin=0, vmax=1, linewidths=.5, linecolor='#000', cmap='Greys', cbar=False, square=True, xticklabels=False, yticklabels=False, ax=ax2)
	plt.show()
