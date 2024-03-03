import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys, math

if __name__=='__main__':
	assert(len(sys.argv) == 2)
	f = open(sys.argv[1], 'r')

	params = f.readline().split(' ')
	if len(params) == 3:
		n1 = int(params[0])
		n2 = int(params[1])	
		n = n1*n2
		d = int(params[2])
	else:
		n = int(params[0])
		n1, n2 = int(math.sqrt(n)), int(math.sqrt(n))
		d = int(100 * float(params[1]) / n)


	S = np.zeros((n1, n2), dtype=int)
	for i in range(n):
		l = f.readline().split(' ')
		S[int(i/n2), i%n2] = int(l[1])

	f.close()
	
	sns.heatmap(S, linewidths=.5, linecolor='#bcbcbc', cmap='Greys', cbar=False, square=True, xticklabels=False, yticklabels=False)
	#plt.savefig(f'images/tai{n}c_{n1}x{n2}_{d}.eps', format='eps', bbox_inches='tight')
	plt.show()
