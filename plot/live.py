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

	def mkplot(ax, vmin=None, vmax=None):
		S = np.zeros((n1, n2), dtype=float)
		for i in range(n):
			v = float(f.readline())
			S[int(i/n2), i%n2] = v
		sns.heatmap(S, vmin=vmin, vmax=vmax, linewidths=.5, linecolor='#000', cmap='Greys', cbar=False, square=True, xticklabels=False, yticklabels=False, ax=ax)


	fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12, 4))
	mkplot(ax1, 0, 1)
	mkplot(ax2, 0, 1)
	mkplot(ax3)
	plt.show()

	f.close()
