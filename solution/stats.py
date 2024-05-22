import pandas as pd
from scipy.stats import gmean
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.special
import numpy as np
import matplotlib.colors as mcolors

df0 = pd.read_csv('build/data.csv')
df0['params'] = df0['params'].fillna('default')

rename_prog = { 'cb': 'kb', 'toff': 'intuitivo', 'custom': 'kb ad-hoc' }
rename_par = { 'default': 'default',
			  'l_0_m_0': 'no-cuts', 'm_0': 'no-m',
			  'p_1': 'p', 'a_1': 'a', 'f_1': 'f', 'r_1': 'r1', 'r_3': 'r3', 'b_1': 'b1', 'b_2': 'b2', 'df_1': 'dfs',
			  'f_1_r_1_b_2': 'f+r1+b2' }
df0['program'] = df0['program'].apply(lambda x: rename_prog[x])
df0['params'] = df0['params'].apply(lambda x: rename_par[x])

df0.sort_values(['n', 'density'], inplace=True)
df0['size'] = pd.Categorical(df0['size'], categories=df0['size'].unique(),ordered=True)

configuration = ['program', 'params']
instance = ['name', 'size', 'density']

df0['configuration'] = df0.apply(lambda row: f"{row['program']} - {row['params']}" if row['program'] == 'kb' else row['program'], axis=1)


# solved config heatmap
df = df0.copy()
df_counts = df.groupby(['size', 'density'], sort=False)['configuration'].apply(len).reset_index()
df_pivot = df_counts.pivot('density', 'size', 'configuration').fillna(0)
print(df_pivot)
plt.figure(figsize=(10, 8))
sns.heatmap(df_pivot, annot=True, fmt=".0f", cbar=False, cmap="Blues_r")
plt.ylabel("Densità")
plt.xlabel("Dimesione")
plt.gca().invert_yaxis()
plt.savefig("build/tot-solved.eps")


# solved time heatmap
df = df0.copy()
df_counts = df.groupby(['size', 'density'], sort=False)['time'].min().reset_index()
df_pivot = df_counts.pivot('density', 'size', 'time')
print(df_pivot)
plt.figure(figsize=(10, 8))
sns.heatmap(df_pivot, annot=True, fmt=".1f", cbar=False, cmap="Blues", norm=mcolors.LogNorm())
plt.ylabel("Densità")
plt.xlabel("Dimesione")
plt.gca().invert_yaxis()
plt.savefig("build/tot-time.eps")


presets = [
	('base', [ ['kb', 'no-cuts'], ['kb', 'no-m'], ['kb', 'default'] ]),
	('variants', [ ['kb', 'default'], ['kb', 'p'], ['kb', 'a'], ['kb', 'f'], ['kb', 'r1'], ['kb', 'r3'], ['kb', 'b1'], ['kb', 'b2'], ['kb', 'dfs'], ['kb', 'f+r1+b2'] ]),
	('custom', [ ['intuitivo', 'default'], ['kb', 'f+r1+b2'], ['kb ad-hoc', 'default'] ]),
]

def pretty_df(df, pconf):
	# sort
	order_df = pd.DataFrame(pconf, columns=['program', 'params'])
	order_df['order'] = range(len(order_df))
	merged_df = df.merge(order_df, on=['program', 'params'], how='left')
	df = merged_df.sort_values('order').drop(columns='order')
	# format
	df['params'] = df['params'].apply(lambda x: "\\code{" + x + "}")
	bt, bn, bs = df['gmean_time'].min(), df['gmean_nodes'].min(), df['solved'].max()
	df['gmean_time'] = df['gmean_time'].apply(lambda x: "%.2f" % x if x != bt else "\\textbf{%.2f}" % x)
	df['gmean_nodes'] = df['gmean_nodes'].apply(lambda x: "%.0f" % x if x != bn else "\\textbf{%.0f}" % x)
	df['solved'] = df['solved'].apply(lambda x: "%.0f" % x if x != bs else "\\textbf{%.0f}" % x)
	# rename
	df = df[['program', 'params', 'gmean_time', 'gmean_nodes', 'solved']]
	df = df.rename(columns={'program': 'Implementazione', 'params': 'Variante', 'gmean_time': 'Tempo medio (s)', 'gmean_nodes': 'Nodi medi', 'solved': 'Istanze risolte'})
	if df['Implementazione'].nunique() == 1:
		df = df.drop(columns=['Implementazione'])
	return df

def set_name(set_configs, x):
	s = sorted(set(x))
	if len(s) == 0: return "irrisolta"
	if len(s) == len(set_configs): return "tutte"
	if len(s) <= 2 or len(s) < len(set_configs)/2: return ', '.join(s)
	return f"eccetto {', '.join(sorted(set_configs - set(x)))}"

for pname, pconf in presets:
	df = df0[df0[configuration].apply(tuple, axis=1).isin(map(tuple, pconf))]

	# solving sets heatmap
	set_configs = set(df['configuration'].drop_duplicates())
	pivot_df = df.groupby(['size', 'density'])['configuration'].apply(lambda x: set_name(set_configs, x)).reset_index()
	pdf = pivot_df.pivot('size', 'density', 'configuration').fillna("irrisolta")
	unique_strings = pd.unique(pdf.values.ravel())
	string_to_num = {string: i for i, string in enumerate(unique_strings)}
	num_to_string = {i: string for string, i in string_to_num.items()}
	numeric_data = pdf.replace(string_to_num)
	colors = sns.color_palette("hsv", len(unique_strings))
	cmap = mcolors.ListedColormap(colors)
	plt.figure(figsize=(7, 4))
	heatmap = sns.heatmap(numeric_data, fmt='', cmap=cmap, cbar=False)
	handles = [plt.Line2D([0], [0], marker='o', color=colors[string_to_num[string]], label=string, linestyle='') for string in unique_strings]
	plt.legend(handles=handles, title='Risolta da', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
	plt.subplots_adjust(right=0.6)
	plt.savefig(f"build/{pname}-solved.eps")

	# num solved
	grouped = df.groupby(configuration)
	num_solved = grouped.size().reset_index(name='solved')
	all_configs = df[configuration].drop_duplicates()
	solved_by_all = df.groupby(instance).filter(lambda x: len(x[configuration]) == len(all_configs))[instance].drop_duplicates()
	print(f"solved_by_all: {len(solved_by_all)}")
	solved_by_all.to_csv(f"build/{pname}-solved-by-all.csv", index=False)

	print()

	# mean values
	merged_df = pd.merge(df, solved_by_all, on=instance)
	grouped = merged_df.groupby(configuration)
	SHIFT = 1
	mean_time = grouped['time'].apply(lambda x: gmean(x+SHIFT)-SHIFT).reset_index(name='gmean_time')
	SHIFT = 1000
	mean_nodes = grouped['nodes'].apply(lambda x: gmean(x+SHIFT)-SHIFT).reset_index(name='gmean_nodes')
	means = pd.merge(mean_time, mean_nodes, on=configuration)
	means = pd.merge(means, num_solved, on=configuration)

	print(means)
	pretty_df(means, pconf).to_csv(f"build/{pname}-means.csv", index=False)

	# instance size vs time plot
	pdf = df[df['density'].isin([20, 30, 40, 50])]
	pdf = pdf[pdf['params'].isin(['default', 'no-cuts', 'no-m', 'f+r1+b2'])]
	plt_binomial = pdf['binomial'] = pdf.apply(lambda row: scipy.special.comb(row['n'], row['m']), axis=1)
	plt_config = pdf['configuration'] = pdf.apply(lambda row: f"{row['program']} - {row['params']}", axis=1)
	plt.figure(figsize=(10, 6))
	sns.scatterplot(x=plt_binomial, y=pdf['time'], hue=plt_config, palette='deep')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Binomiale (n, m)')
	plt.ylabel('Tempo di esecuzione')
	plt.legend(title='Configurazione')
	# linear fit
	for config in pdf['configuration'].unique():
		subset = pdf[pdf['configuration'] == config]
		x = np.log10(subset['binomial'])
		y = np.log10(subset['time'])
		slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
		plt.plot(subset['binomial'], 10**(intercept + slope * np.log10(subset['binomial'])), label=f'{config} fit')
	plt.grid(True, which="both", ls="--")
	plt.savefig(f"build/{pname}-binomial.eps")


	#if pname == "variants":
		#for otherconf in pconf[1:]:
			#pdf = df[df['density'].isin([20, 30, 40, 50])]
			## plot only some configs
			#pdf = pdf[pdf['params'].isin(['default', otherconf[1]])]
			## Calculate binomial coefficient for (n, m)
			#plt_binomial = pdf['binomial'] = pdf.apply(lambda row: scipy.special.comb(row['n'], row['m']), axis=1)
			#plt_config = pdf['configuration'] = pdf.apply(lambda row: f"{row['program']} - {row['params']}", axis=1)
			## Plot
			#plt.figure(figsize=(10, 6))
			#sns.scatterplot(x=plt_binomial, y=pdf['time'], hue=plt_config, palette='deep')
			#plt.title(otherconf)
			#plt.xscale('log')
			#plt.yscale('log')
			#plt.xlabel('Binomiale (n, m)')
			#plt.ylabel('Tempo di esecuzione')
			#plt.legend(title='Configurazione')
			#for config in pdf['configuration'].unique():
				#subset = pdf[pdf['configuration'] == config]
				#x = np.log10(subset['binomial'])
				#y = np.log10(subset['time'])
				#slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
				#plt.plot(subset['binomial'], 10**(intercept + slope * np.log10(subset['binomial'])), label=f'{config} fit')
			#plt.grid(True, which="both", ls="--")
			#plt.show()



