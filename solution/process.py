import pandas as pd
import numpy as np
import json

with open('build/raw.json', 'r') as f:
	data = json.load(f)

df = pd.DataFrame(data)

# Discard rows without a 'status' column
df['status'] = df['status'].fillna('<failed>')
print(df.groupby('status')['status'].count())
df = df.dropna(subset=['status'])
df = df[(df['status'] == '102') | (df['status'] == '101')]
df = df.drop(columns=['filename', 'tot l / m / lp / la / f cuts'])

obj_grouped = df.groupby(['name', 'size', 'density'])

# assert all configuration found the same solution
for name, group in obj_grouped:
	assert len(group['obj'].unique()) == 1, f"Group {name} has multiple 'obj' values: {group['obj']}"

df['density'] = pd.to_numeric(df['density'], errors='coerce')
df['n'] = df['name'].str.extract('tai(\\d+)c').astype(int)
df['m'] = np.round(df['n'] * df['density'] / 100).astype(int)

key_columns = ['program', 'name', 'size', 'density', 'params']

non_key_columns = [col for col in df.columns if col not in key_columns]

for col in non_key_columns:
	df[col] = pd.to_numeric(df[col], errors='raise')

grouped = df.groupby(key_columns)
df_grouped = grouped.mean()

df_grouped['samples'] = grouped.size()

# warning if group size is != 3 -> some files missing
for key, group in grouped:
	if len(group) != 3:
		print(f"Warning: Group {key} has {len(group)} rows")


df_grouped.to_csv("build/data.csv")
