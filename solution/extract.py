import os
import re
import json

dir_path = 'build/data'
output_file = 'build/raw.json'

pattern = r'e_(?P<program>\w+?)_(?P<name>\w+?)_(?P<size>\w+?)_(?P<density>\w+?)_s_(?P<seed>\w+?)(_(?P<params>.*))?'

all_maps = []

for filename in os.listdir(dir_path):
	match = re.match(pattern, filename)
	if not match: 
		print(f"unamtched: {filename}")
		continue
	map = match.groupdict()
	map["filename"] = filename
	if map["params"] is None:
		map["params"] = ""

	with open(os.path.join(dir_path, filename), 'r') as file:
		for line in file:
			# skip lines starting with '#'
			if not line.startswith('#'): continue
			line = line[1:]
			# split on commas and get 'key=value' pairs
			pairs = line.split(',')
			for pair in pairs:
				arr = pair.split('=')
				if len(arr) == 2:
					key, value = arr
					map[key.strip()] = value.strip()

	if map['program']=='custom' and map['params']=='': # backw compat
		map['params'] = 'b_1'

	all_maps.append(map)

with open(output_file, 'w') as file:
	file.write(json.dumps(all_maps).replace("}, {", "},\n{"))
