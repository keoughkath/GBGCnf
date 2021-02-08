import pandas as pd
import sys

har_ss = sys.argv[1]
flank_ss = sys.argv[2]
har = sys.argv[3]

with open(har_ss) as f:
    har_file_lines = f.readlines()
    species_line = har_file_lines[4]
    har_species = species_line.split('= ')[1].split(',')

with open(flank_ss) as f:
    har_file_lines = f.readlines()
    species_line = har_file_lines[4]
    har_flank_species = species_line.split('= ')[1].split(',')

species_to_keep = list(set(har_species).intersection(set(har_flank_species)))

with open(f'species_matches/{har}.txt','w') as f:
	for name in species_to_keep:
		f.write(name + '\n')

