import sys

inhars = sys.argv[1]
outdir = sys.argv[2]

with open(inhars, 'r') as f:
    for line in f.readlines():
        har_name = line.strip().split('\t')[-1]
        with open(f'{outdir}/{inhars.split("/")[-1].split(".bed")[0]}.{har_name}.bed','w') as write:
            write.write(line)