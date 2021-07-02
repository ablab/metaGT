#!/usr/bin/python
import sys
from yaml import load
import yaml

dataset_data = load(open(sys.argv[1]), Loader=yaml.FullLoader)
all_reads = []
for reads_library in dataset_data:
        for key, value in reads_library.items():
            if key.endswith("reads"):
                for reads_file in value:
                    all_reads.append(reads_file)

s = ' '.join(all_reads)
print(s)