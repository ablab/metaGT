#!/usr/bin/env python 

from Bio import SeqIO
import sys


fasta = SeqIO.parse(sys.argv[1], 'fasta')
new=[]
for record in fasta:
    new_name = record.id.split('_')[0:4]
    print('_'.join(new_name))
    record.id='_'.join(new_name)
    new.append(record)
SeqIO.write(new, sys.argv[2], 'fasta')
