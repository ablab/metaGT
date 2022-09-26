#!/usr/bin/env python 

from Bio import SeqIO
import sys

used_contigs = set()
for l in open(sys.argv[1]):
    used_contigs.add(l.strip())

transcripts = SeqIO.parse(sys.argv[2], 'fasta')
new_transcripts = []
for i in transcripts:
    if i.id not in used_contigs and  int(i.id.split("_")[3]) >= 300 and float(i.id.split("_")[5]) >= 1:
        new_transcripts.append(i)

SeqIO.write(new_transcripts, sys.argv[3], "fasta")
