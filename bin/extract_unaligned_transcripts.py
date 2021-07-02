#!/usr/bin/python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import pysam as ps
import sys

align = ps.AlignmentFile(sys.argv[1], "r")
table=[]
for read in align.fetch():
    if read.flag==4:
            table.append(read.qname)
align.close()

new_table = []
for i in table:
    if int(i.split("_")[3]) > 300 and float(i.split("_")[5]) > 5:
        new_table.append(i)

transcripts = SeqIO.parse(sys.argv[2], 'fasta')

new_transcripts = []
for i in transcripts:
    if i.id in new_table:
        new_transcripts.append(i)

SeqIO.write(new_transcripts, sys.argv[3], "fasta")