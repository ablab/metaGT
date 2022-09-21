#!/usr/bin/python
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pysam
import gffutils
import sys
import argparse


def count_coverage(chr_id, start, end, bam):
    gene_len = end - start + 1
    coverage = [0 for _ in range(gene_len)]
    for a in bam.fetch(chr_id, start - 1, end):
        covered_start = max(a.reference_start + 1, start)
        covered_end = min(a.reference_end + 1, end)
        i = -1
        for pos in range(covered_start, covered_end):
            coverage[pos - start] += 1
            i += 1

    return 1 - (coverage.count(0) / gene_len)


def get_gene_stats(gene_db, bam, genome_fasta, threshold, file):
    gene_cov_dict = {}
    gene_records = []
    strand = True
    for g in gene_db.features_of_type('CDS', order_by=('seqid', 'start')):
        gene_name = g.id
        if g.strand == '-':
            strand = False
        #print(g.id, g.seqid, g.start, g.end)
        gene_cov = count_coverage(g.seqid, g.start, g.end, bam)
        gene_cov_dict[gene_name] = gene_cov
        if gene_cov > threshold:
            gene_seq = g.sequence(genome_fasta, use_strand=strand)
            gene_records.append(extract_genes(gene_seq, g))
    print("Genes processed %d" % len(gene_cov_dict))
    return gene_cov_dict, gene_records


def extract_genes(gene_seq, gene):
    name = ''
    if 'Name' in gene.attributes:
        name = gene.attributes["Name"][0]

    rec = SeqRecord(Seq(gene_seq), id=gene.id, description=name)
    return rec

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file", default="gtf_stats.tsv")
    parser.add_argument("--gff", "-g", type=str, help="gene db", required=True)
    parser.add_argument("--bam", "-b", type=str, help="bam to count coverage", required=True)
    parser.add_argument("--genome", "-f", type=str, help="initial file with covered genes", required=True)
    parser.add_argument("--threshold", "-t", type=float, help="threshold", default=0.5)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    gffutils.create_db(args.gff, 'db')
    gffutils_db = gffutils.FeatureDB('db', keep_order=True)
    bam = pysam.AlignmentFile(args.bam, "rb")
    genome= args.genome
    t = float(args.threshold)
    f = open(args.output + '.fasta', 'w')
    gene_cov_dict, gene_records = get_gene_stats(gffutils_db, bam, genome, t, f)
    SeqIO.write(gene_records, f, "fasta")
    with open(args.output+'.csv', 'w') as outf:
        for k in sorted(gene_cov_dict.keys()):
            outf.write("%s\t%.3f\n" % (k, gene_cov_dict[k]))
    f.close()

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)