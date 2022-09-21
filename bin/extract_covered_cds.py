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


def count_coverage(chr_id, start, end, bam, gene_seq):
    gene_len = end - start + 1
    coverage = [0 for _ in range(gene_len)]
    for a in bam.fetch(chr_id, start - 1, end):
        covered_start = max(a.reference_start + 1, start)
        covered_end = min(a.reference_end + 1, end)
        i = -1
        aligned_pairs = a.get_aligned_pairs(with_seq = True)
        for pos in range(covered_start, covered_end):
            coverage[pos - start] += 1
            i += 1
            if aligned_pairs[covered_start-a.reference_start+a.query_alignment_start-1+i][2] is None:
                continue
            else:
                gene_seq[pos - start] = aligned_pairs[covered_start-a.reference_start+a.query_alignment_start-1+i][2]

    return 1 - (coverage.count(0) / gene_len), gene_seq


def get_gene_stats(gene_db, bam, genome_fasta, threshold, file):
    gene_cov_dict = {}
    gene_seq_dict = {}
    strand = True
    for g in gene_db.features_of_type('CDS', order_by=('seqid', 'start')):
        gene_name = g.id
        if g.strand == '-':
            strand = False
        gene_seq = list(g.sequence(genome_fasta, use_strand = strand))
        #print(g.id, g.seqid, g.start, g.end)
        gene_cov, gene_seq = count_coverage(g.seqid, g.start, g.end, bam, gene_seq)
        gene_cov_dict[gene_name] = gene_cov
        if gene_cov > threshold:
            extract_genes(gene_seq, g, file)
    print("Genes processed %d" % len(gene_cov_dict))
    return gene_cov_dict


def extract_genes(gene_seq, gene, file):
    name = ''
    if 'Name' in gene.attributes:
        name = gene.attributes["Name"][0]

    rec = SeqRecord(Seq("".join(str(x) for x in gene_seq)), id=gene.id, description=name)
    SeqIO.write(rec, file, "fasta")


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
    gene_cov_dict = get_gene_stats(gffutils_db, bam, genome, t, f)
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