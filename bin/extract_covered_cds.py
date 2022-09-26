#!/usr/bin/env python 
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
import os
from traceback import print_exc
from multiprocessing import Pool


def count_coverage(chr_id, start, end, bam):
    gene_len = end - start + 1
    coverage = [0 for _ in range(gene_len)]
    contig_list = []
    for a in pysam.AlignmentFile(bam, "rb").fetch(chr_id, start - 1, end):
        contig_list.append(a.query_name)
        covered_start = max(a.reference_start + 1, start)
        covered_end = min(a.reference_end + 1, end)
        i = -1
        for pos in range(covered_start, covered_end):
            coverage[pos - start] += 1
            i += 1

    return 1 - (coverage.count(0) / gene_len), contig_list


def process_single_gene(g, bam, genome_fasta, threshold):
    gene_name = g.id
    strand = True
    if g.strand == '-':
        strand = False

    gene_cov, contig_list = count_coverage(g.seqid, g.start, g.end, bam)
    gene_rec = None
    if gene_cov > threshold:
        gene_seq = g.sequence(genome_fasta, use_strand=strand)
        gene_rec = extract_genes(gene_seq, g)
    return gene_name, gene_cov, gene_rec, contig_list


def process_genes(gene_db, bam, genome_fasta, threshold, threads):
    gene_cov_dict = {}
    gene_records = []
    used_contigs = set()
    pool = Pool(threads)
    results = pool.starmap(process_single_gene,
                           [(g, bam, genome_fasta, threshold)
                            for g in gene_db.features_of_type('CDS', order_by=('seqid', 'start'))])
    pool.close()
    pool.join()
    for res in results:
        if res[2] is not None:
            gene_records.append(res[2])
            used_contigs.update(res[3])
        gene_cov_dict[res[0]] = res[1]
    print("Genes processed %d" % len(gene_cov_dict))
    return gene_cov_dict, gene_records, used_contigs


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
    parser.add_argument("--threshold", "-t", type=float, help="threshold", default=0.25)
    parser.add_argument("--threads", "-p", type=int, help="threads", default=8)

    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    if not os.path.exists('db'):
        gffutils.create_db(args.gff, 'db')
    print("Loading genedb")
    gffutils_db = gffutils.FeatureDB('db', keep_order=True)
    genome= args.genome
    t = float(args.threshold)
    f = open(args.output + '.fasta', 'w')
    gene_cov_dict, gene_records, used_contigs = process_genes(gffutils_db, args.bam, genome, t, args.threads)
    SeqIO.write(gene_records, f, "fasta")
    with open(args.output+'.csv', 'w') as outf:
        for k in sorted(gene_cov_dict.keys()):
            outf.write("%s\t%.3f\n" % (k, gene_cov_dict[k]))
    with open(args.output+'.used_contigs.list', 'w') as outf:
        for contig_id in sorted(used_contigs):
            outf.write("%s\n" % contig_id)
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
