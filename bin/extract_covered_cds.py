#!/usr/bin/python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import pysam as ps
import gffutils
import sys
import argparse


def parse_sam(input_align):
    align = ps.AlignmentFile(input_align, "r")
    table=[]
    for read in align.fetch():
        if read.flag!=4:
            nodes=[]
            name=read.qname.split('_')[0]+'_'+read.qname.split('_')[1]
            nodes.append(name)
            nodes.append(read.reference_name)
            nodes.append(int(read.qlen))
            nodes.append(int(read.qstart))
            nodes.append(int(read.qend))
            nodes.append(int(read.reference_start))
            nodes.append(int(read.reference_end))
            nodes.append(read.flag)
            table.append(nodes)
    align.close()

    df=pd.DataFrame(data=table, columns=['Name','Node', 'Len', 'Qstart', 'Qstop', 'Rstart', 'Rstop', 'Flag'])
    df.sort_values(by=['Node','Rstart', 'Rstop'], inplace = True)
    #print(df.head())
    return df

def check (input_df, start, stop):

    light=False
    
    for i in range(1,len(input_df)):
        if light: 
            if input_df.Rstart.iloc[i] < input_df.Rstop.iloc[i-1] or stop > input_df.Rstart.iloc[i]:
                    loc.add(i-1)
                    loc.add(i)
                    if stop<=input_df.Rstop.iloc[i]:
                        length= (sum(input_df.Rstop.iloc[list(loc)].tolist()[:-1])+stop - sum(input_df.Rstart.iloc[list(loc)].tolist()[1:])+start)/(stop-start)*100
                        light=False
                        print(input_df.Name.iloc[list(loc)].tolist())
                        return length, input_df.Name.iloc[list(loc)].tolist(), 'part'
                    
        else: 
            loc=set()
            
        #skip genes 
        if start>input_df.Rstop.iloc[i-1] or stop<input_df.Rstart.iloc[i-1]:
            continue
        
        if start>=input_df.Rstart.iloc[i-1]:
            if stop<=input_df.Rstop.iloc[i-1]:
                return 100, input_df.Name.iloc[i-1], 'full'
            
            elif input_df.Rstart.iloc[i] < input_df.Rstop.iloc[i-1]: #or stop > input_df.Rstart.iloc[i]:
                    loc.add(i-1)
                    loc.add(i)
                    
                    if stop<=input_df.Rstop.iloc[i]:
                        return 100, input_df.Name.iloc[list(loc)].tolist(), 'full'
                    
            elif stop > input_df.Rstart.iloc[i]:
                    loc.add(i-1)
                    loc.add(i)
                    
                    if stop<=input_df.Rstop.iloc[i]:
                        length=(stop+input_df.Rstop.iloc[i-1] - start - input_df.Rstart.iloc[i])/(stop-start)*100
                        return length, input_df.Name.iloc[list(loc)].tolist(), 'part'
            
                    else: light=True
                    
        if start<input_df.Rstart.iloc[i-1]:
            if stop<=input_df.Rstop.iloc[i-1]:
                return (stop - input_df.Rstart.iloc[i-1])/(stop-start)*100, input_df.Name.iloc[i-1], 'part'
            
            elif input_df.Rstart.iloc[i] < input_df.Rstop.iloc[i-1] or stop > input_df.Rstart.iloc[i]:
                    loc.add(i-1)
                    loc.add(i)
                    
                    if stop>input_df.Rstop.iloc[i]:
                        return (max(input_df.Rstop.iloc[list(loc)].tolist()) - min(input_df.Rstart.iloc[list(loc)].tolist()))/(stop-start)*100, input_df.Name.iloc[list(loc)].tolist(), 'part'
       
        if start>input_df.Rstart.iloc[i-1] and stop>input_df.Rstop.iloc[i-1]:
            return (input_df.Rstop.iloc[i-1] - start)/(stop-start)*100, input_df.Name.iloc[i-1], 'part'
        
                
            
    return 0, '', ''

def check_one (input_df, start, stop):
    for i in range(0,len(input_df)):
        if start>input_df.Rstop.iloc[i] or stop<input_df.Rstart.iloc[i]:
            continue
        if start>input_df.Rstart.iloc[i] and stop>input_df.Rstop.iloc[i]:
            return (input_df.Rstop.iloc[i] - start)/(stop-start)*100, input_df.Name.iloc[i], 'part'
        if start<input_df.Rstart.iloc[i]:
            if stop<=input_df.Rstop.iloc[i]:
                return (stop - input_df.Rstart.iloc[i])/(stop-start)*100, input_df.Name.iloc[i], 'part'
        if start>=input_df.Rstart.iloc[i]:
            if stop<=input_df.Rstop.iloc[i]:
                return 100, input_df.Name.iloc[i], 'full'
    return 0, '', ''

def get_result (df,db):
    results = []
    for i in db.features_of_type('CDS'):
    
        gene = i
        start = gene.start
        stop = gene.stop
        c_name = gene.seqid
    
        subset = df[df.Node == c_name]
        if len(subset) > 1:
            cov, r_name, type_gene = check(subset, start, stop)

        else:
            cov, r_name, type_gene = check_one(subset, start, stop)
   
        if cov != 0:
            results.append([c_name, r_name, start, stop, cov, gene.id, gene.strand, type_gene])
        
    return results

def extract_genes(df, genes, new_genes):
    name_genes = df.Gene.tolist()
    genes = SeqIO.parse(genes, 'fasta')
    new_fasta = []
    for i in genes:
        if i.id in name_genes:
            new_fasta.append(i)
    SeqIO.write(new_fasta, new_genes + '.fasta', "fasta")


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="Output file name")
    parser.add_argument("--sam", "-f", type=str, help="initial file with alignment")
    parser.add_argument("--gff", "-g", type=str, help="initial file with annotation")
    parser.add_argument("--genes", "-r", type=str, help="initial file with covered genes")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    sam = args.sam
    gff = args.gff
    genes = args.genes
    new_genes = args.output
    gffutils.create_db(gff, 'db')
    db = gffutils.FeatureDB('db', keep_order=True)
    all_genes = get_result(parse_sam(sam), db)
    df_results=pd.DataFrame(data=all_genes, columns=['Node', 'Name', 'Start', 'Stop', 'Cov', 'Gene', 'Strand', 'Type'])
    df_results.to_csv(new_genes + '.csv')
    extract_genes(df_results, genes, new_genes)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)