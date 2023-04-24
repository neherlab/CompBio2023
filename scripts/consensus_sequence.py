import argparse
import pandas as pd
import numpy as np
from create_allele_counts import nuc_alpha
from Bio import SeqIO, SeqRecord, Seq

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=str)
    parser.add_argument("--min-coverage", type=int, default=10)
    parser.add_argument("--seq-name", type=str)
    parser.add_argument("--output", type=str)
    args=parser.parse_args()

    ac = np.array(pd.read_csv(args.counts, sep='\t'))
    cov = ac[:, :5].sum(axis=1)
    seq = nuc_alpha[ac.argmax(axis=1)]
    seq[cov<args.min_coverage]='N'
    seq = np.array(seq, dtype='U1')
    SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(''.join(seq)), id=args.seq_name), args.output, 'fasta')
