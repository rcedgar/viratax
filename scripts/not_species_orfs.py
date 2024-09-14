#!/usr/bin/python3

import sys
import fasta

orfs_fn = "q.orfs_relabel.fa"
hits_fn = "top_nt_hit.tsv"

contigs = set()
for line in open(hits_fn):
	flds = line[:-1].split('\t')
	assert len(flds) == 6
	so = flds[4]
	if so == "species":
		continue
	elif so == "other":
		contigs.add(flds[0])
	else:
		assert False

def on_seq(label, seq):
	contig = label.split('/')[0]
	if contig in contigs:
		fasta.WriteSeq(sys.stdout, seq, label)

fasta.ReadSeqsOnSeq(orfs_fn, on_seq)
