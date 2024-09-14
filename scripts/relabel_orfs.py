#!/usr/bin/python3

import sys
import fasta

fn = sys.argv[1]

accs = set()
acc2count = {}

def on_seq(label, seq):
	label = label.split()[0]
	if label.find(";taxid=") > 0:
		flds = label.split(";")
		acc = flds[0]
	else:
		n = label.rfind("_")
		assert n > 0
		acc = label[:n]
	if not acc in accs:
		acc2count[acc] = 0
		accs.add(acc)
	n = acc2count[acc]
	n += 1
	acc2count[acc] = n

	fasta.WriteSeq(sys.stdout, seq, acc + "/" + str(n))

fasta.ReadSeqsOnSeq(fn, on_seq)
