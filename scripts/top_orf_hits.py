#!/usr/bin/python3

import sys

TOPN = 10

acc2node = {}
for line in open("../refdata/acc2node.tsv"):
	flds = line[:-1].split('\t')
	assert len(flds) == 3
	acc = flds[0]
	node = int(flds[1])
	acc2node[acc] = node

def get_acc(orflabel):
	return orflabel.split('/')[0]

def tuple2(x):
	return x[1]

def flush():
	if current_contig is None:
		return
	acc2topid = {}
	acc2toporf = {}
	accset = set()
	for acc, gborf, pctid in hits:
		if not acc in accset:
			accset.add(acc)
			acc2topid[acc] = 0
		topid = acc2topid[acc]
		if pctid > topid:
			acc2topid[acc] = pctid
			acc2toporf[acc] = gborf
	hitlist = []
	for acc in accset:
		topid = acc2topid[acc]
		toporf = acc2toporf[acc]
		hitlist.append((acc, topid, toporf))
	sorted_hitlist = sorted(hitlist, key=tuple2, reverse=True)
	s = current_contig
	for acc, pctid, orf in sorted_hitlist[:TOPN]:
		s += "\t%s,%d" % (orf, pctid)
	print(s)

current_contig = None
hits = []
n = 0
for line in open("diamond_hits.tsv"):
	if n%10000 == 0:
		sys.stderr.write("%d finding top orf hits\r" % n)
	n += 1
	flds = line[:-1].split('\t')
	contigorf = flds[0]
	contig = get_acc(contigorf)
	if contig != current_contig:
		flush()
		current_contig = contig
		hits = []
	gborf = flds[2]
	pctid = int(round(float(flds[9])))
	gbacc = get_acc(gborf)
	if not node is None:
		hits.append((gbacc, gborf, pctid))
flush()

sys.stderr.write("%d finding top orf hits done\n" % n)
