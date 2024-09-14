#!/usr/bin/python3

import sys

f = open("../results/results.csv")
hdrflds = f.readline()[:-1].split(',')
assert len(hdrflds) == 31

MIN_CONF = 0.8

ranks = [ "Realm",
	"Subrealm",
	"Kingdom",
	"Subkingdom",
	"Phylum",
	"Subphylum",
	"Class",
	"Subclass",
	"Order",
	"Suborder",
	"Family",
	"Subfamily",
	"Genus",
	"Subgenus",
	"Species" ]
nr_ranks = len(ranks)
assert 1 + 2*nr_ranks == 31

counts = [0]*nr_ranks
nr_unclassified = 0

for line in f:
	flds = line[:-1].split(",")
	contig = flds[0]
	lowest_rankidx = None
	assert len(flds) == 31
	for rankidx in range(nr_ranks):
		taxon = flds[1 + 2*rankidx]
		if taxon == "":
			continue
		conf = float(flds[2 + 2*rankidx])
		if conf >= MIN_CONF:
			lowest_rankidx = rankidx
	if lowest_rankidx is None:
		print(contig + "\tunclassified")
	else:
		print(contig + "\t" + ranks[lowest_rankidx])
