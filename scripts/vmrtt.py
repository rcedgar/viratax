#!/usr/bin/python3

import sys

default_fn = "../refdata/taxtree.tsv"

rankset = set()
ranks = []
rank2idx= {}
names = []
name2idx = {}
name2rankidx = {}
name2parentname = {}
idx2parentidx = {}

def load(fn = default_fn):
	names.append("(root)")
	for line in open(fn):
		flds = line[:-1].split('\t')
		if flds[0] == "rank":
			assert len(flds) == 3
			idx = int(flds[1])
			rank = flds[2]
			assert idx == len(ranks)
			assert not rank in rankset
			rankset.add(rank)
			ranks.append(rank)
			rank2idx[rank] = idx
			continue
		assert len(flds) == 4
		idx = int(flds[0])
		name = flds[1]
		parentidx = int(flds[2])
		parentname = flds[3]
		assert idx == len(names)
		names.append(name)
		name2idx[name] = idx
		name2parentname[name] = parentname
		idx2parentidx[idx] = parentidx

def split_name(name):
	flds = name.split(':')
	assert len(flds) == 2
	rank = flds[0]
	taxon = flds[1]
	assert rank in rankset
	return rank, taxon

def get_lineage(node):
	lineage = {}
	for rank in ranks:
		lineage[rank] = None
	for sanity_counter in range(len(ranks)+1):
		name = names[node]
		if name == "(root)":
			return lineage
		rank, taxon = split_name(name)
		lineage[rank] = taxon
		parentidx = idx2parentidx.get(node)
		assert not parentidx is None
		node = parentidx
	assert False, "Loop in taxtree"

def get_lineagevec(node):
	lineagevec = []
	for rank in ranks:
		lineagevec.append(rank + ":_unclassified_")
	for sanity_counter in range(len(ranks)+1):
		name = names[node]
		if name == "(root)":
			return lineagevec
		rank, taxon = split_name(name)
		rankidx = rank2idx[rank]
		lineagevec[rankidx] = name
		parentidx = idx2parentidx.get(node)
		assert not parentidx is None
		node = parentidx
	assert False, "Loop in taxtree"

def get_common_lineage(lin1, lin2):
	common_lin = {}
	for rank in ranks:
		taxon1 = lin1[rank]
		taxon2 = lin2[rank]
		if taxon1 == taxon2:
			common_lin[rank] = taxon1
		else:
			common_lin[rank] = None
	return common_lin

def get_lctnode(node1, node2):
	lin1 = get_lineage(node1)
	lin2 = get_lineage(node2)
	lowest_common_rank = None
	lowest_common_name = None
	for rank in ranks:
		taxon1 = lin1[rank]
		taxon2 = lin2[rank]
		if taxon1 is None or taxon2 is None:
			continue
		if taxon1 == taxon2:
			lowest_common_rank = rank
			lowest_common_taxon = taxon1
	lowest_common_name = lowest_common_rank + ":" + lowest_common_taxon
	lctnode = name2idx[lowest_common_name]
	return lctnode

def get_rankidx_node(node):
	name = names[node]
	rank = name.split(":")[0]
	rankidx = rank2idx[rank]
	return rankidx

def get_genus(node):
	lineage = get_lineage(node)
	return lineage.get("Genus")

def get_species(node):
	lineage = get_lineage(node)
	return lineage.get("Species")
