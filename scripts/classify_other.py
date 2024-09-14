#!/usr/bin/python3

import sys
import vmrtt

vmrtt.load()
ranks = vmrtt.ranks
nr_ranks = len(vmrtt.ranks)

acc2node = {}
for line in open("../refdata/acc2node.tsv"):
	flds = line[:-1].split('\t')
	assert len(flds) == 3
	acc = flds[0]
	node = int(flds[1])
	acc2node[acc] = node

gborf2lcts = {}
n = 0
for line in open("../refdata/gb_self_hits_lcts.tsv"):
	if n%10000 == 0:
		sys.stderr.write("%d reading self hits\r" % n)
	n += 1
	flds = line[:-1].split('\t')
	gborf = flds[0]
	lcts = []
	for fld in flds[1:]:
		flds2 = fld.split("=")
		assert len(flds2) == 2
		node = int(flds2[0])
		r = flds2[1]
		if r.find('-') > 0:
			flds3 = r.split('-')
			assert len(flds3) == 2
			lo = int(flds3[0])
			hi = int(flds3[1])
		else:
			lo = int(r)
			hi = lo
		lcts.append((node, lo, hi))
	gborf2lcts[gborf] = lcts
sys.stderr.write("%d reading self hits done\n" % n)

def consensus_species(nodes):
	if len(nodes) == 0:
		return None
	the_species = None
	for node in nodes:
		species = vmrtt.get_species(node)
		if not species is None:
			if the_species is None:
				the_species = species
			elif species != the_species:
				return None
	return the_species

def consensus_genus(nodes):
	if len(nodes) == 0:
		return None
	the_genus = None
	for node in nodes:
		genus = vmrtt.get_genus(node)
		if not genus is None:
			if the_genus is None:
				the_genus = genus
			elif genus != the_genus:
				return None
	return the_genus

def get_weight(pctid, lo, hi):
	if pctid > hi:
		w = 1
	elif pctid >= lo:
		w = pctid/hi
	else:
		w = (pctid/lo)*(lo/hi)*(lo/100)
	return w

def get_weight_acc(rank, pctid):
	if pctid == 100:
		return 1
	if pctid >= 90:
		if rank == "Species":
			return 0.5
		if rank == "Genus" or rank == "Subgenus":
			return 0.8
		return 1.0
	if pctid >= 70:
		if rank == "Species":
			return 0
		if rank == "Genus" or rank == "Subgenus":
			return 0.5
		return 0.8
	if pctid >= 50:
		if rank == "Species":
			return 0
		if rank == "Genus" or rank == "Subgenus":
			return 0.1
		if rank == "Family" or rank == "Subfamily":
			return 0.5
		return 0.8
	return 1

def classify(name2sumweight):
	best_name = None
	best_weight = 0
	confidence = 0
	sumsumweight = 0
	for name, sumweight in name2sumweight.items():
		if sumweight > best_weight:
			best_name = name
			best_weight = sumweight
		sumsumweight += sumweight
	if sumsumweight > 0:
		confidence = best_weight/sumsumweight
	return best_name, confidence

'''
                      0        1     2                                         3                 4                 5...
ICTVTaxoChallenge_763566  species  1.0                      Noxifervirus noxifer  MF063068/621,100   OL674541/256,83
 ICTVTaxoChallenge_69406  species  1.0  Heterocapsa circularisquama DNA virus 01    AB522601/14,92   AY261363/769,48
ICTVTaxoChallenge_359005  species  1.0                Nanovirus necroflaviviciae   GQ274023/10,100     KC979021/9,97
  ICTVTaxoChallenge_5526  species  1.0                          Litunavirus Ab09  HG962375/332,100   FN422399/809,98
 ICTVTaxoChallenge_13245    other  0.0                                        NA  MN175603/514,100  MH271321/232,100
 '''

for line in open("top_orf_hits_species.tsv"):
	flds = line[:-1].split('\t')
	if flds[1] == "species":
		continue
	assert flds[1] == "other"
	nr_hits = len(flds) - 4
	if nr_hits < 1:
		continue
	contig = flds[0]
	rankidx2name2sumweight = []
	nameset = set()
	for rankidx in range(nr_ranks):
		rankidx2name2sumweight.append({})
	toppctid = None
	any = False
	gbaccs = []
	gbacc_pctids = []
	for hit in flds[4:]:
		# MF063068/621,100
		flds = hit.split(',')
		assert len(flds) == 2
		gborf = flds[0]
		pctid = int(flds[1])
		gbacc = gborf.split('/')[0]
		if toppctid is None:
			toppctid = pctid
		elif pctid < toppctid - 10:
			break
		gbaccs.append(gbacc)
		gbacc_pctids.append(pctid)
		lcts = gborf2lcts.get(gborf)
		if lcts is None:
			continue
		any = True
		for node, lo, hi in lcts:
			weight = get_weight(pctid, lo, hi)
			lineagevec = vmrtt.get_lineagevec(node)
		lineagevec = vmrtt.get_lineagevec(node)
		for rankidx in range(nr_ranks):
			name = lineagevec[rankidx]
			if not name in nameset:
				nameset.add(name)
				rankidx2name2sumweight[rankidx][name] = 0
			rankidx2name2sumweight[rankidx][name] += weight
	consensus_lineagevec = [None]*nr_ranks
	s = contig
	if any:
		for rankidx in range(nr_ranks):
			name2sumweight = rankidx2name2sumweight[rankidx]
			if len(name2sumweight) == 0:
				continue
			best_name, confidence = classify(name2sumweight)
			if best_name is None or best_name.find("_unclassified_") > 0:
				continue
			s += "\t%s=%.3f" % (best_name, confidence)
		print(s)
	else:
		for i in range(len(gbaccs)):
			gbacc = gbaccs[i]
			node = acc2node.get(gbacc)
			if node is None:
				continue
			pctid = gbacc_pctids[i]
			lineagevec = vmrtt.get_lineagevec(node)
			for rankidx in range(nr_ranks):
				rank = ranks[rankidx]
				name = lineagevec[rankidx]
				weight = get_weight_acc(rank, pctid)
				if not name in nameset:
					nameset.add(name)
					rankidx2name2sumweight[rankidx][name] = 0
				rankidx2name2sumweight[rankidx][name] += weight
		for rankidx in range(nr_ranks):
			rank = ranks[rankidx]
			if rank == "Species":
				continue
			name2sumweight = rankidx2name2sumweight[rankidx]
			if len(name2sumweight) == 0:
				continue
			best_name, confidence = classify(name2sumweight)
			if best_name is None or best_name.find("_unclassified_") > 0:
				continue
			if rank == "Genus" and confidence > 0.5:
				confidence = 0.5
			s += "\t%s=%.3f" % (best_name, confidence)
		print(s)
