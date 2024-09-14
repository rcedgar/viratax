#!/usr/bin/python3

import sys
import vmrtt

vmrtt.load()

acc2node = {}
for line in open("../refdata/acc2node.tsv"):
	flds = line[:-1].split('\t')
	assert len(flds) == 3
	acc = flds[0]
	node = int(flds[1])
	acc2node[acc] = node

'''
ICTVTaxoChallenge_763566  MF063068/621,100   OL674541/256,83  KR296692/1570,79  MK552327/3111,76   ON932084/271,76
 ICTVTaxoChallenge_69406    AB522601/14,92   AY261363/769,48   AF243438/567,48  JF999965/1961,48   MH509440/330,47
ICTVTaxoChallenge_359005   GQ274023/10,100     KC979021/9,97     KX431388/4,95     AB000923/4,82     KY070240/4,81
  ICTVTaxoChallenge_5526  HG962375/332,100   FN422399/809,98   KX171209/417,97   JX194238/407,97   KR052142/766,97
 '''
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

for line in open("top_orf_hits.tsv"):
	flds = line[:-1].split('\t')
	nr_hits = len(flds) - 1
	if nr_hits < 1:
		continue
	contig = flds[0]
	accs = []
	pctids = []
	nodes_eq100 = []
	nodes_ge90 = []
	nodes = []
	for hit in flds[1:]:
		flds2 = hit.split(',')
		assert len(flds2) == 2
		acc = flds2[0].split('/')[0]
		pctid = int(flds2[1])
		if pctid >= 90:
			node = acc2node.get(acc)
			if node is None:
				continue
			nodes_ge90.append(node)
			if pctid == 100:
				nodes_eq100.append(node)
	if len(nodes_eq100) == 0:
		species = consensus_species(nodes_ge90)
	else:
		species = consensus_species(nodes_eq100)
	genus = consensus_genus(nodes_ge90)
	species_ok = True
	if species is None:
		species_ok = False
	if len(nodes_ge90) > 0 and genus is None:
		species_ok = False

	confidence = 0
	if species_ok:
		if len(nodes_ge90) > 0 or len(nodes_eq100) > 1:
			confidence = 1
		else:
			confidence = 0.9

	s = contig
	if species_ok:
		s += "\tspecies"
		s += "\t%.1f" % confidence
		s += "\t" + species
	else:
		s += "\tother"
		s += "\t0.0"
		s += "\tNA"
	for hit in flds[1:]:
		s += "\t" + hit
	print(s)
