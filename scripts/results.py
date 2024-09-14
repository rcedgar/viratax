#!/usr/bin/python3

import sys
import vmrtt

print("SequenceID,Realm (-viria),Realm_score,Subrealm (-vira),Subrealm_score,Kingdom (-virae),Kingdom_score,Subkingdom (-virites),Subkingdom_score,Phylum (-viricota),Phylum_score,Subphylum (-viricotina),Subphylum_score,Class (-viricetes),Class_score,Subclass (-viricetidae),Subclass_score,Order (-virales),Order_score,Suborder (-virineae),Suborder_score,Family (-viridae),Family_score,Subfamily (-virinae),Subfamily_score,Genus (-virus),Genus_score,Subgenus (-virus),Subgenus_score,Species (binomial),Species_score")

vmrtt.load()

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

labels = set()
classified_labels = set()
for line in open("q_labels.txt"):
	labels.add(line.strip())

acc2node = {}
for line in open("../refdata/acc2node.tsv"):
	flds = line[:-1].split('\t')
	assert len(flds) == 3
	acc = flds[0]
	node = int(flds[1])
	acc2node[acc] = node

'''
					   0          1      2     3       4        5
ICTVTaxoChallenge_100057  NC_038726  100.0  0.98  species  1.0000
ICTVTaxoChallenge_100084   KT099124   54.1  0.96    other  0.5785
ICTVTaxoChallenge_100097   KJ010548  100.0  1.00  species  1.0000
'''

for line in open("top_nt_hit.tsv"):
	flds = line[:-1].split('\t')
	assert len(flds) == 6
	if flds[4] != "species":
		 continue
	label = flds[0]
	assert label in labels
	gbacc = flds[1]
	node = acc2node.get(gbacc)
	if node is None:
		continue
	confidence = float(flds[5])
	lineage = vmrtt.get_lineage(node)
	if lineage.get("Species") is None:
		sys.stderr.write("WARNING %s unclassified at species\n" % gbacc)
		continue

	assert not label in classified_labels
	classified_labels.add(label)
	s = label
	for rank in ranks:
		name = lineage[rank]
		if name is None:
			name = ""
			confstr = ""
		else:
			confstr = "1.0"
		s += "," + name
		if rank == "Species":
			s += ",%.2f" % confidence
		else:
			s += "," + confstr
	print(s)

'''
					   0        1    2                                         3
ICTVTaxoChallenge_763566  species  1.0                      Noxifervirus noxifer  MF063068/621,100  OL674541/256,83  KR>
 ICTVTaxoChallenge_69406  species  1.0  Heterocapsa circularisquama DNA virus 01    AB522601/14,92  AY261363/769,48  JF>
ICTVTaxoChallenge_359005  species  1.0                Nanovirus necroflaviviciae   GQ274023/10,100    KC979021/9,97    >
'''

for line in open("top_orf_hits_species.tsv"):
	flds = line[:-1].split('\t')
	assert len(flds) >= 4
	if flds[1] != "species":
		continue
	label = flds[0]
	assert label in labels
	confidence = float(flds[2])
	species_name = flds[3]
	idx = vmrtt.name2idx.get("Species:" + species_name)
	if idx is None:
		sys.stderr.write("WARNING missing from vmrtt: %s\n" % species_name)
		continue

	assert not label in classified_labels
	classified_labels.add(label)
	s = label
	for rank in ranks:
		name = lineage[rank]
		if name is None:
			name = ""
		s += "," + name
		if rank == "Species":
			s += ",%.2f" % confidence
		else:
			s += ",1.0"
	print(s)

'''
 ICTVTaxoChallenge_13245  Superkingdom:viruses=1.000  Realm:Duplodnaviria=1.000  Kingdom:Heunggongvirae=1.000  Phylum:U>
 ICTVTaxoChallenge_15650  Superkingdom:viruses=1.000  Realm:Duplodnaviria=1.000  Kingdom:Heunggongvirae=1.000  Phylum:U>
ICTVTaxoChallenge_125199  Superkingdom:viruses=1.000  Realm:Duplodnaviria=1.000  Kingdom:Heunggongvirae=1.000  Phylum:U>
'''

for line in open("classify_other.tsv"):
	flds = line[:-1].split('\t')
	label = flds[0]
	assert not label in classified_labels
	rank2taxon = {}
	rank2conf = {}
	for fld in flds[1:]:
		flds2 = fld.split(':')
		assert len(flds2) == 2
		rank = flds2[0]
		taxon_eq_conf = flds2[1]
		flds3 = taxon_eq_conf.split("=")
		assert len(flds3) == 2
		taxon = flds3[0]
		conf = float(flds3[1])
		rank2taxon[rank] = taxon
		rank2conf[rank] = conf

	classified_labels.add(label)
	s = label
	for rank in ranks:
		name = rank2taxon.get(rank)
		if name is None:
			name = ""
			confstr = ""
		else:
			confstr = "%.2f" % rank2conf[rank]
		s += "," + name
		s += "," + confstr
	print(s)

unclassified_labels = set()
for label in labels:
	if label in classified_labels:
		continue
	unclassified_labels.add(label)
	s = label
	for rank in ranks:
		s += ",,"
	print(s)

f = open("unclassified.txt", "w")
for label in sorted(list(unclassified_labels)):
	f.write(label + "\n")
f.close()
