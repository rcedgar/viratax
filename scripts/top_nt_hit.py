#!/usr/bin/python3

import sys

# https://github.com/lh3/miniasm/blob/master/PAF.md

# 1 	string 	Query sequence name
# 2 	int 	Query sequence length
# 3 	int 	Query start (0-based; BED-like; closed)
# 4 	int 	Query end (0-based; BED-like; open)
# 5 	char 	Relative strand: "+" or "-"
# 6 	string 	Target sequence name
# 7 	int 	Target sequence length
# 8 	int 	Target start on original strand (0-based)
# 9 	int 	Target end on original strand (0-based)
# 10 	int 	Number of residue matches
# 11 	int 	Alignment block length
# 12 	int 	Mapping quality (0-255; 255 for missing)

MIN_PCTID = 90
MIN_QUERY_COV = 0.9

def confidence(pctid, qcov):
	conf = qcov*pctid/90
	if conf > 1:
		conf = 1
	return conf

def orfq(label):
	label = label.split()[0]
	n = label.rfind('_')
	return label[:n]

def is_species(pctid, qcov):
	return pctid >= MIN_PCTID and qcov >= MIN_QUERY_COV

def better_hit(pctid, qcov, top_pctid, top_qcov):
	return qcov > top_qcov and pctid - top_pctid >= 0

hits_fn = "minimap2.paf"

qs = set()
q2tophit = {}

species_set = set()
for line in open(hits_fn):
	flds = line[:-1].split('\t')
	assert len(flds) >= 12
	q = flds[0]
	ql = int(flds[1])
	qlo = int(flds[2])
	qhi = int(flds[3])
	qsegl = qhi - qlo
	assert qsegl <= ql

	t = flds[5]
	tl = int(flds[6])
	tlo = int(flds[7])
	thi = int(flds[8])
	if tlo > thi:
		tmp = tlo
		tlo = thi
		thi = tmp
	tsegl = thi - tlo
	assert tsegl <= tl
	
	nrids = int(flds[9])
	nrcols = int(flds[10])
	assert nrids <= nrcols

	pctid = nrids*100/nrcols

	qcov = qsegl/ql
	tcov = tsegl/tl

	if q not in qs:
		qs.add(q)
		q2tophit[q] = (t, pctid, qcov)

	top_t, top_pctid, top_qcov = q2tophit[q]
	if better_hit(pctid, qcov, top_pctid, top_qcov):
		q2tophit[q] = (t, pctid, qcov)

qlist = list(qs)
qlist = sorted(qlist)
not_species_set = set()
nr_species = 0
nr_other = 0
for q in qlist:
	t, top_pctid, top_qcov = q2tophit[q]
	conf = confidence(top_pctid, top_qcov)
	tacc = t.split(";")[0]
	s = q
	s += "\t" + tacc
	s += "\t%.1f" % top_pctid
	s += "\t%.2f" % top_qcov
	if is_species(top_pctid, top_qcov):
		s += "\tspecies"
		nr_species += 1
	else:
		s += "\tother"
		not_species_set.add(q)
		nr_other += 1
	s += "\t%.4f" % conf
	print(s)

sys.stderr.write("nt hits gave %d species, %d other\n" % (nr_species, nr_other))
