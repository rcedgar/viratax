#!/bin/bash -e

contigs=../contigs/ictv_challenge_contigs.fa

if [ ! -s ./RUNME.bash ] ; then
	x=`echo $PWD | sed "-es/.*viratax/scripts/@/"`
	if [ x != @ ] ; then
		echo "Must start ./RUNME.bash from viratax/scripts/ directory"
		exit 1
	fi
fi

PATH=$PATH:$PWD:$PWD/../bin

rm -rf ../tmp ../result
mkdir -p ../tmp ../result
cd ../tmp

grep "^>" $contigs \
  | sed "-es/^>//" \
  > q_labels.txt

getorf \
	-sequence $contigs \
	-outseq q.orfs.fa \
	-minsize 99 \
	-find 0

relabel_orfs.py q.orfs.fa \
	> q.orfs_relabel.fa

minimap2 \
	-t 24 \
	../refdata/gb.fa \
	$contigs \
	> minimap2.paf

top_nt_hit.py \
	> top_nt_hit.tsv

not_species_orfs.py \
	> not_species_orfs.fa

diamond blastp \
  --query not_species_orfs.fa \
  --db ../refdata/gb \
  --sensitive \
  --out diamond_hits.tsv \
  --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue pident

top_orf_hits.py \
  > top_orf_hits.tsv

top_orf_hits_species.py \
  > top_orf_hits_species.tsv

classify_other.py \
  > classify_other.tsv

result.py \
  > ../result/result.csv
