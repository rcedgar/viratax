![Muscle5](http://drive5.com/images/viratax_logo.jpg)

## Viratax pipline for viral metagenomics taxonomy classification

Viratax was developed for the ICTV Challenge (https://ictv-vbeg.github.io/ICTV-TaxonomyChallenge/).

## Pre-calculated results

See `results/results.csv`.

## Running the pipeline
To re-generate `results/results.csv`, choose one of the following methods.

### Using Linux after the repository is made public

Clone the repository and execute `./RUNME.bash` from the `scripts/` directory.

<pre>
git clone https://github.com/rcedgar/viratax.git
cd scripts
./RUNME.bash
</pre>

### Using Linux while the repository is private.

Download the repository from Serratus at https://serratus-public.s3.amazonaws.com/rce/viratax/viratax-repo-2024-09-14.tar.gz, extract the tarball and execute `./RUNME.bash` from the `scripts/` directory.

<pre>
mkdir viratax
cd viratax
wget https://serratus-public.s3.amazonaws.com/rce/viratax/viratax-repo-2024-09-14.tar.gz
tar -zxvf viratax-repo-2024-09-14.tar.gz
cd scripts
./RUNME.bash
</pre>

### Using Docker
There is a Dockerfile in the `docker/` directory. To build the container, run `./build_container.bash`. To run the container, use `./runit.bash`. The generated `results.csv` file will be written to the `docker/` directory.

<pre>
cd viratax/docker
./build_container.bash
./runit.bash
</pre>

### Dependencies

To reproduce the pre-computed results exactly, you need `minimap2` 2.24-r1122, `diamond` v2.0.14 and EMBOSS `getorf` v6.6.0.0. For convenience, these binaries are included in the `bin/` directory.

Also required are `python3` to run scripts (any v3 should work, I used v3.10.12) and `wget` to download data.

### Execution time and RAM

On a 32-core Intel i9-14900K Linux server, the entire pipeline from downloading data through final classification takes about 10 minutes. If you run on stand-alone Linux, the data will only be downloaded once. The second time `RUNME.bash` is executed it will perform classification using the reference data and query sequences previously downloaded. If you use Docker, the file system is not persistent so the query and reference data must be downloaded from scratch each time. RAM requirements are determined by `diamond` and `minimap2` which are memory-efficient by modern standards; how much is needed depends on how many threads you are running but at a guess 8Gb or 16Gb is probably more than enough.

### Reproducibility

The pipeline should give essentially identical results if run twice on the same data there are no randomized procedures like boostrapping. Each query sequence is processed independently, so the results for any given query is reproducible regardless of other contigs present.

### Updating reference data for a new VMR

The Challenge rules specify that ICTV release 39 should be used, so there is no provision here for updating the reference. The process of creating a reference is mostly automated, though some manual curation is needed to deal with a few minor issues such as idiosyncrasies in the VMR spreadsheet (`VMR_MSL39_v1` was used here). If there is interest in using a different VMR, please open an issue in the GitHub repository and I'll document the procedure.

### Design goals

Viratax was designed to be a fast, lightweight classifier which hopefully performs reasonably well on the Challenge. It emphasizes speed, simplicity and reproducibility over maximizing classification accuracy, which looks to be a very hard problem both in terms of implementation (phylogeny and trait prediction from contigs is often very difficult and clade-specific) and benchmarking (how can you verify if the classification of a highly diverged contig is correct?).

### How does it work?

The idea is to generalize rules of thumb relating sequence identity to taxonomic rank. For example, 90% aa or nt identity is often used as a species threshold. Similarly, 75% aa identity is roughly genus, and 50% aa identity is roughly family. These thresholds are gene- and taxon-specific; e.g. 90% aa identity for species is often applied to RNA virus RdRp while other genes diverge more quickly. To address this, every gene in every reference genome (in practice, every ORF) is aligned against all other reference genes to determine ORF-specific thresholds. ORFs in a query contig are aligned to ORFs from reference genomes, and a consensus is taken over aligned ORFs to determine the classification. The confidence estimate reflects identities of the top hits (higher identity gives higher confidence) and consistency (lower consensus gives lower confidence). A fast nucleotide search is performed first to identify hits which are consistent with a species match, i.e. >90% nt identity with >90% coverage of the contig.
