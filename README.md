Nanopore whole-transcriptome cancer classification
==================================================

This respository includes code for building and training a classification model and running predictions from nanopore sequencing-based whole-transcriptome profiles. We also include a pre-trained model for pediatric acute leukemias. 

Installation
------------

Tools to build the gene expression database and perform classification are written in Python (3.6+)

Dependencies:

    numpy
    pandas
    sklearn
    skbio

You will also need to download the ENSEMBL cDNA reference and associated metadata:

Ex. for release 109 (matching the trained pediatric acute leukemia model):

    wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
    wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
    cat Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz > Homo_sapiens.GRCh38.cdna.ncrna.fa.gz

An appropriate ENSEMBL metadata file can be generated by [biomart](https://m.ensembl.org/biomart/martview/) data dump with the following fields in CSV format:

    Gene stable ID
    Gene stable ID version
    Transcript stable ID
    Transcript stable ID version
    Gene name
    Gene Synonym
    Chromosome/scaffold name
    Gene start (bp)
    Gene end (bp)
    Gene description
    Transcript start (bp)
    Transcript end (bp)
    Transcript name
    HGNC symbol

Other dependencies - make sure they are available in your PATH:

  * [minimap2](https://github.com/attractivechaos/minimap2)
  * [minnow](https://github.com/jwanglab/minnow)
  * [fastat](https://github.com/jwanglab/fastat)

Getting started
---------------

Run the example data using the pre-trained model:

    for f in example/*; do
      sh classify.sh $f models/model_240501.default_e109.0.9.Leukemia.pkl
    done

This should take about 5 minutes total and will create several output files in **results/**. Classification results are in files __results/*.classification.tsv__ and should match:

| name | predict\_type | predict\_subtype | B-ALL | AML | TLL | B-ALL/ETV6-RUNX1 | B-ALL/Hyperdiploid / Near haploid | B-ALL/KMT2Ar | B-ALL/Low hypodiploid | B-ALL/Other | B-ALL/Ph / Ph-like | B-ALL/TCF3-PBX1 | AML/Core binding factor | AML/KMT2Ar | AML/Other |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |
| results/NALM6\_1m.ensembl.cts | B-ALL | Other | 0.9944 | 0.0029 | 0.0027 | 0.0023 | 0.0013 | 0.0024 | 0.0008 | 0.9879 | 0.0028 | 0.0025 | 0.0136 | 0.0096 | 0.9769 |
| results/REH\_1m.ensembl.cts | B-ALL | ETV6-RUNX1 | 0.9511 | 0.0077 | 0.0411 | 0.5748 | 0.0128 | 0.0206 | 0.0092 | 0.3228 | 0.0202 | 0.0395 | 0.0986 | 0.0586 | 0.8428 |
| results/SUPB15.ensembl.cts | B-ALL | Ph / Ph-like | 0.9950 | 0.0005 | 0.0045 | 0.0200 | 0.0258 | 0.0816 | 0.0086 | 0.0575 | 0.7469 | 0.0596 | 0.3759 | 0.1353 | 0.4888 |


Training and classifying with a model from scratch:

    sh train.sh metadata.tsv "leukemia_$(date +%Y%m%d)"
    sh classify.sh test.fq "leukemia_$(date +%Y%m%d)_model.0.9.pkl"

Usage
-----

Generate transcript counts from nanopore sequence using [minnow](https://github.com/jwanglab/minnow):

    Usage: minnow [options]
    Options:
      -q: FASTA/Q[.gz] file with reads
      -r: Reference FASTA/Q[.gz]
      -t: Threads (default: 1)
      -c: 'careful': more accurate but slower
      -v, --verbose: verbose
      -h, --help: show this
      --version: show version information

For example:

    minnow -q nanopore_reads.fastq.gz -r Homo_sapiens.GRCh38.cdna.ncrna.fa.gz > nanopore_transcripts.cts

Generate an appropriate metadata file (metadata.tsv) in the following format:

    seq_id	exp_file	sample_type	patient_id	lineage	subtype
    UNIQUE_Sample_ID	FILEPATH/minnow.cts	FROZ	0001	B-ALL	Hyperdiploid
    ...

Build gene expression matrix for training and cross-validation:

    usage: python src/build_training_db3.py [-h] meta out ensembl

    positional arguments:
      meta        Sample metadata (TSV)
      out         Output file for gene expression matrix
      ensembl     Ensembl transcript table

    options:
      -h, --help  show this help message and exit

ex.
    python src/build_training_db3.py metadata.tsv db.pkl ensembl_109_grch38.tsv

Train PLS/SVM classification model on gene expression data:

    usage: python src/train_pls_svm.py [-h] [--use_genes USE_GENES] [--nonzero NONZERO] matrix meta out

    positional arguments:
      matrix                Gene expression matrix
      meta                  Metadata
      out                   Output path/prefix

    options:
      -h, --help            show this help message and exit
      --use_genes USE_GENES
                            File with a list of gene names to use for classification
      --nonzero NONZERO     Proportion of samples in which a gene must be >0 to include it
  
ex.
    python src/train_pls_svm.py db.pkl metadata.tsv model.pkl

Predict/classify samples:

    usage: python src/classify.py [-h] model counts ensembl

    positional arguments:
      model       Pickled trained composite PLS-SVM model
      counts      Gene expression counts from (from minnow)
      ensembl     Ensembl metadata file (transcript and gene names)

    options:
      -h, --help  show this help message and exit

ex.
    python src/classify.py model.pkl nanopore_transcripts.cts ensembl_109_grch38.tsv > classification.tsv

Citation
--------

[Wang J, Bhakta N, Ayer Miller V, Revsine M, Litzow MR, Paietta E, Fedoriw Y, Roberts KG, Gu Z, Mullighan CG, Jones CD, Alexander TB. *Acute Leukemia Classification Using Transcriptional Profiles From Low-Cost Nanopore mRNA Sequencing*. JCO Precision Oncology. 2022 Apr;6:e2100326. doi: 10.1200/PO.21.00326](https://www.ncbi.nlm.nih.gov/pubmed/35442720/)