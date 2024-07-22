#!/bin/bash

metadata=$1
prefix=$2

ensembl_ref=data/Homo_sapiens.GRCh38.cdna.ncrna.fa.gz
if [ ! -s $ensembl_ref ]; then
  wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
  wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
  mkdir -p data/
  cat Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz > $ensembl_ref
  rm Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz
fi

ensembl=data/ensembl_109_grch38.tsv

db=${prefix}_db.gem_df.pkl
if [ ! -s $db ]; then
  python3 src/build_training_db3.py $metadata $db $ensembl
fi

nonzero=0.9
model="${prefix}_model.${nonzero}.pkl"

if [ ! -s $model ]; then
  python3 src/train_pls_svm.py $db $metadata $model --nonzero $nonzero
fi
