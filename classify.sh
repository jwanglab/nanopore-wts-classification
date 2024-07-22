#!/bin/bash

fq=$1
model=$2

ensembl_ref=data/Homo_sapiens.GRCh38.cdna.ncrna.fa.gz
ensembl=data/ensembl_109_grch38.tsv

outdir=results
mkdir -p $outdir
if [[ "$fq" == "" ]]; then
  echo "usage: sh classify.sh <fq> <model>"
  exit
fi
if [ ! -s $fq ]; then
  echo "$fq does not exist!"
  exit
fi

filename=${fq##*/}
prefix=${filename%%.*}
paf=$outdir/${prefix}.hg38transcriptome.paf
cts=$outdir/${prefix}.ensembl.cts

if [ ! -s $paf ]; then
  minimap2 -cx map-ont -t 4 $ensembl_ref $fq > $paf
fi
if [ ! -s $cts ]; then
  minnow -p $paf > $cts
fi

stats=$(fastat -c $fq | tail -n1);
aligned_reads=$(awk '{tot+=$2}END{printf "%d", tot}' $cts);
unique_genes=$(awk '{if($2>0)ct++}END{printf "%d", ct}' $cts);
shannon_entropy=$(python3 src/tx_shannon_entropy.py $cts);
echo -e "seq_id\tfile\tn_seqs\ttotal_size\tavg_size\tmedian\tmaximum\tN50\taligned_reads\tunique_genes\tshannon_entropy" > $outdir/${prefix}.stats.tsv;
echo -e "${prefix}\t${stats}\t${aligned_reads}\t${unique_genes}\t${shannon_entropy}" >> $outdir/${prefix}.stats.tsv;

python3 src/classify.py $model $cts $ensembl > $outdir/${prefix}.classification.tsv
