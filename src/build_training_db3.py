import pandas as pd
import sys
import os.path
import argparse
import numpy as np
import csv
from ensembl_tools import ensembl

def process_file(exp_file, ensembl_table, genes):
  lookup = {genes[i]:i for i in range(len(genes))}
  data = np.zeros(len(genes))
  missing_genes = 0
  if exp_file[exp_file.rindex('.')+1:] in ["txt", "cts", "sf"]: # probably minnow, Salmon, or TARGET data (with header)
    if not os.path.exists(exp_file):
      print(f"ERROR: {exp_file} does not exist! Skipping it")
      return None
    with open(exp_file) as ef:
      reader = csv.reader(ef, delimiter='\t')
      header = next(reader, None)
      #gene.expression (TARGET): gene	raw_counts	median_length_normalized	RPKM
      #minnow: transcript	reads	TPM
      #salmon: Name	Length	EffectiveLength	TPM	NumReads
      #gene.quantification.txt (TARGET/AML): gene	raw_count	rpkm
      for field_name in ["gene", "transcript", "Name"]:
        if field_name in header:
          gene_idx = header.index(field_name)
          break
      else:
        gene_idx = 0
      for field_name in ["RPKM", "TPM", "rpkm", "fpkm_normalized"]:
        if field_name in header:
          ct_idx = header.index(field_name)
          break
      else:
        ct_idx = 1
        # NO HEADER line, so reset the reader to start back at the beginning
        reader = csv.reader(ef, delimiter='\t')
      n_rows = 0
      skipped = 0
      for row in reader:
        n_rows += 1
        # check ensembl tables in order for a valid conversion
        if '|' in row[gene_idx]:
          for part in row[gene_idx].split('|')[::-1]: # usually something like "IL15|ENSG00000164136"
            eid = ensembl_table.get_gene_id(part)
            if eid is not None:
              break
        else:
          eid = ensembl_table.get_gene_id(row[gene_idx])
        if eid is not None:
          gene_name = ensembl_table.get_gene_name(eid)
          if gene_name in lookup:
            data[lookup[gene_name]] += float(row[ct_idx]) # multiple IDs or transcripts might contribute to the same gene, so add them
          else:
            missing_genes += 1
        else:
          skipped += 1
      if skipped > 0:
        print(f"WARNING: skipped {skipped} of {n_rows} ({skipped/n_rows*100:.2f}%) expression rows in {exp_file}")
  else:
    print("ERROR: unknown extension '{}' of file '{}'".format(td[td.rindex('.'):], td))
  if missing_genes > 0:
    print(f"WARNING: {missing_genes} unknown genes in {exp_file} (probably ENSEMBL version mismatch)")
  return data

def main(meta, out, ensembl_file):
  ensembl_table = ensembl(ensembl_file)

  samples = []
  with open(meta) as f:
    for parts in csv.reader(f, delimiter='\t'):
      if len(parts) < 10:
          print(parts)
      seq_id, exp_file, external_id, exclude, sample_type, rna_id, biosample_id, patient_id, lineage, subtype = parts
      if seq_id == "seq_id": # header
        continue
      if os.path.exists(exp_file):
        samples.append((seq_id, exp_file))
      else:
        print(f"ERROR: {exp_file} does not exist!")

  data = []

  for seq_id, exp_file in samples:
    print(seq_id, exp_file)
    data.append(process_file(exp_file, ensembl_table, ensembl_table.gene_names))

  df = pd.DataFrame(data, index=[s[0] for s in samples], columns=ensembl_table.gene_names)
  df.to_pickle(f"{out}")

if __name__ == "__main__":
  parser = argparse.ArgumentParser("Build gene expression matrix")
  parser.add_argument("meta", help="Sample metadata (TSV)")
  parser.add_argument("out", help="Output file for gene expression matrix")
  parser.add_argument("ensembl", help="Ensembl transcript table")
  args = parser.parse_args()
  main(args.meta, args.out, args.ensembl)
