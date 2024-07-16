# ignore ALL sklearn warnings
def warn(*args, **kwargs):
  pass
import warnings
warnings.warn = warn

import time
import pandas as pd
import pickle
import numpy as np
import csv
import argparse
import pls_svm

def main(metadata_file, db_file, out, top_gene_list_file, nonzero_frac):
  lineage_field = "lineage"
  subtype_field = "subtype"
  person_id_field = "patient_id"
  sample_id_field = "seq_id"
  sample_type_field = "sample_type"

  # ------ load metadata ------
  t0 = time.time()
  metadata = {}
  with open(metadata_file) as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader, None)
    for row in reader:
      # assume the first column is the sample name, regardless of the header
      metadata[row[0]] = {header[i]:row[i] for i in range(len(header))}
  df = pd.read_pickle(db_file)
  run_names = list(df.index)
  gene_names = list(df.columns)
  matrix = df.to_numpy()

  t1 = time.time()
  print(f"{t1-t0:.4f} - reading matrix and metadata")
  drop = []
  for i in range(len(run_names)):
    if run_names[i] not in metadata:
      print(f"WARNING: {run_names[i]} missing from metadata, will be dropped!!!")
      drop.append(run_names[i])
      continue

  for a in drop:
    run_names.remove(a)
    df = df.drop(a)

  # ------ load and format gexp matrix ------
  disease_names = list(set([metadata[n][lineage_field] for n in run_names]))
  seq_names = list(set(["nanopore" if "Nano" in metadata[n][sample_id_field] else "illumina" for n in run_names]))
  sample_types = list(set(["ffpe" if "FFPE" in metadata[n][sample_type_field] else "other" for n in run_names]))

  # ------ FILTERS ------

  keep_indices = [idx for idx in list(range(len(run_names))) if metadata[run_names[idx]][lineage_field] != "unknown"]

  matrix = matrix[keep_indices,:]
  orig_names = run_names
  run_names = [run_names[i] for i in keep_indices]
  print(f"Kept {len(run_names)} of {len(orig_names)}")

  if top_gene_list_file is not None:
    with open(top_gene_list_file) as f:
      use_genes = [line.strip().split()[0] for line in f]
      if top_n_genes is None:
        top_n_genes = len(use_genes)
      df = df[use_genes[-top_n_genes:]]
    #print(f"Kept {matrix.shape[1]} genes from top genes list")
  else:
    # all genes must be >0 in EVERY sample...
    #matrix = matrix[:,np.count_nonzero(matrix, 0) == matrix.shape[0]]
    keep_features = np.count_nonzero(matrix, 0) > (matrix.shape[0] * nonzero_frac)
    matrix = matrix[:,keep_features]
    use_genes = [gene_names[k] for k in range(keep_features.shape[0]) if keep_features[k]]
    print(f"Kept {matrix.shape[1]} features (genes) nonzero in >{nonzero_frac*100:.1f}% of samples\n")

  seq_type = np.array([seq_names.index("nanopore" if "Nano" in metadata[n][sample_id_field] else "illumina") for n in run_names])
  # add material type as flag where 0 is FFPE (nanopore) and 1 is everything else (including all Illumina, where we don't know...)
  sample_type = np.array([sample_types.index("ffpe" if "FFPE" in metadata[n][sample_type_field] else "other") for n in run_names])
  matrix_plus_seq = np.zeros((matrix.shape[0], matrix.shape[1]+2), matrix.dtype)
  matrix_plus_seq[:,:-2] =  matrix
  matrix_plus_seq[:,-2] = seq_type
  matrix_plus_seq[:,-1] = sample_type
  matrix = matrix_plus_seq

  # renormalize RPKMs
  for row in range(matrix.shape[0]):
    matrix[row,:] = matrix[row,:] * (1000000 / np.sum(matrix[row,:]))

  t2 = time.time()
  print(f"{t2-t1:.4f} - filtering, normalizing, and otherwise futzing with matrix")

  ans = {}

  m = pls_svm.CompositeModel(matrix, [metadata[n][lineage_field] for n in run_names], [metadata[n][subtype_field] for n in run_names], use_genes)
  pickle.dump(m, open(f"{out}", 'wb'))
  t3 = time.time()
  print(f"{t3-t2:.4f} - training and saving composite model")
  print(f"{t3-t0:.4f} - done")


if __name__ == "__main__":
  parser = argparse.ArgumentParser("Train PLS/SVM classification model on gene expression data")
  parser.add_argument("matrix", help="Gene expression matrix")
  parser.add_argument("meta", help="Metadata")
  parser.add_argument("out", help="Output path/prefix")
  parser.add_argument("--use_genes", help="File with a list of gene names to use for classification")
  parser.add_argument("--nonzero", help="Proportion of samples in which a gene must be >0 to include it", type=float, default=0.99)
  args = parser.parse_args()
  main(args.meta, args.matrix, args.out, args.use_genes, args.nonzero)
