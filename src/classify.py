import time
import pickle
import numpy as np
import argparse
from ensembl_tools import ensembl
from pls_svm import CompositeModel

def main(model_file, counts_file, ensembl_file):
  pls_model = pickle.load(open(model_file, 'rb'))

  # build vector with gene counts as per model, seq type {"nanopore", "illumina"}, and sample type {"ffpe", "other"}
  exp = np.zeros((1,len(pls_model.gene_names)+2))

  ensembl_table = ensembl(ensembl_file)
  for row in open(counts_file):
    #ENST00000644790.2	0.666667	6.550047
    fields = row.strip().split('\t')
    if fields[2] == "TPM":
      continue
    eid = ensembl_table.get_gene_id(fields[0]) # TPM
    gene_name = ensembl_table.get_gene_name(eid)
    if gene_name in pls_model.gene_names:
      exp[0,pls_model.gene_names.index(gene_name)] += float(fields[2])

  # assumed nanopore fresh(frozen) sample
  exp[0,-2] = 0
  exp[0,-1] = 1

  # normalize RPMM
  exp = exp * (1000000 / np.sum(exp))

  pred = pls_model.predict(exp)[0]

  print("name\tpredict_type\tpredict_subtype\t" + "\t".join([f"{t}" for t in pls_model.major_types]) + "\t" + "\t".join(["\t".join([f"{t}/{a}" for a in pls_model.subtype_names[t]]) for t in pls_model.major_types if t in pls_model.subtype_names]))
  print(f"{counts_file}\t{pred[0]}\t{pred[1]}"
      + "\t" + "\t".join([f"{a:.4f}" for a in pred[2][0]])
      + "\t" + "\t".join(["\t".join([f"{a:.4f}" for a in b[0]]) for b in pred[3]])
  )

if __name__ == "__main__":
  parser = argparse.ArgumentParser("Predict type from trained PLS-SVM model")
  parser.add_argument("model", help="Pickled trained composite PLS-SVM model")
  parser.add_argument("counts", help="Gene expression counts from (from minnow)")
  parser.add_argument("ensembl", help="Ensembl metadata file (transcript and gene names)")
  args = parser.parse_args()
  main(args.model, args.counts, args.ensembl)
