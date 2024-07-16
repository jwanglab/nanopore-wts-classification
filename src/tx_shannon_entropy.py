import skbio.diversity.alpha as adiv
import sys

d = [float(line.strip().split('\t')[2]) for line in open(sys.argv[1]) if not line.startswith("transcript")]
print(adiv.shannon(d))
