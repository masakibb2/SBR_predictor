import glob
import os
import sys
import base64
#DIR="cluster1/pssm/*.psm.2"

def iter_pssm(fname):
	with open(fname) as fp:
		for rec in (line.strip().split() for line in iter(fp.readline,"")):
			if len(rec) == 44:
				yield rec[1],[int(i) for i in rec[2:22]]

if __name__ == "__main__":
	if len(sys.argv) < 2:
		sys.exit('Usage: %s pssm directory' % sys.argv[0])
	if not os.path.exists(sys.argv[1]):
		sys.exit('ERROR: Directory %s not found' % sys.argv[1])
	if not os.path.isdir(sys.argv[1]):
		sys.exit('ERROR: %s is not Directory' % sys.argv[1])

	DIR = sys.argv[1]
	if DIR[-1] != "/":
		DIR += "/"
	
	for fname in glob.iglob(DIR + "*.pssm"):
		# remove extension.
		print ">%s" % base64.b64decode(fname.split("/")[-1].replace(".fasta.pssm",""))
		for res,_pssm in iter_pssm(fname):
			print "\t".join([str(i) for i in _pssm])
