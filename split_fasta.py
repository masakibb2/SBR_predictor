import sys
import base64
from core import seq2feature

MAX=1000

def split_fasta(fname,save_dir):
	for cnt,(idch,seq) in enumerate(seq2feature.parse_fasta(fname)):
		if cnt > MAX:
			#break
			pass
		wfname = base64.b64encode(idch.strip().replace("\t",""))
		with open("%s/%s.fasta" % (save_dir,wfname),"w") as fout:
			fout.write(">%s\n" % (idch.strip()))
			fout.write(seq + "\n")

def parse_arg():
	argvs = sys.argv
	argc = len(argvs)
	if argc == 3:
		return argvs[1],argvs[2]
	else:
		raise Exception

if __name__ == "__main__":
	fname,save_dir = parse_arg()
	split_fasta(fname,save_dir)
