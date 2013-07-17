import random
import util

def parse_fasta(fname):
	"""
	>>> {k:v for k,v in parse_fasta(".test.fasta")} == {'148l_E_1_164':'MNIFEMLRID','1ac5_A_1_483':'LPSSEEYKVA','1af6_A_1_421':'VDFHGYARSG'}
	True
	"""
	dic = {}
	with open(fname) as fp:
		for rec in (line.strip() for line in iter(fp.readline,"")):
			if len(rec) == 0:
				continue
			if rec[0] == '>':
				if len(dic) != 0:
					# for bigger window size.
					yield dic.items()[0]
				pdbid = rec[1:].strip()
				dic = {pdbid:''}
			else:
				dic[pdbid] += rec
		yield dic.items()[0]



def fasta2seq(fname):
	"""
	>>> [i for i in fasta2seq(".test.fasta")] == [("148lE",1,'MNIFEMLRID'),("1ac5A",1,'LPSSEEYKVA'),("1af6A",1,'VDFHGYARSG')]
	True
	"""
	# yield idch,start position, sequence
	for seq in parse_fasta(fname):
		#pdbid,ch,start,end = seq[0].split('_')
		rec = seq[0].split('_')
		#idch = pdbid + ch
		idch = rec[0] + rec[1]
		_seq = seq[1]
		#yield idch.replace('_',''),int(start),_seq
		yield idch.replace('_',''),int(rec[2]),_seq


def mkvec(ftr,seq,window,pos):
	"""
	>>> import feature
	>>> mkvec(lambda seq:feature.seq2frq(seq),"AAKDDECCSG",window=10,pos = 4) == {706: 1, 163: 1, 1222: 1, 9: 1, 3243: 1, 844: 1, 436: 1, 862: 1}
	True
	"""
	# ftr is the function make dict like "{pos:value for i in dimension}"
	n = seq[max(0,int(pos - window/2.0)):pos]
	c = seq[pos + 1:min(int(pos + window/2.0) + 1,len(seq))]
	vec = ftr('+'*(int(window/2 - len(n))) + n + seq[pos] + c + '-'*(int(window/2.0 - len(c))))
	return vec

def mkdtst_train(fname,window,ftr,answer):
	"""
	>>> import feature
	>>> [i for i in mkdtst_train('.test.fasta',10,lambda seq:feature.seq2frq(seq),answer='../dataset/answer_monod4.0.cluster1.txt')]
	None
	"""
	# only binding site For make positive dataet
	ans_sgr = util.ans(answer)
	# ans2int is function that return +1 if True else return -1
	ans2int = lambda idch,pos: 1 if ans_sgr.isans(idch,pos) else -1
	
	for idch,start,seq in fasta2seq(fname):
		idch = idch.strip()
		# 2012/1/31 pos -> pos + start
		# !!! now nodyfying !!!
		yield idch,{pos:(ans2int(idch,start + pos),mkvec(ftr,seq,window,pos))
					for pos in range(len(seq)) if ans_sgr.isans(idch,start + pos)}

def mkdtst_test(fname,window,ftr,answer):
	# For test
	"""
	>>> import feature
	>>> [i for i in mkdtst_test('.test.fasta',10,lambda seq:feature.seq2frq(seq),answer='../dataset/answer_monod4.0.cluster1.txt')]
	None
	"""
	
	ans_sgr = util.ans(answer)
	ans2int = lambda idch,pos: 1 if ans_sgr.isans(idch,pos) else -1
	
	for idch,start,seq in fasta2seq(fname):
		idch = idch.strip()
		# 2012/1/31 pos -> pos + start
		# !!! now nodyfying !!!
		yield idch,{pos:(ans2int(idch,start + pos),mkvec(ftr,seq,window,pos))
					for pos in range(len(seq))}


def mkneg_train(fname,window,ftr,size = 5):
	"""
	>>> import feature
	>>> [i for i in mkneg_train('.test.fasta',10,lambda seq:feature.seq2frq(seq))]
	None
	"""
	# To Dataset randomly choiced positions
	for negid,seq in parse_fasta(fname):
		if seq.find('X') > 0:
			continue
		_indx = [i for i in range(1,len(seq))]
		random.shuffle(_indx)
		yield negid,{pos:(-1,mkvec(ftr,seq,window,pos)) for pos in sorted(_indx[:size])}


def mkneg_test(fname,window,ftr):
	# For test
	"""
	>>> import feature
	>>> [i for i in mkneg_test('.test.fasta',10,lambda seq:feature.seq2frq(seq))]
	None
	"""
	# to Dataset whole seqence
	for negid,seq in parse_fasta(fname):
		if seq.find('X') > 0:
			continue
		yield negid,{pos:(-1,mkvec(ftr,seq,window,pos)) for pos in range(len(seq))}

def mkdtst_near(fname,window,ftr,answer,low,up):
	# For test
	"""
	>>> [i for i in mkdtst_test('.test.pssm',10,answer='./.answer.test')]
	None
	"""
	ans_sgr = util.ans(answer)
	ans2int = lambda idch,pos: 1 if ans_sgr.isans(idch,pos) else -1
	
	for idch,start,seq in fasta2seq(fname):
		idch = idch.strip()
		# 2012/1/31 pos -> pos + start
		# !!! now nodyfying !!!
		yield idch,{pos:(ans2int(idch,start + pos),mkvec(ftr,seq,window,pos))
					for pos in range(len(seq)) if low <= ans_sgr.get_dist(pos + start,idch) <= up
					or ans_sgr.isans(idch,start + pos)}

class dataset(object):
	# Interface between dataset and libsvm
	
	def __init__(self,pfasta,nfasta,answer,ftr,window):
		# pfasta is file name of positive dataset (Fasta format)
		self._pfasta = pfasta
		# nfasta is file name of negative dataset (Fasta format)
		self._nfasta = nfasta
		# answer is description of position of binding residue
		self._ans = answer

		# list of id in postive dataset
		if pfasta is not None:
			self.pids = [idch for idch,start,seq in fasta2seq(pfasta)]
		else:
			self.pids = None
		# list of id in negative dataset
		if nfasta is not None:
			self.nids = [rec[0] for rec in parse_fasta(nfasta)]
		else:
			self.nids = None

		self._ftr = ftr
		self._window = window

	def mkTest(self,part_pids,part_nids):
		"""
		>>> import feature
		>>> d = dataset(".test.fasta",".test.fasta","..//dataset/answer_monod4.0.cluster1.txt",lambda seq:feature.seq2frq(seq,2),10)
		>>> d.mkTest(['148lE'],['1af6_A_1_421'])
		None
		"""
		pdatasets = {idch:pvec for idch,pvec in mkdtst_test(self._pfasta,self._window,self._ftr,answer = self._ans) if idch in part_pids}
		ndatasets = {idch:pvec for idch,pvec in mkneg_test(self._nfasta,self._window,self._ftr) if idch in part_nids}
		pdatasets.update(ndatasets)

		return pdatasets.items()
	
	def mkTrain(self,part_pids,part_nids,size = 5):
		"""
		>>> import feature
		>>> d = dataset(".test.fasta",".test.fasta","..//dataset/answer_monod4.0.cluster1.txt",lambda seq:feature.seq2frq(seq,2),10)
		>>> d.mkTrain(['148lE'],['1af6_A_1_421'])
		None
		"""
		pdatasets = {idch:pvec for idch,pvec in mkdtst_train(self._pfasta,self._window,self._ftr,answer = self._ans) if idch in part_pids}
		ndatasets = {idch:pvec for idch,pvec in mkneg_train(self._nfasta,self._window,self._ftr,size = size) if idch in part_nids}
		pdatasets.update(ndatasets)

		return pdatasets.items()


class dataset2(object):
	def __init__(self,seq,answer,ftr,window,low = 5,up = 25):
		# pdata is data source of positive dataset
		self._pseq = seq
		# answer is description of position of binding residue
		self._ans = answer
		self._low = low
		self._up = up

		# list of id in postive dataset
		self.pids = [idch for idch,start,seq in fasta2seq(seq)]
		self._window = window
		self._ftr = ftr
	
	def mkTest(self,part_ids):
		pdatasets = {idch:pvec for idch,pvec in mkdtst_test(self._pseq,self._window,ftr = self._ftr,answer = self._ans) if idch in part_ids}
		return pdatasets.items()
	
	def mkTrain(self,part_ids):
		# negative dataset are generate if low <= distance <= up from binding residue in dataset.
		"""
		>>> d = dataset(".test.seq",".test.pssm",".answer.test",10)
		>>> d.mkTrain(['148lE'],['1af6_A_1_421'])
		None
		"""
		pdatasets = {idch:pvec for idch,pvec in mkdtst_near(self._pseq,self._window,ftr = self._ftr,answer = self._ans,low = self._low,up = self._up) if idch in part_ids}
		return pdatasets.items()

class dataset3(dataset2):
	# For Test dataset whose answer position is Unkonown
	def __init__(self,seq,window):
		# pdata is data source of positive dataset
		self._pseq = seq
		# list of id in postive dataset

		#self.pids = [idch for idch,start,seq in fasta2seq(seq)]
		self.pids = [idch.strip(".fasta.ckp") for idch,seq in parse_fasta(seq)]
		self._window = window

	def mkTest(self,part_ids):
		pdatasets = {idch:pvec for idch,pvec in mkneg_test(self._pseq,self._window,self._ftr) if idch in part_ids}
		return pdatasets.items()
	
if __name__ == "__main__":
	import doctest
	doctest.testmod()
