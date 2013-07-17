import random
import util

class pssm(object):
	def __init__(self,pdbid):
		self.pdbid = pdbid
		self.pssm = []
	def append(self,vec):
		if len(vec) != 20:
			raise Exception,"pssm.append(): Argument length must be 20. len(vec) == %s" % len(vec)
		else:
			self.pssm.append(vec)
			

def parse_pssm(fname):
	
	# parse pssm
	"""
	>>> {k:v for k,v in parse_pssm(".test.pssm")}
	True
	"""
	prot_pssm = None
	with open(fname) as fp:
		for rec in (line.strip() for line in iter(fp.readline,"")):
			if len(rec) == 0:
				continue
			if rec[0] == '>':
				if prot_pssm is not None:
					yield prot_pssm.pdbid,prot_pssm.pssm
				prot_pssm = pssm(rec[1:].strip())
			else:
				try:
					prot_pssm.append([float(i) for i in rec.split()])
				except Exception:
					print rec
					print prot_pssm.idch
					exit()
					
		yield prot_pssm.pdbid,prot_pssm.pssm

def parse_pssm4pos(fname):
	"""
	>>> [i for i in parse_pssm4pos(".test.pssm")]
	True
	"""
	# yield idch,start position, sequence
	for prot_id,pssm in parse_pssm(fname):
		pdbid,ch,start,end = prot_id.split('_')
		idch = pdbid + ch
		yield idch.replace('_',''),int(start),pssm

def _paste(vec):
	out = []
	for i in vec:
		out+=i
	return out

def mkvec(pssm,window,pos):
	"""
	>>> _pssm = [i for i in parse_pssm4pos(".test.pssm")][0][2]
	>>> mkvec(_pssm,10,3)
	None
	"""
	nu = [0]*20
	n = pssm[max(0,int(pos - window/2.0)):pos]
	c = pssm[pos + 1:min(int(pos + window/2.0) + 1,len(pssm))]
	if window%2 == 1:
		nu_n = [nu for i in range(int(window/2 - len(n) + 1))]
	else:
		nu_n = [nu for i in range(int(window/2 - len(n)))]
	nu_c = [nu for i in range(int(window/2.0 - len(c)))]
	vec = nu_n + n + [pssm[pos]] + c + nu_c
	return _paste(vec)

def mkdtst_train(fname,window,answer):
	"""
	>>> import feature
	>>> [i for i in mkdtst_train('.test.pssm',10,answer='.answer.test')]
	None
	"""
	# only binding site For make positive dataet
	ans_sgr = util.ans(answer)
	# ans2int is function that return +1 if True else return -1
	ans2int = lambda idch,pos: 1 if ans_sgr.isans(idch,pos) else -1
	
	for idch,start,pssm in parse_pssm4pos(fname):
		idch = idch.strip()
		# 2012/1/31 pos -> pos + start
		# !!! now nodyfying !!!
		yield idch,{pos:(ans2int(idch,start + pos),mkvec(pssm,window,pos))
					for pos in range(len(pssm)) if ans_sgr.isans(idch,start + pos)}

def mkdtst_test(fname,window,answer):
	# For test
	"""
	>>> [i for i in mkdtst_test('.test.pssm',10,answer='./.answer.test')]
	None
	"""
	if answer is not None:
		ans_sgr = util.ans(answer)
		ans2int = lambda idch,pos: 1 if ans_sgr.isans(idch,pos) else -1
		for idch,start,pssm in parse_pssm4pos(fname):
			idch = idch.strip()
			# 2012/1/31 pos -> pos + start
			# !!! now nodyfying !!!
			yield idch,{pos:(ans2int(idch,start + pos),mkvec(pssm,window,pos))
						for pos in range(len(pssm))}
	else:
		# For negative dataset
		for idch,pssm in parse_pssm(fname):
			idch = idch.strip()
			yield idch,{pos:(-1,mkvec(pssm,window,pos)) for pos in range(len(pssm))}

def mkneg_train(fname,window,size = 8):
	"""
	>>> import feature
	>>> [i for i in mkneg_train('.test.pssm',10)]
	None
	"""
	# To Dataset randomly choiced positions
	for negid,pssm in parse_pssm(fname):
		_indx = [i for i in range(1,len(pssm))]
		random.shuffle(_indx)
		yield negid,{pos:(-1,mkvec(pssm,window,pos)) for pos in sorted(_indx[:size])}


def mkneg_train2(fname,window,size = 8):
	# For one class SVM
	"""
	>>> import feature
	>>> [i for i in mkneg_train('.test.pssm',10)]
	None
	"""
	# To Dataset randomly choiced positions
	for negid,pssm in parse_pssm(fname):
		_indx = [i for i in range(1,len(pssm))]
		random.shuffle(_indx)
		yield negid,[mkvec(pssm,window,pos) for pos in sorted(_indx[:size])]


def mkneg_test(fname,window):
	# For test
	"""
	>>> import feature
	>>> [i for i in mkneg_test('.test.pssm',10)]
	None
	"""
	# to Dataset whole seqence
	for negid,pssm in parse_pssm(fname):
		yield negid,{pos:(-1,mkvec(pssm,window,pos)) for pos in range(len(pssm))}

def mkdtst_near(fname,window,answer,low,up):
	# For test
	"""
	>>> [i for i in mkdtst_test('.test.pssm',10,answer='./.answer.test')]
	None
	"""
	ans_sgr = util.ans(answer)
	ans2int = lambda idch,pos: 1 if ans_sgr.isans(idch,pos) else -1
	
	for idch,start,pssm in parse_pssm4pos(fname):
		idch = idch.strip()
		# 2012/1/31 pos -> pos + start
		# !!! now nodyfying !!!
		yield idch,{pos:(ans2int(idch,start + pos),mkvec(pssm,window,pos))
					for pos in range(len(pssm)) if low <= ans_sgr.get_dist(pos + start,idch) <= up
					or ans_sgr.isans(idch,start + pos)}

def mkdtst_near4svr(fname,window,answer,low,up):
	# For test
	"""
	>>> [i for i in mkdtst_test('.test.pssm',10,answer='./.answer.test')]
	None
	"""
	ans_sgr = util.ans(answer)
	
	for idch,start,pssm in parse_pssm4pos(fname):
		idch = idch.strip()
		# 2012/1/31 pos -> pos + start
		# !!! now nodyfying !!!
		yield idch,{pos:(ans_sgr.get_dist(start + pos,idch),mkvec(pssm,window,pos))
					for pos in range(len(pssm)) if low <= ans_sgr.get_dist(pos + start,idch) <= up
					or ans_sgr.isans(idch,start + pos)}

def mkdtst_test4svr(fname,window,answer):
	# For test
	"""
	>>> [i for i in mkdtst_test('.test.pssm',10,answer='./.answer.test')]
	None
	"""
	ans_sgr = util.ans(answer)

	for idch,start,pssm in parse_pssm4pos(fname):
		idch = idch.strip()
		# 2012/1/31 pos -> pos + start
		# !!! now nodyfying !!!
		yield idch,{pos:(ans_sgr.get_dist(start + pos,idch),mkvec(pssm,window,pos)) for pos in range(len(pssm)) }

def mkdtst(fname,window):
	"""
	>>> import feature
	>>> [i for i in mkdtst_train('.test.pssm',10,answer='.answer.test')]
	None
	"""
	
	for pdbid,pssm in parse_pssm(fname):
		pdbid = pdbid.strip()
		# 2012/1/31 pos -> pos + start
		# !!! now nodyfying !!!
		yield pdbid,{pos:(-1,mkvec(pssm,window,pos)) for pos in range(len(pssm))}

class dataset(object):
	# Interface between dataset and libsvm
	
	def __init__(self,ppssm,npssm,answer,window):
		# pdata is data source of positive dataset
		self._ppssm = ppssm
		# ndata is data sourcde of negative dataset.
		self._npssm = npssm
		# answer is description of position of binding residue
		self._ans = answer

		# list of id in postive dataset
		self.pids = [idch for idch,start,pssm in parse_pssm4pos(ppssm)]
		# list of id in negative dataset 
		self.nids = [idch.strip(".fasta.ckp") for idch,pssm in parse_pssm(npssm)]

		self._window = window

	def mkTest(self,part_pids,part_nids):

		pdatasets = {idch:pvec for idch,pvec in mkdtst_test(self._ppssm,self._window,answer = self._ans) if idch in part_pids}
		ndatasets = {idch:pvec for idch,pvec in mkneg_test(self._npssm,self._window) if idch in part_nids}
		pdatasets.update(ndatasets)

		return pdatasets.items()
	
	def mkTrain(self,part_pids,part_nids,size = 5):
		"""
		>>> d = dataset(".test.pssm",".test.pssm",".answer.test",10)
		>>> d.mkTrain(['148lE'],['1af6_A_1_421'])
		None
		"""
		pdatasets = {idch:pvec for idch,pvec in mkdtst_train(self._ppssm,self._window,answer = self._ans) if idch in part_pids}
		ndatasets = {idch:pvec for idch,pvec in mkneg_train(self._npssm,self._window,size = size) if idch in part_nids}
		pdatasets.update(ndatasets)

		return pdatasets.items()

class dataset2(object):
	def __init__(self,pssm,answer,window,low = 5,up = 25):
		# pdata is data source of positive dataset
		self._ppssm = pssm
		# answer is description of position of binding residue
		self._ans = answer
		self._low = low
		self._up = up

		# list of id in postive dataset
		self.pids = [idch for idch,start,pssm in parse_pssm4pos(pssm)]
		self._window = window

	def mkTest(self,part_ids):
		pdatasets = {idch:pvec for idch,pvec in mkdtst_test(self._ppssm,self._window,answer = self._ans) if idch in part_ids}
		return pdatasets.items()
	
	def mkTrain(self,part_ids):
		# negative dataset are generate if low <= distance <= up from binding residue in dataset.
		"""
		>>> d = dataset(".test.pssm",".test.pssm",".answer.test",10)
		>>> d.mkTrain(['148lE'],['1af6_A_1_421'])
		None
		"""
		pdatasets = {idch:pvec for idch,pvec in mkdtst_near(self._ppssm,self._window,answer = self._ans,low = self._low,up = self._up) if idch in part_ids}
		return pdatasets.items()

class dataset3(dataset2):
	# For Test dataset whose answer position is Unkonown
	def __init__(self,pssm,window,answer,flg_predict = False):
		# pdata is data source of positive dataset
		self._ppssm = pssm
		# list of id in postive dataset
		

		if not flg_predict:
			self.pids = [idch for idch,start,pssm in parse_pssm4pos(pssm)]
			#self.pids = [idch.strip(".fasta.ckp") for idch,pssm in parse_pssm(pssm)]
			self._ans = answer
		else:
			# For not formated Fasta : > pdbid_start_end 
			self.pids = [idch for idch,pssm in parse_pssm(pssm)]
			self._ans = None
			
		self._window = window

	def mkTest(self,part_ids):
		"""
		>>> d = dataset3(".test.pssm",8,".answer.test")
		>>> d.mkTest(['148lE'])
		None
		"""
		pdatasets = {idch:pvec for idch,pvec in mkdtst_test(self._ppssm,self._window,self._ans) if idch in part_ids}
		return pdatasets.items()
	
	def ToPredict(self):
		"""
		>>> d = dataset3(".test.pssm",8,".answer.test")
		>>> d.mkTest(['148lE'])
		None
		"""
		pdatasets = {idch:pvec for idch,pvec in mkdtst(self._ppssm,self._window)}
		return pdatasets.items()
	
class dataset4(object):
	# For One class SVM
	def __init__(self,pssm,window):
		# pdata is data source of positive dataset
		self._ppssm = pssm

		self.pids = [idch for idch,pssm in parse_pssm(pssm)]
		#self.pids = [idch.strip(".fasta.ckp") for idch,pssm in parse_pssm(pssm)]
		self._window = window

	def mkTrain(self,part_ids):
		"""
		>>> d = dataset3(".test.pssm",8,".answer.test")
		>>> d.mkTest(['148lE'])
		None
		"""
		if part_ids is None:
			part_ids = self.pids
		return [pvec for idch,pvec in mkneg_train2(self._ppssm,self._window) if idch in part_ids]
	
	def mkTest(self,part_ids):
		"""
		>>> d = dataset3(".test.pssm",8,".answer.test")
		>>> d.mkTest(['148lE'])
		None
		"""
		if part_ids is None:
			part_ids = self.pids
		return {idch:pvec for idch,pvec in mkneg_test(self._ppssm,self._window) if idch in part_ids}.items()


class dataset5(object):
	def __init__(self,pssm,answer,window,low = 5,up = 25):
		# pdata is data source of positive dataset
		self._ppssm = pssm
		# answer is description of position of binding residue
		self._ans = answer
		self._low = low
		self._up = up

		# list of id in postive dataset
		self.pids = [idch for idch,start,pssm in parse_pssm4pos(pssm)]
		self._window = window

	def mkTest(self,part_ids):
		pdatasets = {idch:pvec for idch,pvec in mkdtst_test4svr(self._ppssm,self._window,answer = self._ans) if idch in part_ids}
		return pdatasets.items()
	
	def mkTrain(self,part_ids):
		# negative dataset are generate if low <= distance <= up from binding residue in dataset.
		"""
		>>> d = dataset(".test.pssm",".test.pssm",".answer.test",10)
		>>> d.mkTrain(['148lE'],['1af6_A_1_421'])
		None
		"""
		pdatasets = {idch:pvec for idch,pvec in mkdtst_near4svr(self._ppssm,self._window,answer = self._ans,low = self._low,up = self._up) if idch in part_ids}
		return pdatasets.items()

if __name__ == "__main__":
	import doctest
	doctest.testmod()
