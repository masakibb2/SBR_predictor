import random
import util

class pssm(object):
	def __init__(self,fname,length = 20):
		self._length = length
		self._fname = fname
	
	class _pssm(object):
		def __init__(self,pdbid,length):
			self.pdbid = pdbid
			self.pssm = []
			self.length = length
		
		def append(self,vec):
			if len(vec) != self.length:
				print "pssm.append(): Argument length must be %s. len(vec) == %s" % (len(vec),len(vec))
				raise Exception,"pssm.append(): Argument length must be %s. len(vec) == %s" % (len(vec),len(vec))
			self.pssm.append(vec)
			

	def parse_pssm(self):
		# parse pssm
		"""
		>>> {k:v for k,v in parse_pssm(".test.pssm")}
		True
		"""
		prot_pssm = None
		with open(self._fname) as fp:
			for rec in (line.strip() for line in iter(fp.readline,"")):
				if len(rec) == 0:
					continue
				if rec[0] == '>':
					if prot_pssm is not None:
						yield prot_pssm.pdbid,prot_pssm.pssm
					prot_pssm = pssm._pssm(rec[1:].strip(),self._length)
				else:
					#if rec.count(":") > 0:
					#	prot_pssm = pssm._pssm(rec[1:].strip(),length = self._length)
					#else:
					prot_pssm.append([float(i) for i in rec.split()])
				
			else:
				yield prot_pssm.pdbid,prot_pssm.pssm

	def parse_pssm4pos(self):
		"""
		>>> [i for i in parse_pssm4pos(".test.pssm")]
		True
		"""
		# yield idch,start position, sequence
		for prot_id,v_pssm in self.parse_pssm():
			pdbid,ch,start,end = prot_id.split('_')
			idch = pdbid + ch
			#idch,start,end = prot_id.split('_')
			#idch,start = prot_id.split('_')
			yield idch.replace('_',''),int(start),v_pssm

	def mkvec(self,v_pssm,window,pos):
		"""
		>>> _pssm = [i for i in parse_pssm4pos(".test.pssm")][0][2]
		>>> mkvec(_pssm,10,3)
		None
		"""
		def _paste(vec):
			out = []
			for i in vec:
				out+=i
			return out

		nu = [0]*self._length
		n = v_pssm[max(0,int(pos - window/2.0)):pos]
		c = v_pssm[pos + 1:min(int(pos + window/2.0) + 1,len(v_pssm))]
		if window%2 == 1:
			nu_n = [nu for i in range(int(window/2 - len(n) + 1))]
		else:
			nu_n = [nu for i in range(int(window/2 - len(n)))]
		nu_c = [nu for i in range(int(window/2.0 - len(c)))]
		vec = nu_n + n + [v_pssm[pos]] + c + nu_c
		return _paste(vec)


	def mkdtst_train(self,window,answer):
		"""
		>>> import feature
		>>> [i for i in mkdtst_train('.test.pssm',10,answer='.answer.test')]
		None
		"""
		# only binding site For make positive dataet
		ans_sgr = util.ans(answer)
		# ans2int is function that return +1 if True else return -1
		ans2int = lambda idch,pos: 1 if ans_sgr.isans(idch,pos) else -1
	
		for idch,start,v_pssm in self.parse_pssm4pos():
			idch = idch.strip()
			idchs = [(idch,pos + start) for pos in range(len(v_pssm))]
			dataset = [self.mkvec(v_pssm,window,pos) for pos in range(len(v_pssm)) if ans_sgr.isans(idch,start + pos)]
			label = [ans2int(idch,start + pos) for pos in range(len(v_pssm)) if ans_sgr.isans(idch,start + pos)]

			yield idch,idchs,label,dataset

	def mkdtst_neg(self,window,answer,size = 1):
		"""
		>>> import feature
		>>> [i for i in mkdtst_neg('.test.pssm',10,'.answer.test',1)]
		None
		"""

		ans_sgr = util.ans(answer)
		# To Dataset randomly choiced positions
		for idch,start,pssm in self.parse_pssm4pos():
			_indx = [i for i in range(1,len(pssm)) if not ans_sgr.isans(idch,i + start)]
			random.shuffle(_indx)
			selected_pos = _indx[:size*len(ans_sgr.get_pos(idch))]
			idch = idch.strip()
			idchs = [(idch,pos + start) for pos in selected_pos]
			dataset = [self.mkvec(pssm,window,pos) for pos in selected_pos]
			label = [-1 for pos in selected_pos]
			yield idch,idchs,label,dataset

	
	def mkdtst_test(self,window,answer):
		# For test
		"""
		>>> [i for i in mkdtst_test('.test.pssm',10,answer='./.answer.test')]
		None
		"""
	
		ans_sgr = util.ans(answer)
		ans2int = lambda idch,pos: 1 if ans_sgr.isans(idch,pos) else -1
	
		for idch,start,v_pssm in self.parse_pssm4pos():
			idch = idch.strip()
			idchs = [(idch,pos + start) for pos in range(len(v_pssm))]
			dataset = [self.mkvec(v_pssm,window,pos) for pos in range(len(v_pssm))]
			label = [ans2int(idch,start + pos) for pos in range(len(v_pssm))]
			yield idch,idchs,label,dataset

	def mkdtst_near(self,window,answer,low,up):
		# For test
		"""
		>>> [i for i in mkdtst_test('.test.pssm',10,answer='./.answer.test')]
		None
		"""
		ans_sgr = util.ans(answer)
		ans2int = lambda idch,pos: 1 if ans_sgr.isans(idch,pos) else -1
	
		for idch,start,v_pssm in self.parse_pssm4pos():
			idch = idch.strip()
			idchs = [(idch,pos + start) for pos in range(len(v_pssm))
					 if low <= ans_sgr.get_dist(pos + start,idch) <= up or ans_sgr.isans(idch,start + pos)]
			dataset = [self.mkvec(v_pssm,window,pos) for pos in range(len(v_pssm))
			         if low <= ans_sgr.get_dist(pos + start,idch) <= up or ans_sgr.isans(idch,start + pos)]
			label = [ans2int(idch,start + pos) for pos in range(len(v_pssm))
			         if low <= ans_sgr.get_dist(pos + start,idch) <= up or ans_sgr.isans(idch,start + pos)]
			yield idch,idchs,label,dataset


class dataset(object):
	# Interface between dataset and libsvm
	
	def __init__(self,ppssm,npssm,answer,window,length = 20):
		# pdata is data source of positive dataset
		self._ppssm = pssm(ppssm,length)
		# ndata is data sourcde of negative dataset.
		self._npssm = pssm(npssm,length)
		# answer is description of position of binding residue
		self._ans = answer

		# list of id in postive dataset
		self.pids = [idch for idch,start,v_pssm in self._ppssm.parse_pssm4pos()]
		# list of id in negative dataset 
		self.nids = [idch.strip(".fasta.ckp") for idch,v_pssm in self._npssm.parse_pssm()]

		self._window = window

	def mkTest(self,part_pids,part_nids):

		pdatasets = {idch:pvec for idch,pvec in self._ppssm.mkdtst_test(self._window,answer = self._ans) if idch in part_pids}
		ndatasets = {idch:pvec for idch,pvec in self._npssm.mkneg_test(self._window) if idch in part_nids}
		pdatasets.update(ndatasets)

		return pdatasets.items()
	
	def mkTrain(self,part_pids,part_nids,size = 5):
		"""
		>>> d = dataset(".test.pssm",".test.pssm",".answer.test",10)
		>>> d.mkTrain(['148lE'],['1af6_A_1_421'])
		None
		"""
		pdatasets = {idch:pvec for idch,pvec in self._ppssm.mkdtst_train(self._window,answer = self._ans) if idch in part_pids}
		ndatasets = {idch:pvec for idch,pvec in self._npssm.mkneg_train(self._window,size = size) if idch in part_nids}
		pdatasets.update(ndatasets)

		return pdatasets.items()

class dataset2(object):
	def __init__(self,ppssm,answer,window,low = 5,up = 25,length = 20):
		# pdata is data source of positive dataset
		self._ppssm = pssm(ppssm,length)
		# answer is description of position of binding residue
		self._ans = answer
		self._low = low
		self._up = up

		# list of id in postive dataset
		self.pids = [idch for idch,start,ppssm in self._ppssm.parse_pssm4pos()]
		self._window = window

	def mkTest(self,part_ids):
		data_label = []
		labels = []
		dataset = []
		for idch,idchs,label,data in self._ppssm.mkdtst_test(self._window,answer = self._ans):
			if idch in part_ids:
				data_label += idchs
				labels += label
				dataset += data

		return data_label,labels,dataset

	def mkTrain(self,part_ids):
		# negative dataset are generate if low <= distance <= up from binding residue in dataset.
		"""
		>>> d = dataset(".test.pssm",".test.pssm",".answer.test",10)
		>>> d.mkTrain(['148lE'],['1af6_A_1_421'])
		None
		"""
		data_label = []
		labels = []
		dataset = []
		
		for idch,idchs,label,data in self._ppssm.mkdtst_near(self._window,answer = self._ans,low = self._low,up = self._up):
			if idch in part_ids:
				data_label += idchs
				labels += label
				dataset += data
		
		return data_label,labels,dataset



class dataset3(dataset2):
	def __init__(self,ppssm,answer,window,length = 20):
		# Under sammpling dataset for Random Forest
		# pdata is data source of positive dataset
		self._ppssm = pssm(ppssm,length)
		# answer is description of position of binding residue
		self._ans = answer
		# list of id in postive dataset
		self.pids = [idch for idch,start,ppssm in self._ppssm.parse_pssm4pos()]
		self._window = window
	
	def mkTrain(self,part_ids):
		data_label = []
		labels = []
		dataset = []
		
		for idch,idchs,label,data in self._ppssm.mkdtst_train(self._window,answer = self._ans):
			if idch in part_ids:
				data_label += idchs
				labels += label
				dataset += data
		
		for idch,idchs,label,data in self._ppssm.mkdtst_neg(self._window,answer = self._ans,size = 1):
			if idch in part_ids:
				data_label += idchs
				labels += label
				dataset += data

		return data_label,labels,dataset

if __name__ == "__main__":
	import doctest
	doctest.testmod()

