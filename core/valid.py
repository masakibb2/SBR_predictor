import time
import itertools
import os
import random
import sqlite3

import svm
import svmutil

from core import util
from core import seq2feature
from core import pssm2feature2

def crossed(groups,ngroups):
	while len(groups) > 1:
		yield groups.pop() + ngroups.pop(),groups + ngroups

def split_dtst(lid,k = 5,seed = None):
	"""
	>>> len([i for i in range(10,1000) if len(split_dtst(range(i))) != 5])
	0
	
	## PASS:split dataset into k part.

	>>> len(_recur_list(split_dtst(range(144))))
	144

	>>> print [len(i) for i in split_dtst(range(144))]
	None
	>>> print [len(i) for i in split_dtst(range(143))]
	None
	>>> print [len(i) for i in split_dtst(range(142))]
	None
	>>> print [len(i) for i in split_dtst(range(141))]
	None
	>>> print [len(i) for i in split_dtst(range(140))]
	None
	"""
	if seed is not None:
		# initalize seed.
		print "Using seed: %f" % seed
		random.seed(seed)

	random.shuffle(lid)
	n,m = divmod(len(lid),k)
	# amari no kazu dake +1 suru.
	_indx = [i*n + min(cnt + 1,m) for cnt,i in enumerate(range(1,k+1))]
	#_indx = [i*n for i in range(k + 1 - m)] + [i*n + 1 for i in range(k + 1 - m,k + 1)]
	#_indx[k] += m - 1 
	# !!! list[start:end] is not contained last element. 
	_indx[-1]+= 1
	start = 0
	divided = []
	for end in _indx:
		divided.append(lid[start:end])
		start = end
	return divided

def _recur_list(rlist):
	"""
	>>> _recur_list([[1,2,3],4,[5,6,[7,8,9,[10]]]])
	>>> [1,2,3,4,5,6,7,8,9,10]
	"""
	out = []
	for i in rlist:
		if isinstance(i,list):
			out += _recur_list(i)
		else:
			out.append(i)
	return out

def fold(ids,k = 5,seed = None):
	# !!! must be same length groups and ngroups !!!
	"""
	>>> [(i,j) for i,j in fold([[i] for i in range(5)])]
	None
	"""
	# split dataset into k part
	groups = split_dtst(ids,k,seed)
	
	for i in range(k):
		test = groups[i]
		train = groups[:i] + groups[i+1:]
		if isinstance(test,list):
			test = _recur_list(test)
		if isinstance(train,list):
			train = _recur_list(train)
		yield test,train

def lbl_iter(test):
    for key,i in test:
        for j in i.values():
            yield j[0]

def data_iter(test):
    for key,i in test:
        for j in i.values():
            yield j[1]


def id_iter(test):
    for key,i in test:
        for j in i.keys():
            yield key + ':' + str(j)


def test2svm_prob(test):
    ltest = [l for l in lbl_iter(test)]
    dtest = [d for d in data_iter(test)]
    itest = [i for i in id_iter(test)]
	
    return ltest,dtest,itest


def train2svm_prob(train):
	ltrain = [l for l in lbl_iter(train)]
	dtrain = [d for d in data_iter(train)]
	itrain = [i for i in id_iter(train)]

	return svm.svm_problem(ltrain,dtrain),itrain
	


class valid(object):
	# For 2 class SVM.
	def __init__(self,name,dirname):
		# name is used as log file name
		# dirname is location saving log file
		self._name = name
		self._dir = dirname

	def _save_log(self,itest,plbl,pval,cnt):
		"""
		>>> import feature
		>>> import seq2feature
		>>> v = valid(".test","./.test")
		>>> d = seq2feature.dataset(".test.fasta",".test.neg.fasta","..//dataset/answer_monod4.0.cluster1.txt",lambda seq:feature.seq2frq(seq,2),10)
		>>> v.valid(d,opt = "-c %f -g %f" % (1.0,1.0),opp = "")
		"""
		
		# create directory of saving logfile
		with open("%s/%s.log.%s" % (self._dir,self._name,cnt),'a') as fp:
			for idch,iter_result in itertools.groupby(zip(itest,plbl,pval),key = lambda x:x[0].split(':')[0]):
				for pos,lbl,val in iter_result:
					fp.write("%s\t%s\t%s\n" % (pos,lbl,val[0]))

	def valid(self,datasets,opt,opp,method = fold,seed = None):
		# Should groups and ngroups be idch ?
		groups = [(test,train) for test,train in method(datasets.pids,seed = seed)]
		ngroups = [(test,train) for test,train in method(datasets.nids,seed = seed)]
		
		for cnt,(pdtsts,ndtsts) in enumerate(zip(groups,ngroups)):
			# cnt is number of cluster.
			ltest,dtest,itest = test2svm_prob(datasets.mkTest(pdtsts[0],ndtsts[0]))
			ptrn,itrain = train2svm_prob(datasets.mkTrain(pdtsts[1],ndtsts[1]))
			
			print "start %s validation" % (cnt)
			#opt = svm.svm_parameter(opt)
			model = svmutil.svm_train(ptrn,opt)
			plbl,pacc,pval = svmutil.svm_predict(ltest,dtest,model,opp)


			# create saving direcotry
			#self._mkdir(cnt)
			# create log files
			self._save_log(itest,plbl,pval,cnt)
			model_name = "%s/model/%s.model.%s" % (self._dir,self._name,cnt)
			#svmutil.svm_save_model(model_name, model)


class mkDB(object):
	def __init__(self,iter_record,mktable,table_name):
		# iter_record is iterator of yielding records.
		# mktable is sql statement of create table
		# table_name is created table name by mktable().
		
		self._iter_rec = iter_record
		self._mkTable = mktable
		self._tbl = table_name
	
	def mkTable(self,con):
		con.execute(self._mkTable)
		
	def updtDB(self,con):
		for record in self._iter_rec:
			sql = lambda record,table : "insert into %s values ( " % (table) + ", ".join(["?" for i in range(len(record))]) + ");"
			con.execute(sql(record,self._tbl),record)


class valid2(valid):
	"""
	
	save cross validation log using sqlite3.
	
	"""
	def __init__(self,name,dirname,fname):
		valid.__init__(self,name,dirname)
		self._starts = {idch:start for idch,start,seq in seq2feature.fasta2seq(fname)}
		
	def _iter_result(self,itest,plbl,pval,cnt):
		for idch,iter_result in itertools.groupby(zip(itest,plbl,pval),key = lambda x:x[0].split(':')[0]):
			if self._starts.has_key(idch):
				start = self._starts[idch]
			else:
				start = 0
			for pos,lbl,val in iter_result:
				pos = start + int(pos.split(":")[1])
				yield idch,pos,val[0],cnt
	
	def _save_log(self,itest,plbl,pval,cnt):
		saving_db = "%s/log/%s.log.db" % (self._dir,self._name)
		with sqlite3.connect(saving_db) as con:
			#is_ans bool, answer is written in answer database.
			# We can summrize below sql statement.
			# ( where valid.idch = answr.idch and valid.pos = answer.bp + answer.start )
			
			mktbl = """
			create table valid (
			idch text,
			pos interger,
			dec_val real,
			cnt interger

			);"""
			
			man_db = mkDB(self._iter_result(itest,plbl,pval,cnt),mktbl,"valid")
			if cnt == 0:
				man_db.mkTable(con)
			
			man_db.updtDB(con)
			con.commit()
		
class valid3(valid2):
	# change self.valid() method.
	def valid(self,datasets,opt,opp,method = fold,part_ids = None,seed = None,test_data = None):
		if seed is None:
			# If seed is not set. UNIX time is used as seed.
			seed = time.time()
		saving_seed = "%s/log/%s.log.seed" % (self._dir,self._name)
		with open(saving_seed,"w") as fp:
			# Save used seed value.
			fp.write("seed:%f\n" % seed)
		
		if part_ids is None:
			part_ids = datasets.pids
		groups = [(test,train) for test,train in method(part_ids,seed = seed)]
		
		for cnt,pdtsts in enumerate(groups):
			# cnt is number of cluster.
			if test_data is None:
				test = False
				ltest,dtest,itest = test2svm_prob(datasets.mkTest(pdtsts[0]))
			else:
				test = True
				ltest,dtest,itest = test2svm_prob(test_data.mkTest(test_data.pids))

			print "start %s validation" % (cnt)
			ptrn,itrain = train2svm_prob(datasets.mkTrain(pdtsts[1]))
			#opt = svm.svm_parameter(opt)
			model = svmutil.svm_train(ptrn,opt)
			
			plbl,pacc,pval = svmutil.svm_predict(ltest,dtest,model,opp)

			# create saving direcotry
			#self._mkdir(cnt)
			# create log files
			self._save_log(itest,plbl,pval,cnt,test)
			model_name = "%s/model/%s.model.%s" % (self._dir,self._name,cnt)
			#svmutil.svm_save_model(model_name, model)
	
	def create_model(self,datasets,opt,opp,part_ids = None):
		# Should groups and ngroups be idch ?
		if part_ids is None:
			part_ids = datasets.pids
		ptrn,itrain = train2svm_prob(datasets.mkTrain(part_ids))
		print "create model ..."
		#opt = svm.svm_parameter(opt)
		model = svmutil.svm_train(ptrn,opt)
		# create saving direcotry
		#self._mkdir(cnt)
		# create log files
		#self._save_log(itest,plbl,pval,cnt)
		model_name = "%s/model/%s.model" % (self._dir,self._name)
		svmutil.svm_save_model(model_name, model)

class valid4(valid3):
	# save distance from binding site.
	def __init__(self,name,dirname,fname,fans):
		valid3.__init__(self,name,dirname,fname)
		if fans is not None:
			self._ans = util.ans(fans)
		
	def _iter_result(self,itest,plbl,pval,cnt,test):
		for idch,iter_result in itertools.groupby(zip(itest,plbl,pval),key = lambda x:x[0].split(':')[0]):
			for pos,lbl,val in iter_result:
				if not test:
					start = self._starts[idch]
					pos = start + int(pos.split(":")[1])
					dist = self._ans.get_dist(pos,idch)
				else:
					dist = 0
				yield idch,pos,dist,val[0],cnt
	
	def _save_log(self,itest,plbl,pval,cnt,test = False):
		saving_db = "%s/log/%s.log.db" % (self._dir,self._name)
		with sqlite3.connect(saving_db) as con:
			#is_ans bool, answer is written in answer database.
			# We can summrize below sql statement.
			# ( where valid.idch = answr.idch and valid.pos = answer.bp + answer.start )
			
			mktbl = """
			create table valid (
			idch text,
			pos interger,
			dist interger,
			dec_val real,
			cnt interger

			);"""
			
			man_db = mkDB(self._iter_result(itest,plbl,pval,cnt,test = test),mktbl,"valid")
			if cnt == 0:
				man_db.mkTable(con)
			
			man_db.updtDB(con)
			con.commit()



class valid5(valid4):
	def _iter_result(self,tst_d_lbl,dec_vals,cnt):
		for idch_pos,dec_val in zip(tst_d_lbl,dec_vals):
			idch,pos = idch_pos
			dec_val = dec_val[0]
			dist = self._ans.get_dist(pos,idch)
			yield idch,pos,dist,dec_val,cnt

	def _save_log(self,tst_d_lbl,dec_vals,cnt):
		saving_db = "%s/log/%s.log.db" % (self._dir,self._name)
		with sqlite3.connect(saving_db) as con:
			#is_ans bool, answer is written in answer database.
			# We can summrize below sql statement.
			# ( where valid.idch = answr.idch and valid.pos = answer.bp + answer.start )
			
			mktbl = """
			create table valid (
			idch text,
			pos interger,
			dist interger,
			dec_val real,
			cnt interger

			);"""
			
			man_db = mkDB(self._iter_result(tst_d_lbl,dec_vals,cnt),mktbl,"valid")
			if cnt == 0:
				man_db.mkTable(con)
			
			man_db.updtDB(con)
			con.commit()
	
	def valid(self,datasets,C,gamma,class_weight=None,method = fold,filter = None):
		# Should groups and ngroups be idch ?
		groups = [(test,train) for test,train in method(datasets.pids)]
		
		for cnt,pdtsts in enumerate(groups):
			# cnt is number of cluster.
			tst_d_lbl,tst_lbl,tst_dtst = datasets.mkTest(part_ids=pdtsts[0])
			trn_d_lbl,trn_lbl,trn_dtst = datasets.mkTrain(part_ids=pdtsts[1])
			
			print "start %s validation" % (cnt)
			clf = SVC(C= C ,gamma = gamma,class_weight = class_weight)

			tst_lbl,tst_dtst = numpy.array(tst_lbl),numpy.array(tst_dtst)
			trn_lbl,trn_dtst = numpy.array(trn_lbl),numpy.array(trn_dtst)
			
			if filter is not None:
				# For Univarable Feature Selection.
				filter.fit(trn_dtst,trn_lbl)
				tst_dtst = filter.transform(tst_dtst)
				trn_dtst = filter.transform(trn_dtst)
			clf.fit(trn_dtst,trn_lbl)
			dec_vals = clf.decision_function(tst_dtst)
			self._save_log(tst_d_lbl,dec_vals,cnt)
	
class valid5_2(valid5):
	# save distance from binding site.
	def __init__(self,name,dirname,fname,fans):
		self._name = name
		self._dir = dirname
		self._ans = util.ans(fans)
		pssmpp = pssm2feature2.pssm(fname,length = 1)
		self.starts = {idch:start for idch,start,v_pssmpp in pssmpp.parse_pssm4pos()}
			
if __name__ == "__main__":
	import doctest
	doctest.testmod()
