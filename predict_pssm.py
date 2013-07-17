import argparse
import svm
import svmutil
from core import pssm2feature
from core import valid

PATH = "/home/ubuntu/galaxy-dist/tools/bilab/"

class predictor(object):
	# Contain Information of preditor.

	def __init__(self,model,reverse = False):
		self.model = svmutil.svm_load_model(model)
		if reverse:
			self.comp = lambda x,thr: x <= thr
		else:
			self.comp = lambda x,thr: x > thr

	def predict(self,dataset):
		ltest,dtest,itest = valid.test2svm_prob(dataset.ToPredict())
		plbl,pacc,pval = svmutil.svm_predict(ltest,dtest,self.model,"")
		return itest,pacc,pval
		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Cross Validation')
	# Input PSSM
	parser.add_argument('-i',action="store",dest="i")
	# Therehold
	#parser.add_argument('-thr',action="store",dest="thr",type = float)
	# reverse flag (Optional)
	parser.add_argument('-r',action="store_true",dest="r")
	
	#thr = parser.parse_args().thr
	data = parser.parse_args().i
	
	#model = svmutil.svm_load_model(model)
	d_cluster1 = pssm2feature.dataset3(data,4,None,flg_predict = True)
	d_acd = pssm2feature.dataset3(data,6,None,flg_predict = True)

	model_cluster1 = predictor(PATH + "models/pssm.d4.0.cluster1.last.4.0.-11.model")
	model_acd = predictor(PATH + "models/pssm.d4.0.acd.last.6.-1.-10.model",reverse = True)

	itest_cls,pacc_cls,pval_cls = model_cluster1.predict(d_cluster1)
	itest_acd,pacc_acd,pval_acd = model_acd.predict(d_acd)
	
	for i_cls,i_acd,val_cls,val_acd in zip(itest_cls,itest_acd,pval_cls,pval_acd):
		# dvalue is decision value of SVM
		val_cls = val_cls[0]
		val_acd = val_acd[0]
		result_cls = model_cluster1.comp(val_cls,0.588)
		# Best Parameter is 0.802. but model label is "-1" . so threshold is used as -0.802
		result_acd = model_acd.comp(val_acd,-0.802)
		if i_cls != i_acd:
			raise Exception
		# For the case using ":" in header.
		i_cls = i_cls.split(":")
		if len(i_cls) == 2:
			prot_id,pos = i_cls
		else:
			pos = i_cls[-1]
			prot_id = ":".join(i_cls[:2])
		result = [prot_id,pos,val_cls,val_acd,result_cls,result_acd]
		print "\t".join([str(i) for i in result])
