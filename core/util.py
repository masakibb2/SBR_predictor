import sys

try:
	import yaml
except:
	sys.stderr.write("Cannot import yaml.\n")
	pass

class ans(object):
	# For make positive dataset
	def __init__(self,fans):
            with open(fans) as fp:
				self._ans = {i[0]:[int(pos) for pos in i[1:]]
                             for i in (line.strip().split() for line in iter(fp.readline,""))}

	def isans(self,pdbid,nsq):
            if not self._ans.has_key(pdbid):
                # not in ID => Negative set
                return False
            if nsq in self._ans[pdbid]:
                return True
            else:
                return False
	def get_pos(self,pdbid):
		if not self._ans.has_key(pdbid):
			print "## warnings pdbid = %s doesn't exists" % (pdbid)
			return None
		else:
			return self._ans[pdbid]
	
	def get_dist(self,pos,idch):
		# most near residue number for each binding site.
		if self.get_pos(idch) is None:
			return None
		near = min(self.get_pos(idch),key = lambda x:abs(x - pos))
		return abs(pos - near)

class surf(object):
	# For make positive dataset
	def __init__(self,fsurf):
		with open(fsurf) as fp:
			self.surf = yaml.load(fp.read())

	def get_surf(self,pos,idch):
		# most near residue number for each binding site.
		if not self.surf.has_key(idch):
			return -2
		if self.surf[idch].has_key(pos):
			return self.surf[idch][pos]
		else:
			return -1
	def get_pos(self,idch):
		if not self.surf.has_key(idch):
			return None
		return self.surf[idch]

class blust(object):
	"""
	>>> idch2clust = "../dataset/idch2clust_num.txt"
	>>> clust2idch = "../dataset/clust2idch.txt"
	>>> mapper = blust(idch2clust,clust2idch,3)

	# if input idch not have clust. return None
	>>> mapper.get_clust("148lE") is None
	True
	>>> mapper.get_clust("3gpbA")
	0

	# if idch is same cluster, retrun same number.
	>>> mapper.get_clust("2av6B")
	0
	>>> mapper.get_clust("3ic3B")
	372

	# if input idch not have clust. return None
	>>> mapper.conv_clust("148lE",2) is None
	True

	# convert idch to idch of input class.
	>>> mapper.conv_clust("2av6B",0)
	'3gpbA'

	# return same idch if idch is same as idch of class.
	>>> mapper.conv_clust("2av6B",2)
	'2av6B'
	
	#  return 'Nan' if idch cannot map to input class_ID.
	>>> mapper.conv_clust("3ic3B",0)
	'Nan'

	# return True if input idch in input class_ID.
	>>> mapper.is_clust("3ic3B",0)
	False
	>>> mapper.is_clust("148lE",0)
	False
	>>> mapper.is_clust("3ic3B",1)
	True

		"""
	clust_ID = {"acd":0,"cls":1,"mono":2}
	
	def __init__(self,idch2clust,clust2idch,class_size):
		# Input number of cluster
		self.class_size = class_size
		# mapper of idch -> clust_ID
		self.idch2clust = {}
		with open(idch2clust) as fp:
			for line in iter(fp.readline,""):
				record = line.strip().split()
				idch = record[0]
				clust = int(record[1])
				self.idch2clust.update({idch:int(clust)})

		# mapper of clut_ID -> idch
		self.clust2idch = {}
		with open(clust2idch) as fp:
			for line in iter(fp.readline,""):
				record = line.strip().split()
				clust = record[0]
				idchs = record[1:]
				self.clust2idch.update({int(clust):idchs})

	def get_clust(self,idch):

		if self.idch2clust.has_key(idch):
			return self.idch2clust[idch]
		else:
			return None

	def conv_clust(self,idch,class_ID):
		# Error if class_ID > class size -1.
		if class_ID >= self.class_size:
			raise Exception,"Unvalid class_ID : %s, this dataset class size is %s" % (class_ID,self.class_size)
		clust_ID = self.get_clust(idch)

		if clust_ID is None:
			sys.stderr.write("%s is not in class_ID = %s\n" % (idch,class_ID))
			return None
		else:
			# if input class_ID have not same cluster. return "Nan"
			if idch != self.clust2idch[clust_ID][class_ID]:
				sys.stderr.write("convert %s to %s\n" % (idch,self.clust2idch[clust_ID][class_ID]))
			return self.clust2idch[clust_ID][class_ID]

	
	def is_clust(self,idch,class_ID):
		result = self.conv_clust(idch,class_ID)
		if result is None or result == "Nan":
			return False
		else:
			return result == idch

if __name__ == "__main__":
	import doctest
	doctest.testmod()
