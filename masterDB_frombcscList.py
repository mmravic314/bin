# from blastClust list of single chains clustered

import sys, os, subprocess as sp, cPickle as pic
from prody import *
from collections import defaultdict


#pre_dir 	= sys.argv[2]
#finlaDir 	= sys.argv[3]
# python ~/bin/masterDB_frombcscList.py entries_xray_res.txt pdb_seqres-byCh.txt 152901_bc_30-scPDB.txt 


clusters = []

class PDBentry:

	def __init__(self, name, res = 0.0):
		self.res	= res
		self.name	= name

	chains			= {}	  ## tuple of chains per length
	def __repr__(self):
		return self.name

	def getinfo(self, ch):
		return ( self.res, self.chains[ch] )


class PDBchain( PDBentry ):

	def __init__(self, length):
		self.length = length



### Pre-compute section... for all PDB chains, determine and sort the size of each chain along with its XRAY data resolution
# skip non-xray proteins

topLevelDict	 = {}

# Look through files in following order... to gain info about each PDB chain
# 1) diffraction type and resolution. dump all not X-ray or Neutron
with open ( sys.argv[1] ) as file:
	for i in file:

		name, exp = i[:4], i.split()[-5:]

		# skip NMR or non-atomic resoultion methods
		if exp[0] == "NOT" or exp[2]=='NOT':
			continue

		# Collect joint neutron-xray
		if exp[0][0].isdigit() and exp[-2] == 'NEUTRON' or exp[0][0].isdigit() and exp[1] == 'NEUTRON':
			topLevelDict[ name ]	 = PDBentry(name) 
			topLevelDict[ name ].res = float( exp[0] )
#			print "neutron2", name, topLevelDict[ name ].res
			continue

		# collect neutron exps
		if exp[-2] == 'NEUTRON':  
			topLevelDict[ name ] 	 = PDBentry(name) 
			topLevelDict[ name ].res = float( exp[2] )
#			print "neutron", name, topLevelDict[ name ].res
			continue


		# collect xrays
		if exp[-2] == 'X-RAY' and exp[-1] == 'DIFFRACTION' :
			
			topLevelDict[ name ] 	 = PDBentry(name) 
			topLevelDict[ name ].res = float( exp[2] )
#			print 'xray', name, topLevelDict[ name ].res
			continue


# 2) fasta file with chain length... store chain length
with open( sys.argv[2] ) as file:
	for i in file:

		if i[0] != '>' or i.split()[1][4:].strip()== 'na': continue
		name, ch, length = i[1:5].upper(), i[6:8].strip(), int( i.split()[2][7:].strip() )
		if length < 10:
				continue
		try:

			topLevelDict[ name ].chains[ch] = length
#			print name, ch, topLevelDict[ name ].getinfo(ch)
		except KeyError:
			pass

#sys.exit()
count = 0
# The finally look through the cluster list and find best resolution and largest proteni in each cluster
with open( sys.argv[3] ) as file:
	for i in file:
		# in each cluster, rank chains first by resolution 
		
	#	best = [   topLevelDict[ x[:4] ].chains[ x[5:] ]
		valid = []

		for x in i.split():
			try: 
				valid.append( (  x, topLevelDict[ x[:4] ].getinfo( x[5:] )[0], topLevelDict[ x[:4] ].getinfo( x[5:] )[1] ) )
			except KeyError:
				pass

		best = ''
		for val in sorted( valid, key=lambda x: (x[1], x[2]) ):
			if val[1] > 2.8:
				break
			else:
				best = val[0]
				break

		if len( best ) > 1:
			print best, count 
			count += 1

