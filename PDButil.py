## PDB file utilities ##
##################### Written By: Marco Mravic  ###  Degrado Lab UCSF Biophysics Aug 2015

import numpy as np
from prody import * 
from collections import Counter

###################### Globally defined helpers  ####################

divalent_metals = ['BA','BE', 'CD', 'CA', 'CR', 'CO', 'CU', 
'EU', 'GD', 'GE', 'FE', 'LA', 'PD', 'MG', 'MN', 'HG', 'NI', 
'OS', 'PT', 'RU', 'SR', 'SN','U', 'V', 'Y', 'ZN' ]

natAA = {}
natAA["ALA"] = 'A'; natAA["CYS"] = 'C'; natAA["ASP"] = 'D'; natAA["GLU"] = 'E'; natAA["PHE"] = 'F';
natAA["GLY"] = 'G'; natAA["HIS"] = 'H'; natAA["ILE"] = 'I'; natAA["LYS"] = 'K';
natAA["LEU"] = 'L'; natAA["MET"] = 'M'; natAA["ASN"] = 'N'; natAA["PRO"] = 'P'; natAA["GLN"] = 'N';
natAA["ARG"] = 'R'; natAA["SER"] = 'S'; natAA["THR"] = 'T'; natAA["VAL"] = 'V'; natAA["TRP"] = 'W'; natAA["TYR"] = 'Y';

# Arbitrary amino acid distribution from a database of PDBs
aaProp = {}
aaProp["ALA"] = 7.73; aaProp["CYS"] = 1.84; aaProp["ASP"] = 5.82; aaProp["GLU"] = 6.61; aaProp["PHE"] = 4.05;
aaProp["GLY"] = 7.11; aaProp["HIS"] = 2.35; aaProp["HSD"] = 2.35; aaProp["ILE"] = 5.66; aaProp["LYS"] = 6.27;
aaProp["LEU"] = 8.83; aaProp["MET"] = 2.08; aaProp["ASN"] = 4.50; aaProp["PRO"] = 4.52; aaProp["GLN"] = 3.94;
aaProp["ARG"] = 5.03; aaProp["SER"] = 6.13; aaProp["THR"] = 5.53; aaProp["VAL"] = 6.91; aaProp["TRP"] = 1.51; aaProp["TYR"] = 3.54;

## Unnatural/PTM Amino acid converter
unNatAA = { 'ABA':'ALA', 
'CSO':'CYS' , 'CSD':'CYS', 'CME':'CYS', 'OCS':'CYS', 
"HSD":'HIS',
'KCX':'LYS', 'LLP':'LYS', 'MLY':'LYS', 'M3L':'LYS', 
'MSE':'MET', 
'PCA':'PRO', 'HYP':'PRO',
'SEP':'SER', 'TPO':'THR', 'PTR':'TYR'
  }

# Same Three to One res ID converter but includes/converts unnatural/modified amino acids
UnNatAA = {}
UnNatAA["ALA"] = 'A'; UnNatAA["CYS"] = 'C'; UnNatAA["ASP"] = 'D'; UnNatAA["GLU"] = 'E'; UnNatAA["PHE"] = 'F';
UnNatAA["GLY"] = 'G'; UnNatAA["HIS"] = 'H'; UnNatAA["ILE"] = 'I'; UnNatAA["LYS"] = 'K';
UnNatAA["LEU"] = 'L'; UnNatAA["MET"] = 'M'; UnNatAA["ASN"] = 'N'; UnNatAA["PRO"] = 'P'; UnNatAA["GLN"] = 'N';
UnNatAA["ARG"] = 'R'; UnNatAA["SER"] = 'S'; UnNatAA["THR"] = 'T'; UnNatAA["VAL"] = 'V'; UnNatAA["TRP"] = 'W'; UnNatAA["TYR"] = 'Y';
UnNatAA['ABA'] = 'A'; UnNatAA['CSO'] = 'C'; UnNatAA['CSD'] = 'C'; UnNatAA['CME'] = 'C';
UnNatAA['OCS'] = 'C'; UnNatAA["HSD"] = 'H'; UnNatAA['KCX'] = 'K'; UnNatAA['LLP'] = 'K';
UnNatAA['MLY'] = 'K'; UnNatAA['M3L'] = 'K'; UnNatAA['MSE'] = 'M'; UnNatAA['PCA'] = 'P'; UnNatAA['HYP'] = 'P';
UnNatAA['SEP'] = 'S'; UnNatAA['TPO'] = 'T'; UnNatAA['PTR'] = 'Y'


#################### End of Global Definitions ####################



############### Class def's ##############


# mapped PDB chain containing which pdb residue corresponds to each uniprot entry 
#Has optiono to initialize or update a chain with a DBREF line (pdbfile) or list of segments from EMBL SIFTS mapping .tsv file 
class PdbChain():
	
	# does mapping at initialization given a DBREF line from a PDB file
	def __init__(self, init_flg, init_info):
		self.name 		= ''
		self.uniprot 	= ''
		self.pdbRes 	= np.array( [] )
		self.uniRes  	= np.array( [] )
		self.exceptFlg 	= False				# This flag is True if the   

		if init_flg == 'pdb':
			self.pdbLine_init(init_info)	# initialize with DBREF lines
		elif init_flg == 'SIFTS':
			self.SIFTS_init(init_info)	# initlize with list of segments in SIFTS .tsv file format
		else:						#Stay empty
			pass


	def __repr__(self):
		return '%s, %s' % (self.name, self.uniprot) 

	# For each chain in pdb entry, should have uniprot mapping at chain level (except large structures or synthesized constructs)
	def pdbLine_init(self, pdbLine):
		self.name = pdbLine[7:11] + pdbLine[12].lower() 
		self.uniprot = pdbLine[33:39]
		self.pdbRes = np.arange( int( pdbLine[14:18].strip() ) , int ( pdbLine[19:24].strip() ) + 1 )
		self.uniRes  = np.arange( int( pdbLine[55:60].strip() ) , int ( pdbLine[62:68].strip() ) + 1 )
		return

	# reads in format from EMBL-EBI mapping https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html, pdb to uniprot
	def SIFTS_init(self, lineList):
		if len(lineList) == 0: print 'WHY', self.name, self.uniprot
		for i in lineList:
			ln 			= i.rstrip().split(',')
			chain 		= ''.join( [ln[0].upper(), ln[1].lower() ] )
			uniprot 	= ln[2]
			try:
				p_res 	= np.arange(  int( ln[5] ) , int( ln[6] ) +1 )
			except ValueError:
				p_res 	= np.array( [] )
				self.exceptFlg = True
			try:	
				u_res 		= np.arange(  int( ln[7] ) , int( ln[8] ) +1 )
			except ValueError:
				u_res 	= np.array( [] )
				self.exceptFlg = True

			self.pdbRes = np.append( self.pdbRes, p_res )
			self.uniRes = np.append( self.uniRes, u_res )

		self.name		= chain
		self.uniprot 	= uniprot
		return


	def pdb2uni(self, resNum):
		val = 'X'
		try:
			val = list( self.pdbRes ).index( resNum )
		except ValueError:
			print "Requested pdb residue %d not in chain %s" % (resNum, self.name)
			pass

		try:
			val = self.uniRes[val]
		except IndexError:
			print  "Requested pdb residue %d not in uniprot mapping of chain %s" % (val, self.name) 
			pass
		
		return val


	def uni2pdb(self, resNum):
		val = 'X'
		try:
			val = list( self.uniRes ).index( resNum )
		except ValueError:
			print "Requested uniprot residue %d not in uniprot entry %s" % (resNum, self.uniprot)
			pass

		try:
			val = self.pdbRes[val]
		except IndexError:
			print  "Requested uniprot residue %d not in pdb mapping of chain %s" % (val, self.name) 
			pass
		
		return val


########## Metal pair class  ##########
#Initialize with prody atom objects for metals and list of prody residue objects contacting
# Tihs class slims down info each atom carries to only store essentials (reduce filesize/RAM used)
class biMsite():
	def __init__(self, atom1, atom2, residue_matesList, pdb):
		self.name = "%s_%s+%s_%s" % ( atom1.getElement(), str( atom1.getSerial() ), atom2.getElement(), str( atom2.getSerial() ) )
		self.pdb = pdb
		self.contacts = []
		#For coordinating residues, loss of metal-res link, but bidentate should be store as 'A_ASP180-OD2+OD1'
		#self.contacts = [ '%s_%s%s-%s' % ( mate.getChid(), mate.getResname(), str( mate.getResnum() ), mate.getName() ) for mate in residue_matesList ] 
		for mate in residue_matesList:
			biFlg = 0
			fullID 	= '%s_%s%s-%s=%s' % ( mate.getChid(), mate.getResname(), str( mate.getResnum() ), mate.getName(), str(int(mate.getBeta()) ) )	 
			if len( self.contacts ) == 0:
				self.contacts.append( fullID )
			else:
				rBaseId = '%s_%s%s' % ( mate.getChid(), mate.getResname(), str( mate.getResnum() ) ) 

				for c in self.contacts:
					if c.split('-')[0] == rBaseId:
						new = c.split('=')[0] + '+' + mate.getName() + '=' + str(int(mate.getBeta()) ) + '=' + c.split('=')[-1]
						old = c
						biFlg +=1
				if biFlg > 0 :
					self.contacts.remove( old )
					self.contacts.append( new )
				else:
					self.contacts.append( fullID )


				

	def __repr__(self):
		return self.name




###### class end ######################



####################### Class definitions end #####################


############### Functions  ########################################

## Return the count of each amino acid found (includes non-natural) in a PDB model
def cntAA( pdbPath ):
	aaList = []
	print '*'
	pdbF = parsePDB( pdbPath )
	print '*'

	for res in pdbF.iterResidues():
		try:
			name = unNatAA[ res.getResname() ]
		except KeyError:
			name = res.getResname()
		if name in natAA.keys():
			aaList.append( name )

	final = Counter ( aaList )	

	return final

# determine the frequency of each amino acids in percent (e.g. 15.6) for a pdb database, given a list file with a path to each file
def freqAA (pathListFile):
	freq = Counter([])

	with open( pathListFile ) as file:
		for i in file:
			freq.update(   cntAA( i.rstrip() ) )
			print freq.items()
			print 

	freq2 = {}		
	for k,v in freq.items():
		val = 100 * float(v) / sum( freq.values() )
		freq2[k] = round(val, 2)

	print freq2.items()
	print sum( freq2.values() )/100.0

	return freq2

#freqAA( 'localZNdbFiles.txt' )


# Structural clustering using pre-user-defined symmetric distance matrix  
# Input a threshold level (presumably RMSD-cut off) such that all clusters are within threshold value

def heirarchy_cluster( matrix, threshold_RMSD = 2.0):
	import scipy.cluster.hierarchy as sp_clust

	# Complete linkage clustering
	link_matrix = sp_clust.complete( matrix )
	# Make flat clusters at the tree point where distance threshold is met: default RMSD < 2.0
	clusters 	= sp_clust.fcluster( link_matrix, threshold_RMSD, criterion = 'distance' )

	return clusters



def kmedoid_clustering( matrix, threshold_RMSD, lookupHash):
	from Bio.Cluster.cluster import kmedoids
	from collections import Counter
	import time

	cost = 4.0
	c = 135
	while cost > threshold_RMSD:
		trials	= np.size( matrix )
		start 	= time.time()
		clust 	= kmedoids( matrix , c, trials)
		print time.time()-start, 'elapsed for ', c, 'clusters', clust[-1], 'identical trajectory(s) from total',trials,
		# for this clustering scheme, calculate mean RMSD of each cluster
		clusterStats = {}
		index = 0
		# Look for each RMSD to the cluster centroid, store max within the cluster
		for p in clust[0]:
			#print p, index, matrix[p][index]
			try:
				clusterStats[p] = max( matrix[p][index], clusterStats[p] )
			except KeyError:
				clusterStats[p] = matrix[p][index]

			index += 1

		cost = max( clusterStats.values() )*4.0
		print 'maximum rmsd: ', cost
	print Counter( clust[0] )
	return clust



