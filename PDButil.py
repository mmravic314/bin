## PDB file utilities ##
##################### Written By: Marco Mravic  ###  Degrado Lab UCSF Biophysics Aug 2015

import numpy as np

divalent_metals = [
'BA',
'BE', 
'CD', 
'CA', 
'CR', 
'CO', 
'CU', 
'EU', 
'GD', 
'GE', 
'FE', 
'LA', 
'PD', 
'MG', 
'MN', 
'HG', 
'NI', 
'OS', 
'PT', 
'RU', 
'SR', 
'SN',
'U', 
'V', 
'Y', 
'ZN' ]

natAA = {}
natAA["ALA"] = 'A'; natAA["CYS"] = 'C'; natAA["ASP"] = 'D'; natAA["GLU"] = 'E'; natAA["PHE"] = 'F';
natAA["GLY"] = 'G'; natAA["HIS"] = 'H'; natAA["HSD"] = 'H'; natAA["ILE"] = 'I'; natAA["LYS"] = 'K';
natAA["LEU"] = 'L'; natAA["MET"] = 'M'; natAA["ASN"] = 'N'; natAA["PRO"] = 'P'; natAA["GLN"] = 'N';
natAA["ARG"] = 'R'; natAA["SER"] = 'S'; natAA["THR"] = 'T'; natAA["VAL"] = 'V'; natAA["TRP"] = 'W'; natAA["TYR"] = 'Y';

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


####################### Class definitions end #####################
