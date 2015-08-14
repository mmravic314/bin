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


# mapped PDB chain containing which pdb residue corresponds to each uniprot entry from DBREF lines 
# Also holds uniprot with download option
class PdbChain():
	# does mapping at initialization given a DBREF line from a PDB file
	def __init__(self, pdbLine):
		self.name = pdbLine[7:11] + pdbLine[12].lower() 
		self.uniprot = pdbLine[33:39]
		self.pdbRes = np.arange( int( pdbLine[14:18].strip() ) , int ( pdbLine[19:24].strip() ) + 1 )
		self.uniRes  = np.arange( int( pdbLine[55:60].strip() ) , int ( pdbLine[62:68].strip() ) + 1 )
		self.pdbMetalCoords = []

	def __repr__(self):
		return '%s, %s' % (self.name, self.uniprot) 

	def pdb2uni(self, resNum):

		try:
			val = list( self.pdbRes ).index( resNum )
		except ValueError:
			print "Requested pdb residue %d not in chain %s" % (resNum, self.name)
			sys.exit()

		try:
			val = self.uniRes[val]
		except IndexError:
			print  "Requested pdb residue %d not in uniprot mapping of chain %s" % (val, self.name) 
			sys.exit()
		
		return val


	def uni2pdb(self, resNum):
		
		try:
			val = list( self.uniRes ).index( resNum )
		except ValueError:
			print "Requested uniprot residue %d not in uniprot entry %s" % (resNum, self.uniprot)
			sys.exit()

		try:
			val = self.pdbRes[val]
		except IndexError:
			print  "Requested uniprot residue %d not in pdb mapping of chain %s" % (val, self.name) 
			sys.exit()
		
		return val
####################### Class definitions end #####################
