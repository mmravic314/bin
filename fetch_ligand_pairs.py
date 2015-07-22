### Marco Mravic UCSF DeGrado Lab July 2015 

## Heavy script that reads in a directory of pdb files with ligands
#  Blasts the sequence of each protein and find all pdb entries with 100% identity
#  Downloads them all and outputs a list grouping pdb entries

# looking to find pairs/groups of structures with the same sequence: w/ ligand & apo  


import sys,os


##
# Read .csv file and print list of eligible pdbs
pdbList = []
bioLigandF = open( sys.argv[1], 'r')
for i in bioLigandF:

	#Ignore peptide and nucleic acid ligands
	line = i.split(',')
	if line[4].strip('\"') == 'III' or line[4].strip('\"') == 'NUC':
		continue

	if line[1].strip('\"')  not in pdbList:
		pdbList.append( line[1].strip('\"') )
		print line[1].strip('\"'), line







