# Marco Mravic UCSF DeGrado lab Sep 2015
# for a list of PDB files and a path to the directory that pdb file is in, print the PDB ID and title lines from file

import sys,os

with open(sys.argv[1]) as file:
	for p in file:
		path 	=  os.path.join( os.path.abspath( sys.argv[2] ) , '%s.pdb'%(p.rstrip()) )
		pF 		=  open( path, 'rU' )
		print p
		for l in pF:
			if 'TITLE' in l or 'HEADER' in l:
				print l
			if 'REMARK' in l:
				break 
		print '\n'
		pF.close() 


