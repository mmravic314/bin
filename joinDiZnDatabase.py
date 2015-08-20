#Marco Mravic UCSF Biophysics August 2015 DeGrado lab
# Quick script that takes a list (output from PISCES protein culling server) and a directory of PDB files
# makes directory and copies all PDB file on this list to that directory.

# python joinDiZnDatabase/py   ~/tertBuilding/cleanDiZn_nonRedund_res3-XRAY.txt   ~/pdb_080615/  ZN_db/

import sys, os, shutil

inF 	= open( sys.argv[1], 'rU' )
inF 	= inF.readlines()[1:]
if not os.path.exists( os.path.abspath( sys.argv[3] ) ):
		os.mkdir( os.path.abspath( sys.argv[3] ) )

for i in inF:
	oldPath = os.path.join( os.path.abspath( sys.argv[2] ), '%s.pdb' % i[:4] ) 
	newPath = os.path.join( os.path.abspath( sys.argv[3] ), '%s.pdb' % i[:4] )
	if not os.path.exists( newPath ):
		shutil.copy( oldPath, newPath )
