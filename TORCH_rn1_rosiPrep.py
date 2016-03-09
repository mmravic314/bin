##  given input list of good models
#   unzip & move PDBs their own folder
#   Write resfile for system
import sys, os, shutil, gzip

if not os.path.exists( sys.argv[2] ):
	os.mkdir( sys.argv[2] )



with open( sys.argv[1] ) as file:
	for i in file:

		if i[0] == '#': continue

		oldPath = os.path.join( sys.argv[3], 'model_%s.pdb.gz' % (i.split()[0] ) )
 		
 		newDir 	= os.path.join( sys.argv[2], 'model_%s/' % (i.split()[0] ) )
 		
 		newPath = os.path.join( sys.argv[2], 'model_%s.pdb' % (i.split()[0] ) )

		with gzip.open(oldPath, 'r') as f_in, open(newPath, 'w') as f_out:
			shutil.copyfileobj(f_in, f_out)


