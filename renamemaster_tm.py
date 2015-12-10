
import os,sys

listFile = open( os.path.join( sys.argv[1], 'list_tm.txt' ), 'w' )
for f in os.listdir( sys.argv[1] ):
		print f
		if f[-4:] == '.pds':
			listFile.write( os.path.join( sys.argv[1], f ) +'\n')
