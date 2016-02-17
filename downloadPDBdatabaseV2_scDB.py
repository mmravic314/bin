# quick download pdbs
# From list of files like below:
# 3CAG 245

# Download whole pdb into input directory
# pull out chain with prody, by atom group. write atom group of htat chain to file
# Run a create PDS to make PDS file for data base
# also cath output PDB (parsed from createPDS) into another directory

import sys,os,subprocess as sp
from prody import *


with open( sys.argv[1] ) as file:
	for i in file:
		call, pdb, chain = i.split()[0], i[:4], i.split()[0][5:]

		targDir = sys.argv[2]
		pdsDir	= sys.argv[3]
		oPdbDir	= sys.argv[4]


		sp.call( [ 'wget', 'www.rcsb.org/pdb/files/%s.pdb' % (pdb) , '-P', targDir])

		pdbPath = os.path.join( targDir, '%s.pdb' % ( pdb ) )
		

		## Some subunits of HUGE proteins have different file format. JUST PASS!
		try:
			subset 	= parsePDB( pdbPath ).select( 'chain %s' % (chain) )
		except AttributeError:
			os.remove( pdbPath )
			continue


		newPath = os.path.join( targDir, call + '.pdb' )
		writePDB( newPath, subset )
		os.remove( pdbPath )

		sp.call( [ '/home/xray/termanal/createPDS', '--type', 'target', '--pdb', newPath, '--pds', os.path.join( pdsDir, '%s.pds' % ( call ) ) , '--opdb', os.path.join(oPdbDir, '%s.pdb' % ( call ) ) ] )
		#cmd = [ '/home/xray/termanal/createPDS', '--type', 'query'] #'--pdb', os.path.join( targDir, '%s.pdb' % ( pdb ) ), '--opdb', os.path.join( targDir, '%s_.pdb' % ( pdb ) ) ]	

		#sp.call( cmd )
