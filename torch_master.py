#!/bin/python

# given input dir full of fragment file, check if appropriate number of files (if not, abort and print null score file)
# make pds for each
# do master search 
# print score for each score to file

# maybe take some pairwise sequence stats, given residue identities on the sheet with output sequence file

import sys, os, subprocess as sp, shutil, time
from PDButil import UnNatAA			# this is just a three to one residue identitiy converter than handles unnatural/modified amino acids

# example command line
# python ~/bin/torch_master.py ~/peptideAmyloid/tmp_model_1/ ~/peptideAmyloid/scores_rnd1/ ~/termanal/ ~/termanal/support.default/152901_bc_30-scPDB/list.txt
###

start = time.time()

fragsDir	= sys.argv[1]


if fragsDir[-1] == '/':
	fragsDir  = fragsDir[:-1]

scorePath 	= os.path.join( sys.argv[2], fragsDir.split('/')[-1][4:] + '.sc' )

createPDSp	= os.path.join( sys.argv[3], 'createPDS' 	)
masterp		= os.path.join( sys.argv[3], 'master' 		)
dbList 		= sys.argv[4]

# Jump out if all files not find
if len( [ x for x in os.listdir( sys.argv[1] ) if x[-3:] == 'pdb' ] ) < 4:

	print 'not enough fragments made, abort!'
	shutil.rmtree( fragsDir )

	txt = 'HH SH1 SH2 SH3\n0 0 0 0 0'
	oFile = open(scorePath, 'w')
	oFile.write( txt  )

	sys.exit()


txt = 'HH SH0 SH1 SH2\n'
for f in sorted( os.listdir( fragsDir ) ):

	if f[-3:] != 'pdb': continue

	
	if 'SH' in f:
		rmsd = '1.1'
	else:
		rmsd = '0.85'

	fragPath 	= 	os.path.join( fragsDir, f )
	pdsPath 	=	os.path.join( fragsDir, f[:-1] + 's' )
	seqOutPath	=	fragPath[:-4] + '.seq'

	cmdPDS 	= [ createPDSp, '--type', 'query', '--pdb', fragPath ] 

	cmdMstr	= [ masterp, '--query', pdsPath, '--targetList', dbList, '--bbRMSD', '--rmsdCut', rmsd, '--seqOut', seqOutPath ]

	sp.call( cmdPDS )
	sp.call( cmdMstr )

	# search files for unique sequences  n*(n-1) complexity
	uniq = []
	with open( seqOutPath ) as file:
		for i in file:
			if float( i.split()[0] ) > float( rmsd ): break

			seq = ''.join( [ UnNatAA[x] for x in i.split()[2:] ] )
			if seq not in uniq:
				uniq.append( seq )
	txt += '%d ' % len( uniq )


print 'SUCESSFUL wall clock time used (s)', time.time() - start
print txt
print scorePath
oFile = open(scorePath, 'w')
oFile.write( txt  )
#shutil.rmtree( fragsDir )
