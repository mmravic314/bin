#!/usr/bin/env python2

#$ -S /usr/bin/python
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G
#$ -l h_rt=00:30:00
#$ -cwd
#$ -j y


import os, sys, subprocess as sp
from prody import * 

## Marco Mravic DeGrado lab Dec 2015 UCSF 
## input 1: FULL/ABSOLUTE PATH of each file
## input 2: path to directory with both createPDS and master binaries
## input 3: path to text file list of pathes to TM-helix-only pds (database)
## input 4: path to text list to append 4 column summary of (path & number of matches below 1.5 A to each fragment/interface)

## SAMPLE
## python ~/bin/breakUpInterfaces.py ~/splayBundle/tmBBsplay/1m0l_1E_bb/1m0l_1E_bb_01000.pdb ~/termanal/support.default/TM_Database/list_tm.txt


## make 3 pdb files of 12 residue fragments from each interface in 4 helix coiled coils: A-B, B-C, C-D
## convert to PDS and submit each to master search of transmembrane data base
## record number of matches under 2 A for each fragment in .match files, in a single match directory
selDict = { 
			'AB':'chain A B resnum 14 to 25',
			'BC':'chain B C resnum 14 to 25',
			'CD':'chain C D resnum 14 to 25'
			}


def matchCnt( lineList ):
	cnt = 0

	return cnt

pdb 			= parsePDB( sys.argv[1] )
filename 		= os.path.basename( sys.argv[1] )
outList 		= [ sys.argv[1] ]
### for each pdb fragment (aka helix-helix interface)
### 1) write a temp file, 2) convert to pds, 3) search master, 4) record matches, 5) remove temp files
for k,v in selDict.items():

		# 1) write temp file of pdb fragment
		newPath 	= os.path.join( os.path.dirname( sys.argv[1] ) , '_'.join( [ k ,filename.split('_')[0], filename.split('_')[1], filename.split('_')[-1] ] ) ) 
		writePDB( newPath, pdb.select( v )  )

		# 2) convert to pds
		pdsPath		= newPath[:-4] + '.pds'
		cPDSpath 	= os.path.join( sys.argv[2], 'createPDS' )
		cPDScmd 	= [ cPDSpath, '--type', 'query', '--pdb', newPath  ]
		print cPDScmd
		#sp.call( cPDScmd )

		# submit master, putting output results into file
		masterPath 	= os.path.join( sys.argv[2], 'createPDS' )
		pdsPath		= newPath[:-4] + '.pds'
		matchOut 	= newPath[:-4] + '.match'
		masterCmd 	= [ masterPath, '--query', pdsPath, '--targetList', sys.argv[3], '--matchOut',  ]

		# sleep until finish and then delete tmp files
		while not os.path.exists( matchOut ):
			time.sleep( 30 )
		print "DONE WITH MATCHING"
		os.remove( pdsPath )
		os.remove( newPath )

		# summarize match data into 3 column text file
		mFile = open( matchOut, 'rU' ).readlines()
		outList.append( matchCnt( mFile ) ) 

		print newPath
		print outList



		sys.exit()
print outList
	