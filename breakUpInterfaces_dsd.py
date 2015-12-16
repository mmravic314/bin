
import os, sys, subprocess as sp, time
from prody import * 
## Marco Mravic DeGrado lab Dec 2015 UCSF 
## input 1: FULL/ABSOLUTE PATH to directory containing backbone files
## input 2: path to directory with both createPDS and master binaries
## input 3: path to text file list of pathes to pds database to search
## input 4: path to text list to append 4 column summary of (path & number of matches below 1 A to each fragment/interface)

## SAMPLE
## > python ~/bin/breakUpInterfaces_dsd.py  ~/tertBuilding/DSD_ligation/dsd_hetero/ ~/termanal/ ~/termanal/support.default/database/list designScore.txt

## make pdb files of e helix residue frag at core
## convert to PDS and submit each to master search of transmembrane data base
## record number of matches under 2 A for each fragment in .match files, in a single match directory
selStr = 'chain A B C resnum 10 to 19'

# Summary string to append recordings
sumStr = ''


for filename in os.listdir( sys.argv[1] ):

	# enter each reformatted pdb file
	if filename[0] == '0' or filename[-4:] != '.pdb': continue
	
	try:
		pdb 			= parsePDB(  os.path.join( sys.argv[1], filename ) )
	except:
		continue


	print filename

### 1) write a temp file, 2) convert to pds, 3) search master, 4) record matches, 5) remove temp files


	# 1) write temp file of pdb fragment
	newPath 	= os.path.join(  sys.argv[1] , filename[:-4] + '_frag.pdb' )
		
	print newPath


	try: 
			writePDB( newPath, pdb.select( selStr )  )
	except:
			break
	# 2) convert to pds
	cPDSpath 	= os.path.join( sys.argv[2], 'createPDS' )
	cPDScmd 	= [ cPDSpath, '--type', 'query', '--pdb', newPath  ]
	sp.call( cPDScmd )

	# submit master, putting output results into file
	masterPath 	= os.path.join( sys.argv[2], 'master' )
	pdsPath		= newPath[:-4] + '.pds'
	matchOut 	= newPath[:-4] + '.match'
	masterCmd 	= [ masterPath, '--query', pdsPath, '--targetList', sys.argv[3], '--matchOut', matchOut,
						'--rmsdCut', '1.5'
					  ]
	
	sp.call( masterCmd )
	time.sleep( 0.1 )

	# summarize match data into 3 column text file, delete intermediate files
	numMatch	 = str( len( open( matchOut, 'rU' ).readlines() ) )

	time.sleep( 0.1 )
	os.remove( pdsPath )
	os.remove( newPath )
	os.remove( matchOut )


	sumStr += os.path.join( sys.argv[1], filename ) + str( numMatch ) + '\n' 
print sumStr
# Summary file to append recordings
sumFile = open( sys.argv[4], 'a' )	
sumFile.write( sumStr )
sumFile.close()
