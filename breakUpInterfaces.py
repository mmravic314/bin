#!/usr/bin/env python2
                 
#$ -S /usr/bin/python                    #-- the shell for the job
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=1G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=1G,scratch=1G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)



import os, sys, subprocess as sp, time
from prody import * 
## Marco Mravic DeGrado lab Dec 2015 UCSF 
## input 1: FULL/ABSOLUTE PATH to directory containing backbone files
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

# Summary string to append recordings
sumStr = ''



## SKIP ALREADY DONE DIRECTORY 
### Note this is ghetto way to parallelize using 4 cores locally
lists = [
'/home/xray/splayBundle/tmBBsplay/1v55_1D_bb/',
'/home/xray/splayBundle/tmBBsplay/1v55_6C_bb/',
'/home/xray/splayBundle/tmBBsplay/4dji_1C_bb/',
'/home/xray/splayBundle/tmBBsplay/4gc0_1B_bb/'
]


if sys.argv[1] not in lists:
	sys.exit()

#if sys.argv[1] in [ lists[0],  lists[5],  lists[6],  lists[7],  lists[8] ]:
#	sys.exit()


#### End ghetto section #############


for filename in os.listdir( sys.argv[1] ):

	# enter each reformatted pdb file
	if filename[0] == '0' or filename[-4:] != '.pdb' or filename[:2] in selDict.keys(): continue
	
	try:
		pdb 			= parsePDB(  os.path.join( sys.argv[1], filename ) )
	except:
		continue

	outList 		= [ os.path.join( sys.argv[1], filename ) ]

	print filename
### for each pdb fragment (aka helix-helix interface)
### 1) write a temp file, 2) convert to pds, 3) search master, 4) record matches, 5) remove temp files
	for k,v in selDict.items():

		# 1) write temp file of pdb fragment
		newPath 	= os.path.join(  sys.argv[1] , '_'.join( [ k ,filename.split('_')[0], filename.split('_')[1], filename.split('_')[-1] ] ) ) 
		
		print newPath
		try: 
			writePDB( newPath, pdb.select( v )  )
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
						'--rmsdCut', '1'
					  ]
	
		sp.call( masterCmd )
		time.sleep( 0.1 )

		# summarize match data into 3 column text file, delete intermediate files
		numMatch	 = str( len( open( matchOut, 'rU' ).readlines() ) )
		outList.append( numMatch ) 
		time.sleep( 0.1 )
		os.remove( pdsPath )
		os.remove( newPath )
		os.remove( matchOut )


	sumStr += ' '.join( outList )  + '\n' 

# Summary file to append recordings
sumFile = open( sys.argv[4], 'a' )	
sumFile.write( sumStr )
sumFile.close()
