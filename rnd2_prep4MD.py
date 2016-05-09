# rosetta to MD

# prep a pdb file for a equilibration and v short production run in explicit water
# This means solvating, 

import sys, os, subprocess as sp
from prody import *

# input 1: directory to put all work and input files in... makes new working dir for each trajecctory
# input 2: directory of input pdb files to  move and prep for MD
# input 3: path to quickPrep.tcl (Thomas wrote)
# input 4: path to setUp.tcl (Thomas wrote)
# input 5: path to template equilibration ACEMD script inputEQ
# input 6: path to production run script inputPROD (No changes, just copy)
# input 7: path to parameters file (just copt to wrking directory)

## REQUIRES:: VMD installation & alias to 'vmd' in .bashrc profile

# example command line
#
# python ~/bin/rnd2_prep4MD.py ~/peptideAmyloid/rnd2_design/MD_runs/ ~/peptideAmyloid/rnd2_design/MD_inputsptideAmyloid/test/quickPrep.tcl ~/peptideAmyloid/rnd2_design/tmp_work/setUp.tcl ~/peptideAmyloid/test/inputEQ ~/peptideAmyloid/test/inputPROD ~/peptideAmyloid/test/parameters 
#
#

## I/O for scripts to subcall

quickPrep 	= sys.argv[3]
setUp 		= sys.argv[4]
subScrEQ	= sys.argv[5]
subscrPROD	= sys.argv[6]
parameters	= sys.argv[7]

# read in EQ script once to edit each time
print subScrEQ
inEQ 		= open( subScrEQ, 'rU' ).readlines()


##



originalDir = os.getcwd()
# working dir for MD work
if not os.path.exists( sys.argv[1] ):
	os.mkdir( sys.argv[1] )

# work through input files one at a time
for f in os.listdir( sys.argv[2] ):

	if f[-3:] != 'pdb': continue

	# Make new working dir to put MD inputs/outputs
	oldPath 	= os.path.join( sys.argv[2], f )
	wrkSubDir	= os.path.join( sys.argv[1], os.path.splitext( f )[0] )
	if not os.path.exists( wrkSubDir ):
		os.mkdir( wrkSubDir )

	# clean up PDB file... No hydrogens, remove non-symmetrically repeated strands, rename things to CHARMM
	inPDB 		= parsePDB( oldPath, chain='XYZCDEFGH' ).select( 'not element H' ).copy()

	for k in inPDB.iterAtoms( ):

		if k.getChid() in [ 'C','D','E','F','G','H' ] and k.getResnum() == 13:
			k.setResname('HSP')
			continue
		if k.getResname() == 'ILE' and k.getName() == 'CD1':
			k.setName('CD')
			continue
		if k.getResname() == 'HID':
			k.setResname('HSD')
		if k.getResname() in [ 'HIE', 'HIS' ]:
			k.setResname('HSE')
		
	newPath = os.path.join( wrkSubDir, 'design.pdb' )
	writePDB( newPath, inPDB )

	## Move into working dir.... copy setup scripts into working dir
	
	os.chdir( wrkSubDir )

	sp.call( [ 'cp', quickPrep, wrkSubDir ] )

	sp.call( [ 'vmd', '-dispdev', 'text', '-e', setUp ] )

	## Write submission script from template

	# read in min-max data to print into inputEQ script
	mF 			= open( 'minmax.dat', 'rU' ).readlines()[0]
	x ,y, z 	= tuple( [ str( round( float( val ), 2) ) for val in mF.split() ] )
	# change template 
	newTxt 		= ''
	for line in inEQ:

		if line[:3] == 'set': 
			if line[4] == 'x':
				#print 'set x %s\n' % x 
				newTxt += 'set x %s\n' % x 
			elif line[4] == 'y':
				#print 'set y %s' % y 
				newTxt += 'set y %s\n' % y 
			elif line[4] == 'z':
				#print 'set z %s' % z 
				newTxt += 'set z %s\n' % z
			else:
				newTxt += line
		else:
			newTxt += line
	# print to file
	oFile = open( os.path.join( wrkSubDir, 'inputEQ' ), 'w' )
	oFile.write( newTxt )
	oFile.close()

	# copy parameters and prudction run script
	sp.call( [ 'cp', subscrPROD, wrkSubDir ] )
	sp.call( [ 'cp', parameters, wrkSubDir ] )

	# return to starting working directory... not really necessary 
	os.chdir( originalDir )

	#sys.exit()






'vmd -dispdev text -e setUp.tcl'
#find and replace cell size minmax.dat
'acemd inputEQ > logEQ'
'acemd inputPROD > logPROD'