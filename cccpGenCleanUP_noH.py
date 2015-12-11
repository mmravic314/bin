####### INPUT FULL/ABSOLUTE PATH TO DIRECTORY OF CCCP PDB's TO MODIFY

################# Clean up all-ala coiled-coil pdb files ##################

import os, sys, numpy as np

## reformat alanine coiled-coils from cccp gen (before hydrogens are added)
def rosCCprep( lines ):
	elements 	= [ 'C', 'N', 'O', 'H', ]
	cleanStr 	= ''


	for l in lines:
		#print l.rstrip()
		# IGNORE non-'ATOM' lines and second conformations 
		if l[:3] 	== 'TER':

			outStr 	=  'TER{:>8}{:>9}{:>2}{:>4}\n'.format( l[6:11], "GLY", l[67:].strip(), str(resi) )
			cleanStr+= outStr
			chInd	+= 1  
			continue

		if l[:6] == 'REMARK':
			cleanStr+=l

		if l[:4] != 'ATOM':
			continue

		if l[12:16].strip() == 'CB': continue

		
		lineList = [
				l[0:6],
				l[6:11],
				l[12:16],
				' ',
				"GLY",
				l[67:].strip(),
				l[22:26],
				' ',
				l[30:38],
				l[38:46],
				l[46:54],
				l[54:60],
				l[60:66],
				[x for x in l[12:16] if x in elements ][0], 
				'    '
				]

		outStr = '{:<6}{:>5} {:<4}{:<1}{:<3} {:<1}{:>4}{:<1}   {:>8}{:>8}{:>8}{:>6}{:>6}          {:>2}{:>2}\n'.format( *lineList )
		#print outStr
		cleanStr+= outStr

	return cleanStr+'END'



###################   MAIN   ##############################


for f in os.listdir( sys.argv[1] ):
	if f[0] != '0': continue

	oldPath = os.path.join( sys.argv[1], f ) 
	## Handle path of directory whether ends in '/' or not (dumb way, there is probably an app for that)
	if sys.argv[1].strip()[-1] == '/':
		dirname = os.path.split( sys.argv[1].strip()[:-1] )[-1]
	else:
		dirname = os.path.split( sys.argv[1].strip() )[-1]

	newPath = os.path.join( sys.argv[1], '%s_%s.pdb' % ( dirname , f.split('.')[0] ) )


	inFile 		= open( oldPath, 'rU' )
	outStr 		= rosCCprep( inFile.readlines() )
	inFile.close( )

	newfile		= open( newPath, 'w' )
	newfile.write( outStr )
	newfile.close()
	

	