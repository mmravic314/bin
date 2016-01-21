#Marco Mravic UCSF Biophysics

# input .seq file and output text multiple sequence alignment for that file
# and remember the residue length of each helix in this design. 
import sys, os
from PDButil import UnNatAA

# loop through seq files  in DEC15 dir
for i in os.listdir( sys.argv[1] ):

	path 	= os.path.join( sys.argv[1], i )
	oPath 	= os.path.join( sys.argv[1], i[:-4] + '.mseq' )
	seqList = []
	singles = []
	if i[-4:] == 'mseq': continue
	with open( path ) as file:	
		for l in file:
			H1 = ''.join( [ UnNatAA[x] for x in l.split()[2:14]] )
			H2 = ''.join( [ UnNatAA[x] for x in l.split()[14: ]] )
			seq= H1 + '-' + H2
			print l

			# prevent symmetric matches being counted twice
			if H1 in singles: continue
			if H2 in singles: continue

			if seq not in seqList:
				seqList.extend( ['>%s' % ( l.split()[1] ), seq ] )

			# add every sequence to a list so chain A of match one isnt counted for chain B in match 2
			if H1 in singles: continue
			if H2 in singles: continue



	oFile = open( oPath, 'w' )
	for s in seqList:
		oFile.write( s + '\n' )
		print s
		
	
	print 
	print