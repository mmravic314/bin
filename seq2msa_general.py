#Marco Mravic UCSF Biophysics

# input .seq file and output text multiple sequence alignment for that file
# and remember the residue length of each helix in this design. 
import sys, os
from PDButil import UnNatAA

recall = []
# loop through seq files  in DEC15 dir
with open( sys.argv[1] ) as fin:
	for l in fin:
#		print l

		seq 	= ''.join( [ UnNatAA[ x ] for x in l.split()[1:] ] )
		rmsd 	= l.split()[0]
		
		# skip redundants
		if seq in recall: continue
		recall.append( seq )

		print 	'>%s\n%s\n' % ( rmsd, seq )

#print len( recall )