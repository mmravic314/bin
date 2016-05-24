import os, sys, numpy as np
from prody import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from collections import defaultdict

##
#python ~/bin/screen_TMpdbs.py ~/jan2016TMpdb/TM_redundant-FASTA/

import time

start = time.time()


matrix = matlist.blosum62
gap_open = -10
gap_extend = -2

def Align( s1, s2 ):
	return pairwise2.align.localxs( s1, s2, gap_open, gap_extend )[0]

#print aln_1+'\n'+aln_2
#print score, begin, end
#sys.exit()

# look through all TM_only pdbs/fastas and rank/sort all by pdb size (unique segments only)
# put all chains in memory with a dict within a dir
pLen = []
pdbDir = defaultdict( list )

for f in os.listdir(sys.argv[1]):
	
	if f[-6:] !='.fasta': continue

	seen 	= []
	tlen	= 0
	fpath 	= os.path.join( sys.argv[1], f )
	pdb 	= os.path.splitext( f ) [0]

	
	with open( fpath ) as fin:
		for i in fin:
			
			if i[0] =='>': 	continue


			if i not in seen:

				tlen += len(i)
				pdbDir[ pdb ].append(i)
				seen.append(i)

	pLen.append( ( tlen, pdb ) )


## store dir in pickle 

from itertools import product, combinations, permutations

dumped = []
print 'doing all pairwise fasta compairsons...'
for k in combinations( sorted( pLen, reverse=True), 2 ) :
	f1, f2, size1, size2 = k[0][1], k[1][1], k[0][0], k[1][0]

	if f1 in dumped: continue
	if f2 in dumped: continue




#	print f1, f2, size1, size2
	for s1 in pdbDir[f1]:
		for s2 in pdbDir[f2]:

			# dont compare short segments
			if min( len(s1), len(s2) ) < 8: continue
#			print s1, s2
#			print score( s1, s2 )
			aln_1, aln_2, score, begin, end = Align( s1, s2 )

			percID = score / min( len(s1), len(s2) )

#			print

#			print 
			if percID > 0.90: 
				#print score, min( len(s1), len(s2) ), percID
				#print aln_1 + '\n' + aln_2
				
				break

		if percID > 0.90: 

			if size2 > size1:
				dumped.append( f1 )
				#print 'rejected', f1
			else:
				dumped.append( f2 )
				#print 'rejected', f2
			break


	#print k 
	#sys.exit()
	#break

print time.time() - start, 's elapsed'
print 
good = []
# save
for pdb in pLen:

	if pdb[1] not in dumped:
		good.append( pdb[1] )
		print pdb[1]

print len( good ), 'total PDBs in non-redundant directory'

	# do pairwise comparison of each helix in list, quit and shu







	#break	 


