# parse rosettag output file, dump sequence and filter stats (hard coded here)

import sys, os
from PDButil import UnNatAA		## mmravic local module
from prody import *
#python ~/bin/parseRosiOutCMP.py ~/tertBuilding/CMP_bobo/redesign_inputs/outputs/ > MMredesignSeqs.txt

#packstat 0.66498
#score_eval -18.3668
#uhb 4

designs = {}

for i in os.listdir( sys.argv[1] ):
	if i[-3:] != 'pdb': continue

	path = os.path.join( sys.argv[1], i )

	# find sequence
	seq = ''
	inPDB = parsePDB( path, subset = 'ca' )
	for a in inPDB.iterAtoms():
		seq += UnNatAA[ a.getResname() ]    ## 3 to 1 residue ID code conversion

	# change sequence to the single peptide representation
	seq = seq[9:] + 'C' + seq[:9]


	# pull stats from file lines
	with open( path ) as file:
		for l in file:
			if l[:4] == 'ATOM': continue

			if len(l) < 1: continue

			first = l.split()[0]

			if first == 'packstat':
				ps = l.split()[1]

			elif first == 'score_eval':
				sc = l.split()[1]

			elif first == 'uhb':
				uhb= l.split()[1]

			else: continue

	#print i, seq, sc, ps, uhb

	try: 
		if designs[seq][0] < ps:
			designs[seq] = [ps, sc, uhb, i]

	except KeyError:
		designs[seq] = [ps, sc, uhb, i]

#	sys.exit()
print
for k, v in designs.items():
	print '>', v[0], v[1], v[2], v[3]
	print k
	print