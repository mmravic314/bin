# Use MSL/MASTER/TERMs contact finder to build all rotamers for 3 helix bundles

'''
# on directory of 3 helix bundle pdb files, with unique segnames per helix
xray@dgl-xray:~/splayBundle/6_helix/TM_3helix$ for i in *.pdb; do ~/termanal/confind --p $i --o ./${i::-4}.cmap --rLib ~/termanal/support.default/rotlib/RR2000.rotlib --rout ./${i::-4}.rout ; done

then this script (in ~/splayBundle/6_helix/TM_3helix/  )
> python ~/bin/extract_HH_rots.py 1h2s_Af-Ag-BbCLEAN.pdb 1h2s_Af-Ag-Bb.cmap 1h2s_Af-Ag-Bb.rout

'''

# Now use this script to extract uncharged H-bonding rotamers only (.rout file), and only from residues having contacts (.cmap file)

import sys, os, numpy as np, time
from prody import * 
from collections import defaultdict

#from PDButil import NatAA

# I/O and global definition of rotamers searching for 
pdb 		= sys.argv[1]
cmap 		= sys.argv[2]
rotFile 	= sys.argv[3]

eligible 	= [ 'SER', 'THR', 'TYR', 'TRP', 'ASN', 'GLN', 'HIS', 'GLY' ]

ts = time.time()
#



### determine which residues are in each segment... to allow inter-helical contacts to be identified
resiHash = {}
inPDB = parsePDB( pdb ) # Prody atom group object

for r in inPDB.iterResidues():
	resiHash[ '%s,%d' % ( r.getChid(), r.getResnum() ) ] = r.getSegname()

###

### read in cmap file & define list of residues to pull rotamers for
## Include Ala to get entire structure's backbone coordinates
resi 	= []
cutoff 	= 0.03		# cmap rotemer interaction probability for contacts

#byRes 	= defaultdict( list ) 
with open( cmap ) as fin: 
	for l in fin:
		if l[:7] == 'contact':
			cnt = l.split()
			r1 	= cnt[1]
			r2 	= cnt[2]
			freq= float( cnt[3] )

			if resiHash[r1] != resiHash[r2] and freq > cutoff:
				if not r1 in resi:
					resi.append( r1 )
				if not r2 in resi:
					resi.append( r2 )
###


### extract rotamers for
rotamer_Cnts = {}
for k in resi:
	rotamer_Cnts[k] = 0 


txt = ''
with open( rotFile ) as fin:
	for i in fin:

		if i[:3] == 'REM':
			info 			= i.split()
			resID, name 	= info[1][:-5], info[1][-4:-1] 
			if name == 'ALA':
				txt += i
				continue

			if resID in resi and name in eligible:
				txt += i

		elif i[:3] == 'TER':
			info 			= i.split()
			resID, name 	= ','.join ( info[-2:] ), info[2] 
			if name == 'ALA':
				txt += i
				continue

			if resID in resi and name in eligible:
				txt += i + 'END\n'

		elif i[:3] == 'END': continue

		else:

			####
			# NOTE: it might be smart to include all GLY for all residues, as the representative backbone atoms 
			# and not include the backbone atoms of the side chains 
			# and not include side-chain Hydrogens... 
			# and rewrite the pdb lines to include the element...
			####

			resID 	= '%s,%s' % ( i[21] , i[22:26].strip() )
			name 	= i[17:20]
			atom 	= i[12:16].strip()

			# optional section for counting rotamers total
#			if atom == 'CA': 
#				try:
#					rotamer_Cnts[resID] += 1
#				except KeyError: 
#					pass
			# 

			if name == 'ALA':
				txt += i
				continue

			if resID in resi and name in eligible:
				txt += i

oPath = os.path.join( os.path.dirname( sys.argv[3] ), 'tmp_' + os.path.basename( sys.argv[3] ) )
oFile = open( oPath, 'w' )
oFile.write( txt )
print sorted( resi )

cnt 	= np.array( rotamer_Cnts.values() )
print 'rotamers per resi (mean & array):', np.mean(cnt), sorted( cnt )  

print 'time elapsed (s):', time.time() - ts

###
