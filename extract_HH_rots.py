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

## preparation for counting rotamers per residue later on
rotamer_Cnts = {}
for k in resi:
	rotamer_Cnts[k] = 0
## 

## defining the renaming scheme to be used later & write a re-formatted pdb file
segs 	= set( resiHash.values() )
ChASCII = 64 
prvSeg	= ''
reName 	= {} 
for r in inPDB.iterResidues():
	seg 	= r.getSegname()
	resID 	= '%s,%d' % ( r.getChid(), r.getResnum() )

	if seg != prvSeg:
		resN 	= 1000
		ChASCII += 1
	r.setResnum(resN)
	ch = chr( ChASCII )
	r.setChids( ch )
	reName[ resID ] = '%s,%d' % ( ch, resN )

	resN += 250
	prvSeg = seg

##

### extract rotamers for
txt 	= ''
prvID 	= '' 
rotNum 	= 0
with open( rotFile ) as fin:
	for i in fin:

		if i[:4] == 'ATOM':

			####
			# NOTE: renumber each residue and helix 
			# do NOT include the backbone atoms of the side chains 
			# do NOT include side-chain Hydrogens... 
			# DO ewrite the pdb lines to include the element...
			####

			## gather info
			resID 	= '%s,%s' % ( i[21] , i[22:26].strip() )
			name 	= i[17:20]
			atom 	= i[12:16].strip()
			newBase = reName[resID]
			
			if atom == 'N': 

				# at each new residues, calculate ID
				if resID != prvID:
					newBase = reName[resID]
					print newBase


				prvID = resID

				# optional section for counting rotamers total
				try:
					rotamer_Cnts[resID] += 1
				except KeyError: 
					pass
				# 

			if name == 'ALA':
				txt += i
				print i.rstrip()

				continue


			continue



			if resID in resi and name in eligible:
				txt += i

sys.exit()

### Done extracting rotamers.... now write the big PDB file

oPath = os.path.join( os.path.dirname( sys.argv[3] ), 'tmp2_' + os.path.basename( sys.argv[3] ) )
oFile = open( oPath, 'w' )
oFile.write( txt )
print sorted( resi )

cnt 	= np.array( rotamer_Cnts.values() )
print 'rotamers per resi (mean & array):', np.mean(cnt), sorted( cnt )  

print 'time elapsed (s):', time.time() - ts

###
