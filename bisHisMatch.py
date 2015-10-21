##############################################################
## Marco Mravic DeGrado Lab Sept 2015
## From bis Zn site with 2 HIS each in bacterial PTE, 
## write rotamers & their rotations satisfying geometric constraints to Bis Zn with activated hydroxyl
# Want 1 PDB with each of 4 target HID residues, indexed 0000 to 7777 for each combo of 8 rotamers
# Then rotate each HID around NE-ZN axis, indexed by 40 degree increments 1111 to 9999
# complexity:  8^4 * 9^4  
##############################################################

# python ~/bin/bisHisMatch.py ~/tertBuilding/His_library/HSD-Richardson_rotamersFinal.pdb ~/tertBuilding/paraoxon-bisZNOH.pdb ~/tertBuilding/paraoxon-bisZNOHBare.pdb  ~/tertBuilding/bisHisPossibleSites/



import sys, os, cPickle as pic, numpy as np, gzip
from itertools import product
from prody import *
from math import *
from numpy.linalg import norm

def rotation_matrix( axis, angle ):
    """Generate the rotation matrix from the axis-angle notation.

    Conversion equations
    ====================

    From Wikipedia (http://en.wikipedia.org/wiki/Rotation_matrix), the conversion is given by::

        c = cos(angle); s = sin(angle); C = 1-c
        xs = x*s;   ys = y*s;   zs = z*s
        xC = x*C;   yC = y*C;   zC = z*C
        xyC = x*yC; yzC = y*zC; zxC = z*xC
        [ x*xC+c   xyC-zs   zxC+ys ]
        [ xyC+zs   y*yC+c   yzC-xs ]
        [ zxC-ys   yzC+xs   z*zC+c ]


    @param matrix:  The 3x3 rotation matrix to update.
    @type matrix:   3x3 numpy array
    @param axis:    The 3D rotation axis.
    @type axis:     numpy array, len 3
    @param angle:   The rotation angle.
    @type angle:    float
    """

    # Trig factors.
    ca = cos(angle)
    sa = sin(angle)
    C = 1 - ca

    # Depack the axis.
    x, y, z = tuple( axis )

    # Multiplications (to remove duplicate calculations).
    xs = x*sa
    ys = y*sa
    zs = z*sa
    xC = x*C
    yC = y*C
    zC = z*C
    xyC = x*yC
    yzC = y*zC
    zxC = z*xC

    # Update the rotation matrix.
    matrix  	 = np.zeros( (3,3) )
    matrix[0, 0] = x*xC + ca
    matrix[0, 1] = xyC - zs
    matrix[0, 2] = zxC + ys
    matrix[1, 0] = xyC + zs
    matrix[1, 1] = y*yC + ca
    matrix[1, 2] = yzC - xs
    matrix[2, 0] = zxC - ys
    matrix[2, 1] = yzC + xs
    matrix[2, 2] = z*zC + ca
    return matrix



#Naming for each his coordinating nitrogen:
# G1 (coming 90 degrees to ZNa-ZNb plane) below plane, to Zna 
# H1 Above the ZN-ZN-OH, coordinated (~120 degrees from A1 via ZNa) trigonal bipyrimdal
# I1 (coming 90 degrees to ZNa-ZNb plane) below plane, to Znb
# J1 Above the ZN-ZN-OH, coordinated (~100 degrees from B1 via ZNa) distorted Tetrahedral
### transform each HID and HIE rotamer to each of the 4 nitrogens 
##  Must place the ring such that it is in a reasonable angle to the ZN - N vector (160 degrees between )
##  For each rotamer, rotate each 33 degrees at a time

#ProDy format selection string for epsilon coordinating HID tautomer
selStrHID 		= [ 'chain G resnum 1', 'chain H resnum 1', 
		'chain I resnum 1', 'chain J resnum 1' ]

selStrMobile 	= [ 'chain G resnum 1', 'chain H resnum 1', 
		'chain I resnum 1', 'chain J resnum 1', 'resname EPL' ]

selStrNE2 		= [ 'chain G resnum 1 name  NE2', 'chain H resnum 1 name NE2', 
		'chain I resnum 1 name  NE2', 'chain J resnum 1 name NE2' ]

selStrRing 		= [ 'chain G resnum 1 name  CD2 CE1 NE2', 'chain H resnum 1 name CD2 CE1 NE2', 
		'chain I resnum 1 name   CD2 CE1 NE2', 'chain J resnum 1 name  CD2 CE1 NE2' ]

#Remember which Zn is coordinated to which HIS, by chain
znHash = { 'G':'G', 'H':'G', 'I':'I', 'J':'I' }

# CHARMM topology file syntax atom names
hsdCharmm = [
'N', 
'CA',
'C',
'O',
'CB',
'CG',
'CD2',
'ND1',
'CE1',
'NE2',
'HN',
'HA',
'HB1',
'HB2',
'HD2',
'HD1',
'HE1' ] 

# read in pdb file with multiple models of each HIS rotamer (Lovell et al 2000 Richardson)
# Rename atom names to match charmm topology
hisObjs = [ parsePDB(sys.argv[1], model = x) for x in np.arange(1,9) ]
for h in hisObjs:
	h.setResnames( [ 'HIS' for x in h ] )
	h.setNames( hsdCharmm )

## Loading in pdb bis his site ('Core'), and site stripped of HIS residues ('Bare')
## transform coordinates to center HOH at origin
znCore 	= parsePDB( sys.argv[2] )
znBare 	= parsePDB( sys.argv[3] )
moveAtoms( znCore.select( 'resname OH name O1'), to = np.array( [0,0,0] ) , ag = True )
moveAtoms( znBare.select( 'resname OH name O1'), to = np.array( [0,0,0] ) , ag = True )


# Quick hashing of vector between NE and coordinating Zn2+ indexed by his selection string
# if that zn has coords of (0,0,0). Store this vector and the NE2's coordinates when HOH is at origin
vectorAxes = {}
for NE, his in zip( selStrNE2, selStrHID ) :
	ZNStr 				= 'name ZN chain %s' % ( znHash[ NE.split()[1] ] )
	copy 				= znCore.copy()
	moveAtoms( copy.select( ZNStr ), to = np.array( [0,0,0] ) , ag = True )
	vectorAxes[ his ] 	= ( copy.select( NE ).getCoords(), znCore.select( NE ).getCoords() )
	

# Superposition vector matrix, heavily weighing the super positions of ring members, mostly NE2
weightArray = np.array( [
 [ 0.00],		# Atom N (index 0)
 [ 0.00],		# Atom CA (index 1)
 [ 0.00],		# Atom C (index 2)
 [ 0.00],		# Atom O (index 3)
 [ 0.00],		# Atom CB (index 4)
 [ 0.26],		# Atom CG (index 5)
 [ 0.08],		# Atom CD2 (index 6)
 [ 0.08],		# Atom ND1 (index 7)
 [ 0.08],		# Atom CE1 (index 8)
 [ 0.50],  		# Atom NE2 (index 9)
 [ 0.00],		# ATOM     10  HN
 [ 0.00],		# ATOM     13  HA  
 [ 0.00],		# ATOM     15  HB2 
 [ 0.00],		# ATOM     16  HB3 
 [ 0.00],		# ATOM     17  HD2 
 [ 0.00],		# ATOM     18  HD1 
 [ 0.00],		# ATOM     19  HE1
		])


# directory to store results
path2sites = sys.argv[4]
if not os.path.exists(  path2sites ):
	os.mkdir( path2sites )
		
# Make pdb for each combination of 8 rotamers fit to ring of the 4 ZN ligands, HIS (HID)
for rotSet in product( '01234567', repeat=4 ):
	ID 		= ''.join( list( rotSet ) )
	bisSite = znBare.copy()
	step 	= 1
	cnt 	= 1
	for rot, target in zip( rotSet, selStrHID ):
		
		hCopy 		= hisObjs[ int( rot ) ] 									## working copy of His rot
							
		## reset residue numbers, chain, & atom serials
		#print rot, target
		newResnum 	= target.split()[3]
		hCopy.setResnums( [ newResnum for x in hCopy ] )															
		newCh 		= target.split()[1]
		hCopy.setChids( [ newCh for x in hCopy ] )
		hCopy.setSerials( np.arange( cnt, len( hCopy ) + cnt ) )			

		#Super impose imidozols with weight-matrix driven super position
		hCopy 		= superpose( hCopy, znCore.select( target ), weights = weightArray )[0]
		
		# Add positioned histidine rotamer to metals and PTE substrate;
		bisSite 	+= hCopy

		step 		+= 1
		cnt 		+= len( hCopy )

	# write final prody atomGroup object pdb file
	moveAtoms( bisSite.select( 'resname OH name O1'), to = np.array( [0,0,0] ) , ag = True )
	print "\nwriting", os.path.join( path2sites, 'bisHis%s-0000.pdb.gz' % (ID))
	writePDB(  os.path.join( path2sites, 'bisHis%s-0000.pdb.gz' % (ID)), bisSite )
	bisSite.setTitle( ID ) 
	break


# For each of His rotamers made, go back and rotate each HIS atomGroup by 33 degrees
# Label each file by the rotation done & store path to each in .txt list file
# Label by index of rotation: for index n, degrees to rotate is 360*n/9. Label is 'nnnn', for each for 4 HIS
# Ignore potantial site if any HIS atom built is within 2.1 Angstroms of another HIS or Paraoxon ligand

print
pathList 			= []
originals2delete 	= []
for f in os.listdir( sys.argv[4] ):
	if '0000.pdb.gz' not in os.path.basename(f):
		continue
	path 		= os.path.join( path2sites, f )
	initPDB  	= parsePDB( path )

	# collect atomGroups (no easy command for this yet) for each his in array
	hisSet 		= []
	for resi in selStrHID:
		saveAtoms(  initPDB.select( resi ), 'tmp' )
		hisSet.append( loadAtoms( 'tmp.ag.npz' ) )


	# Write PDB for each combination of n*33 degree rotations (n = 1, 2, ...9) for each 4 HIS
	for rotIndex in product( '0123456789', repeat=4 ):


		# If initial state, check if
		if  not rotIndex == ('0', '0', '0', '0'): 
			rotIndex = ('3', '5', '9', '1')
		ID = ''.join( list(rotIndex) )

		# Begin with ligand + ZN-OH-ZN for each rotation
		rotPDB 	= znBare.copy()
		
		newPath = os.path.join( path2sites, f.split('-')[0] + '-%s.pdb.gz' % ( ID ) )  
		
		print ID
		print newPath


		# Decide to build or not based on clashes
		killFlg = 0

		# Perform each rotaion on each HIS object from superposition files from above
		for n, his, key, neStr, ring in zip( rotIndex, hisSet, selStrHID, selStrNE2, selStrRing):


			# No rotation condition, add un-rotated HIS if no clashes present 
			if n == '0':
				try:
					killFlg = len( rotPDB.select('not element Zn H not resname HOH and (within 2.3 of hisRot)', hisRot = his).getNames() ) 
					print 'died 1', n
					break
				except AttributeError:	
					try:
						killFlg = len( rotPDB.select('not element Zn not resname HOH and (within 1.4 of hisRot)', hisRot = his).getNames() ) 
						print 'died 2', n
						break
					except AttributeError:
						pass 
				rotPDB 	+= his 


				continue

			# I'm too dumb to figure out geometry... so instead line rotation axis up to 
			moveAtoms( his.select( neStr ), to = vectorAxes[ key ][0], ag = True )
			moveAtoms( rotPDB.select( 'name ZN chain %s' % ( znHash[neStr.split()[1]] )), to = np.array( [ 0, 0, 0 ] ) , ag = True)

			vector      = np.array( calcCenter( his.select( ring ) ) ) 
			vector      = vector / norm( vector )
			theta 		= int( n ) * 4 * pi / 10 

			# Calculate rotation matrix from vector of axis
			rotMatrix 	= rotation_matrix( vector , theta )

			# store & apply 4x4 matrix of rotation matrix + empty translation vector Prody object
			rotation 	= Transformation( rotMatrix, np.array( [0,0,0] )  )
			newHis		= applyTransformation( rotation, his.copy() )

			# move this rotated rotamer back to align with NE2 position of HOH-centered core
			moveAtoms( newHis.select( neStr ), to = vectorAxes[ key ][1], ag = True )
			moveAtoms( rotPDB.select( 'resname OH name O1' ), to = np.array( [ 0, 0, 0 ] ) , ag = True)

			#Check if new rotamer clashes any object in site, if not, kill building of this his-Zn site
			#clash = rotPDB.select('not element Zn H not resname HOH and (within 2.4 of hisRot)', hisRot = newHis)
			try:
				killFlg = len( rotPDB.select('not element Zn H not resname HOH and (within 2.4 of hisRot)', hisRot = newHis).getNames() ) 
				print 'died 3', n
				break
			except AttributeError:	
				try:
					killFlg = len( rotPDB.select('not element Zn not resname HOH and (within 1.4 of hisRot)', hisRot = newHis).getNames() ) 
					print 'died 4', n
					break
				except AttributeError:
					pass
			# Add rotamer if no clashes found
			rotPDB += newHis


		print "writing", newPath, '\n'
		if killFlg == 0:
			rotPDB.setTitle( ID ) 
			writePDB( newPath , rotPDB )
			pathList.append( path )
		else:
			if ID == '0000':
				originals2delete.append( newPath )
			break


	print


for f in originals2delete:
	os.remove( f )

pathListFile = open( 'hisSiteList.txt', 'w' )
for f in pathList:
	pathListFile.write( f + '\n' )

os.remove('tmp.ag.npz')

