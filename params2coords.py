

# example command line
# python ~/bin/params2coords.py params.txt ~/peptideAmyloid/parameterization/helix_prep/18_2ztaCEN.pdb ~/peptideAmyloid/OsakaModels/ ~/peptideAmyloid/parameterization/FullInterfaceCENTERED.pdb
# or
# python ~/bin/params2coords.py ~/peptideAmyloid/parameterization/params_OG.txt ~/peptideAmyloid/parameterization/helix_prep/18_Ideal_allALA.pdb ~/peptideAmyloid/OsakaModels/ ~/peptideAmyloid/parameterization/FullInterfaceCENTERED.pdb

# From a set of parameters for the placement of a helix relative to (0,0,0) on a flat amyloid surface 

# input file to read parameters from
#dBeta Theta N_t Z_n Z_c W_n; indexed by python array 0,1,2,3,4,5 ... and Phi added as 6th

from prody import *
import sys, os, numpy as np


##
#Shady  parameter break up for parallelization
# 349160 done so far
##

#### GLOBAL PARAMETERS
hLen 	= 25.71			# length of helix, based on 18 residue "idealized" alpha helix 
#hLen 	= 25.8			# length of helix, based on 18 residue fragment of GCN4 (7-24)
R_1		= 2.26		
# Degrees around the omega (W_n) axis between the first CA and the last CA (e.g. -120 from -pi to pi, or 240) 
rot 	= 102.8 		# Ideal alpha helix (poly ala GCN4 internal coords/angles/distances averaged & minimized)
#rot		= 119.4			# GCN4
offset	= np.array( [ 0.23, 9.86, 0.28 ] )

#### GLOBAL PARAMS END

########### Functions defined ########

def ang2rad( ang ):
	return ang*np.pi/180.0

def rnd( np_vector ):
	return [ round( x, 3 ) for x in np_vector ]


# Write a pdb formatted line
def pdb_line( coords, serial, resnum, chain):
		##  '{:<6}{:>5} {:<4}{:<1}{:<3} {:<1}{:>4}{:<1} {:>8}{:>8}{:>8}{:>6}{:>6} {:>2}{:>2}\n'.format( *lineList )
	lineList = [
		'ATOM',				#1	"ATOM " or "HETATM"	6	%-6s	01-06	[0:6]
		serial, 			#2	atom serial number	5	%5d	07-11	[6:11]
		'CA', 				#3	atom name	4	%4s	13-16	[12:16]
		'', 				#4	alternate location indicator	1	%1s	17	[16:17]
		'HIS', 				#5	residue name	3	%3s	18-20	[17:20]
		chain, 				#6	chain identifier	1	%1s	22	[21:22]
		resnum, 			#7	residue sequence number	4	%4d	23-26	[22:26]
		'', 				#8	code for insertion of residues	1	%1s	27	[26:27]
		str( coords[0] ),	#9	orthogonal coordinates for X (in Angstroms)	8	%8.3f	31-38	[30:38]
		str( coords[1] ),	#10	orthogonal coordinates for Y (in Angstroms)	8	%8.3f	39-46	[38:46]
		str( coords[2] ),	#11	orthogonal coordinates for Z (in Angstroms)	8	%8.3f	47-54	[46:54]
		'1.00', 			#12	occupancy	6	%6.3f	55-60	[54:60]
		'30.00', 			#13	temperature factor	6	%6.3f	61-66	[60:66]
		'C', 				#14	element symbol	2	%2s	77-78	[76:78]
		'' 					#15	charge on the atom	2	%2s	79-80	[78:80]
	]
	return '{:<6}{:>5} {:<4}{:<1}{:<3} {:<1}{:>4}{:<1}   {:>8}{:>8}{:>8}{:>6}{:>6}          {:>2}{:>2}\n'.format( *lineList )

class Helix:
	
	def __init__( self, dBeta, Theta, N_t, Z_n, Z_c, W_c, Phi ):

		# Coords of N and C terminal axis points
		self.X_n 	= N_t 	+ (.5*hLen) * (np.sin(Theta)) * (np.cos(Phi)) 
		self.Y_n 	= dBeta - (.5*hLen) * (np.cos(Theta)) * (np.cos(Phi)) 
		self.Z_n 	= Z_n 

		self.X_c 	= N_t 	- (.5*hLen) * (np.sin(Theta)) * (np.cos(Phi)) 
		self.Y_c 	= dBeta + (.5*hLen) * (np.cos(Theta)) * (np.cos(Phi)) 
		self.Z_c 	= Z_c 

		self.N_axis = np.array( [ self.X_n , self.Y_n , self.Z_n ]  )
		self.C_axis = np.array( [ self.X_c , self.Y_c , self.Z_c ]  )

		# tranformation vector (and reverse) placing c-terminus at origin
		self.trans 		= -1 * self.C_axis
		# normal vector through helical axix, when c-termini centered at origin
		self.axisVect 	= self.N_axis + self.trans
		self.axisNorm 	= self.axisVect / np.linalg.norm( self.axisVect )

		# position of reference point for C-alpha, given C-term axis point at origin
		# When rotation W_n = 0, then this point is C-alpha
		self.dX_c 	= R_1 * ( np.sin( Phi ) * np.sin( Theta ) )
		self.dY_c 	= R_1 * ( np.sin( Phi ) * np.cos( Theta ) * (-1) )
		self.dZ_c 	= R_1 *   np.cos( Phi ) * (-1)

		self.refV_b = np.array( [ self.dX_c, self.dY_c, self.dZ_c ] )
		# Take norm of vector, should be orthogonal to axis normal vector
		self.norm_b = self.refV_b / np.linalg.norm( self.refV_b )

		# Find final basis vector for new orthonormal frame (@ origin)
		self.norm_a = np.cross(  self.norm_b, self.axisNorm )

		# at center, find CA coords for a certain rotation round helix axis W_i, at helical radius R_1
		# given rotation angle W_n, define the C-terminal coordinates, given parametric circle arc eqs.
		## Details center = 0,0,0; orthogonal bases vectors b and a, normal to axis, rotation param = W_i
		self.CA_c_X = 0 + R_1 * self.norm_b[0] * np.cos( W_c ) + R_1 * self.norm_a[0] * np.sin( W_c )
		self.CA_c_Y = 0 + R_1 * self.norm_b[1] * np.cos( W_c ) + R_1 * self.norm_a[1] * np.sin( W_c )
		self.CA_c_Z = 0 + R_1 * self.norm_b[2] * np.cos( W_c ) + R_1 * self.norm_a[2] * np.sin( W_c )

		self.CA_c 	= np.array( [self.CA_c_X, self.CA_c_Y, self.CA_c_Z ] )

		# Same for N term coords (at center), given W_n = W_c + 120 degress
		W_n = W_c + ang2rad( rot )
		self.CA_n_X = 0 + R_1 * self.norm_b[0] * np.cos( W_n ) + R_1 * self.norm_a[0] * np.sin( W_n )
		self.CA_n_Y = 0 + R_1 * self.norm_b[1] * np.cos( W_n ) + R_1 * self.norm_a[1] * np.sin( W_n )
		self.CA_n_Z = 0 + R_1 * self.norm_b[2] * np.cos( W_n ) + R_1 * self.norm_a[2] * np.sin( W_n )

		self.CA_n 	= np.array( [self.CA_n_X, self.CA_n_Y, self.CA_n_Z ] )

		# Transform CA vectors back into onto termini of the helixal axis
		self.CA_n = self.CA_n + self.axisVect - self.trans
		self.CA_c = self.CA_c - self.trans

		self.coords_Axis = ( self.N_axis 	, self.C_axis 	)
		self.coords_CA	 = ( self.CA_n 		, self.CA_c 	)
		self.alignMarks  = np.array( [ self.CA_n, self.CA_c, self.C_axis, self.N_axis], dtype=float )


	### Function to write out some sanity check about the helical coordinates and omega (W_n) rotation
	def test_W( self ):
		print 'Axis Points: N, C', rnd( self.N_axis ) , rnd( self.C_axis )
		print 'axis vector & unit norm', rnd( self.axisVect ) , rnd( self.axisNorm )
		print 'ref vect b  & unit norm', rnd( self.refV_b ), rnd( self.norm_b )
		print 'dot & cross of these norms', round( np.dot(  self.norm_b, self.axisNorm ), 3 ), rnd( np.cross(  self.norm_b, self.axisNorm ) )
		print 'CA c-term coords', rnd( self.CA_c )
		print 'CA N-term coords', rnd( self.CA_n )
		print 'CA-axes distances', round( np.linalg.norm( self.CA_n - self.N_axis ) , 3 ), round( np.linalg.norm( self.CA_c - self.C_axis ) , 3 )
		return

	def __repr__(self):
		return str( self.coords_Axis )

########### Functions END ############


#################### MAIN ###########################################


# Make folder to store models
if not os.path.exists( sys.argv[3] ):
	os.mkdir( sys.argv[3] )


# for each line in param.txt input, calc array and make helix (calculate coordinates for end points)
params 		= []
inPDB 		= parsePDB( sys.argv[2] ) 
	# the order of selected atoms from PDB file matches those in helix object: Helix.alignMarks


index 		= 0
osaka 		= parsePDB( sys.argv[4] )
with open( sys.argv[1] ) as file:
	for i in file:
		if i[0] == '#':
#			if i[1].isdigit(): 
#				index = int( i[1:].split()[0] )
				continue
		print index,
		#sys.exit()
		#if len( i.split() ) != 6: continue
		
		params 	= [] 
		
		path 	= os.path.join( sys.argv[3], 'model_%d.pdb.gz' % index )
		alignMobile	= inPDB.select( 'ca resnum 7 24' ).copy()
		cnt = 0
		for x in i.split():
			if cnt in [ 1,5 ]:		# pick angle parameters and convert to radians
				params.append( ang2rad( float( x ) ) )
			else:
				params.append( float( x ) ) 
			cnt += 1

		phi = np.arcsin( ( params[3] - params[4] ) / hLen )
		params.append( phi )    								# Add Phi parameter, Z_n to Z_c angle
		print '\n'
		print 'model', index, 'params', params

		
		# derive coordinates for axis, given input parameters
		helix = Helix( *params )

		## unhash this line to print out sanity check for coordinates and axis points
		####helix.test_W() 

		# Align GCN4 helix, given N- and C-terminal CA's and helical Axes
		target 		= helix.alignMarks
		marks, transMat 	= superpose( alignMobile, target )
		print 'RMSD to alignment marks (n=4)', round( calcRMSD(target, alignMobile), 5 )

		model 	= applyTransformation( transMat, inPDB.select( 'chain A' ).copy() )
		model.setChids( [ 'Y' for x in model.iterAtoms() ] )
		model.setTitle( ' '.join(  i.split() ) )

		# Make the previous repeat and next repeat
		modelPrv= model.copy()
		modelNxt= model.copy()
		modelPrv.setChids( [ 'X' for x in model.iterAtoms() ] )
		modelNxt.setChids( [ 'Z' for x in model.iterAtoms() ] )
		modelPrv.setCoords( np.array( [ x - offset for x in model.getCoords() ] ) )
		modelNxt.setCoords( np.array( [ x + offset for x in model.getCoords() ] ) )

		model = modelPrv + model + modelNxt + osaka.copy() 


		##### Unhash this section to include the alignment marks in the PDB files
		#tar = AtomGroup('target Marks')
		#tar.setCoords( np.array( target ) )
		#tar.setResnums( [x for x in np.arange( 4 )] )
		#tar.setChids(['Z' for x in np.arange( 4 )]  )
		#tar.setElements( ['C' for x in np.arange( 4 )]  )
		#tar.setNames( ['CA' for x in np.arange( 4 )]  )
		#tar.setResnames( ['MSE' for x in np.arange( 4 )]  )
		#model += tar

		# write pdb for this parameter set
		writePDB( path, model )

		index += 1

#################### MAIN END #######################################