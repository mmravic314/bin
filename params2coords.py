## Marco Mravci Feb 2016 DeGrado Lab

# inputs required ( in order of arguments to script )
# input 1: a 5 character string of parameter values 
# input 2: 'ideal helix' all-ala backbone pdb to use as torch template. global variable defaults (length, pitch) for 18_Ideal_allALA.pdb
# input 3: path to save file to. usually ~/peptideAmyloid/OsakaModels/model_2319.pdb
# input 4: path to centered fibril PDb file

# outputs: model at path given in sys.argv[3], input 3

# example command line
## python ~/bin/params2coords.py '0 90 0 5 5 180' ~/peptideAmyloid/parameterization/helix_prep/18_Ideal_allALA.pdb ~/peptideAmyloid/OsakaModels/model_0.pdb ~/peptideAmyloid/parameterization/FullInterfaceCENTERED.pdb

# From a set of parameters for the placement of a helix relative to (0,0,0) on a flat amyloid surface 
#dBeta Theta N_t Z_n Z_c W_n; indexed by python array 0,1,2,3,4,5 ... and Phi added as 6th

from prody import *
import sys, os, numpy as np


#### GLOBAL PARAMETERS
hLen 	= 25.71			# length of helix, based on 18 residue "idealized" alpha helix 
#hLen 	= 25.8			# length of helix, based on 18 residue fragment of GCN4 (7-24)
R_1		= 2.26		
# Degrees around the omega (W_n) axis between the first CA and the last CA (e.g. -120 from -pi to pi, or 240) 
rot 	= 102.8 		# Ideal alpha helix (poly ala GCN4 internal coords/angles/distances averaged & minimized)
#rot		= 119.4			# GCN4
## Vector to translate each torch helix for wall-paper symmetry, given average vector linking "identical" atoms between fibril repeat unit (2 strands)
offset	= np.array( [ 0.23, 9.86, 0.28 ] )

#### GLOBAL PARAMS END

########### Functions defined ########

def ang2rad( ang ):
	return ang*np.pi/180.0

def rnd( np_vector ):
	return [ round( x, 3 ) for x in np_vector ]



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

# for each line in param.txt input, calc array and make helix (calculate coordinates for end points)
params 		= []
inPDB 		= parsePDB( sys.argv[2] ) 
# the order of selected atoms from PDB file matches those in helix object: Helix.alignMarks



osaka 		= parsePDB( sys.argv[4] )
i 			= sys.argv[1]					# parameter string, given at input 


if len( i.split() ) != 6:
				print 'Error in params string input, incorrect number of parameters, 6 needed: dBeta Theta N_t Z_n Z_c W_n'
				sys.exit()
else:

		params 	= [] 
		
#		path 	= os.path.join( sys.argv[3], 'model_%d.pdb.gz' % index )
		path 	= sys.argv[3]	# debugging option

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
	#	print '\n'
		print 'model', os.path.basename( sys.argv[3] ), 'params', params

		
		# derive coordinates for axis, given input parameters
		helix = Helix( *params )

		## unhash this line to print out sanity check for coordinates and axis points
		####helix.test_W() 

		# Align ideal helix, given N- and C-terminal CA's and helical Axes
		target 		= helix.alignMarks
		marks, transMat 	= superpose( alignMobile, target )
		#print 'RMSD to alignment marks (n=4)', round( calcRMSD(target, alignMobile), 5 )

		model 	= applyTransformation( transMat, inPDB.select( 'chain A' ).copy() )
		model.setChids( [ 'Y' for x in model.iterAtoms() ] )

		# Make the previous repeat and next repeat
		modelPrv= model.copy()
		modelNxt= model.copy()
		modelPrv.setChids( [ 'X' for x in model.iterAtoms() ] )
		modelNxt.setChids( [ 'Z' for x in model.iterAtoms() ] )
		modelPrv.setCoords( np.array( [ x - offset for x in model.getCoords() ] ) )
		modelNxt.setCoords( np.array( [ x + offset for x in model.getCoords() ] ) )

		model = modelPrv + model + modelNxt + osaka.copy() 
		model.setTitle( ' '.join(  i.split() ) )

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


#################### MAIN END #######################################