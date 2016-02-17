#params2coords.py
# From a set of parameters for the placement of a helix relative to (0,0,0) on a flat amyloid surface 

# input file to read parameters from
#dBeta Theta N_t Z_n Z_c W_n; indexed by python array 0,1,2,3,4,5 ... and Phi added as 6th

from prody import *
import sys, os, numpy as np


#### GLOBAL PARAMETERS
hLen 	= 25.8			# length of helix, based on 18 residue "idealized" alpha helix 
R_1		= 2.26		
rot		= -120			# Degrees around the omega (W_n) axis between the first CA and the last CA (-120, from -pi to pi, or 240) 

#### GLOBAL PARAMS END

########### Functions defined ########

def ang2rad( ang ):
	return ang*np.pi/180.0

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
	
	def __init__( self, dBeta, Theta, N_t, Z_n, Z_c, W_n, Phi ):

		# Coords of N and C terminal axis points
		self.X_n 	= round( N_t 	+ (.5*hLen) * (np.sin(Theta)) * (np.cos(Phi)) , 3 )
		self.Y_n 	= round( dBeta 	+ (.5*hLen) * (np.cos(Theta)) * (np.cos(Phi)) , 3 )
		self.Z_n 	= round( Z_n , 3)

		self.X_c 	= round( N_t 	- (.5*hLen) * (np.sin(Theta)) * (np.cos(Phi)) , 3 )
		self.Y_c 	= round( dBeta 	- (.5*hLen) * (np.cos(Theta)) * (np.cos(Phi)) , 3 )
		self.Z_c 	= round( Z_c , 3 )

		self.N_axis = np.array( [ self.X_n , self.Y_n , self.Z_n ]  )
		self.C_axis = np.array( [ self.X_c , self.Y_c , self.Z_c ]  )

		# mods to each Axis end point for approx to CA atom coords
		self.dX_n 	= round( R_1 * np.sin( Theta ) * np.sin( Phi ) * np.cos( W_n ) , 3 )
		self.dY_n 	= round( R_1 * np.cos( Theta ) * np.sin( Phi ) * np.sin( W_n ) , 3 )
		self.dZ_n 	= round( R_1 * np.cos( Phi ) * np.cos( W_n )  , 3 )

		self.dX_c 	= round( R_1 * np.sin( Theta ) * np.sin( Phi ) * np.cos( W_n + ang2rad( rot )) , 3 )
		self.dY_c 	= round( R_1 * np.cos( Theta ) * np.sin( Phi ) * np.sin( W_n + ang2rad( rot )) , 3 )
		self.dZ_c 	= round( R_1 * np.cos( Phi ) * np.cos( W_n + ang2rad( rot ) )  , 3 )

		self.dW		= ( [ self.dX_n , self.dY_n , self.dZ_n ] , [ self.dX_c , self.dY_c , self.dZ_c ]  )

		self.N_w = self.N_axis + np.array( self.dW[0] )
		self.C_w = self.C_axis + np.array( self.dW[1] )

		self.coords_Axis = ( self.N_axis , self.C_axis )
		self.coords_CA	 = ( self.N_w 	 , self.C_w    )

	def test_W( self ):
		print self.dW
		print 'N dist', np.linalg.norm( self.N_w - self.N_axis )
		print 'C dist', np.linalg.norm( self.C_w - self.C_axis )
		return	

	def __repr__(self):
		return str( self.coords_CA )






########### Functions END ############


#################### MAIN ###########################################


# Make folder to store models
if not os.path.exists( sys.argv[3] ):
	os.mkdir( sys.argv[3] )

# for each line in param.txt input, calc array and make helix (calculate coordinates for end points)
params 	= []
inPDB 	= parsePDB( sys.argv[2] ) 
osaka 	= parsePDB( sys.argv[4] )
with open( sys.argv[1] ) as file:
	for i in file:
		if i[0] == '#' or len( i.split() ) != 6: continue
		params 	= [] 
		index 	= 0
		path 	= os.path.join( sys.argv[3], 'model_%d.pdb' % index )

		for x in i.split():
			if index in [ 1,5 ]:		# pick angle parameters and convert to radians
				params.append( ang2rad( float( x ) ) )
			else:
				params.append( float( x ) ) 
			index += 1

		params.append(  np.arcsin( ( params[4] - params[3] ) / hLen ) )    # Add Phi parameter
		
		print params
		helix = Helix( *params )
		print 

		print "CA-distance", np.linalg.norm(  helix.coords_CA[0] - helix.coords_CA[1] )
		print "\nCA coords (N,C)", helix
		print "\nAxis (N,C)", helix.coords_Axis
		print

		helix.test_W()
		# 
		sys.exit()
		

		## Print out test PDB file format
		print
		print pdb_line( helix.coords_CA[0], '1', '1', 'P')
		print pdb_line( helix.coords_CA[1], '2', '2', 'P')
		print pdb_line( helix.coords_Axis[0], '3', '1', 'X')
		print pdb_line( helix.coords_Axis[1], '4', '2', 'X')



		target = np.array( [ helix.coords_CA[0], helix.coords_CA[1] ] )
		mobile = inPDB.copy().select(' ca resnum 7 24')			# hard coded indices, but can rewrite this for variable peptide termini
		# align GCN4 helix to termini markers
		print calcRMSD(target,mobile)
		transformMat = superpose(mobile, target )[1] 
		print calcRMSD(target,mobile)	
		mobile = applyTransformation( transformMat, inPDB.copy() ) + osaka.copy()
		mobile.setTitle( ' '.join( [str(x) for x in params] ) )
		#writePDB( path, mobile )



#################### MAIN END #######################################