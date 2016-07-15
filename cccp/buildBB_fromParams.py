## call Gevorg's matlab code from octave, with input parameters from a text file
# eventually... get C-alpha coords and build full backbone coordiantes

import sys, os, subprocess as sp, numpy as np
from collections import defaultdict

########### Global definitions for internal coordinates ##############

Nt_dist = 1.45
Nt_ang	= 95*2*np.pi/360
Nt_dih	= -65*2*np.pi/360

Ct_dist = 1.52
Ct_ang	= 21*2*np.pi/360
Ct_dih	= 79*2*np.pi/360


########## Global Defs END ###########################################



#########  functions to calculate remaining backbone atoms ##############

# input: 3 consecuative 3 CA coords (ca_0,ca_1+,ca_2+) relative to point N (resi 0), an distance (ca_0,N), angle (ca_1+,ca_0,N), dihedral (ca_2+,ca_1+,ca_0,N), then initial coordinate guess
# return: coordinates to that free point N, in R3 when N is a terminal nitrogen
def calcCoordsNt( ca_arr, l_b3, ang, dihed, init , nt=True):

	# Set up known info/vectors for points i,j,k... solving for point l (x_l, y_l, z_l)
	i 		= ca_arr[2]
	j 		= ca_arr[1]
	k 		= ca_arr[0]

	b_1 		= j - i
	b_2 		= k - j
	l_b1	= np.linalg.norm( b_1 )
	l_b2	= np.linalg.norm( b_2 )

	n1 		= np.cross( b_1, b_2 )

	# Set initial guess

	## prepare vector valued functions, functions of (x_l, y_l, z_l)
	## And below them, list their partial derivatives for each: x_l, y_l, z_l

	# distance restraint aka kl distance
	def f1( x_l, y_l, z_l ):
		return (x_l - k[0])**2 + (y_l - k[1])**2 + (z_l - k[2])**2 - l_b3**2

	def df1_dx( x_l, y_l, z_l  ):
		return 2*( x_l - k[0] )
	def df1_dy( x_l, y_l, z_l  ):
		return 2*( y_l - k[1] )
	def df1_dz( x_l, y_l, z_l ):
		return 2*( z_l - k[2] )

	# Angle restraint aka j,k,l angle... b2, b3 vector/dot product
	def f2( x_l, y_l, z_l ):
		return b_2[0]*( x_l - k[0] ) + b_2[1]*( y_l - k[1] ) + b_2[2]*( z_l - k[2] ) + np.cos( ang ) * l_b2 * l_b3

	def df2_dx(x_l, y_l, z_l  ):
		return b_2[0]
	def df2_dy( x_l, y_l, z_l ):
		return b_2[1]
	def df2_dz( x_l, y_l, z_l ):
		return b_2[2]


	# dihedral restraint aka angle between planes b1 x b2 and b2 x b3, via dot product of plane's normal vectors
	#   coefficients for b3_vector, used in function 3 and it's partial derivatives...
	a = n1[1] * b_2[2] - n1[2] * b_2[1]
	b = n1[2] * b_2[0] - n1[0] * b_2[2]
	c = n1[0] * b_2[1] - n1[1] * b_2[0]
	#if nt:
	#	sign = -1
	#else:
	#	sign = 1

	def f3( x_l, y_l, z_l ):
		return a * ( x_l - k[0] ) + b * ( y_l - k[1] ) + c * ( z_l - k[2] )  +np.cos( dihed ) * l_b1 * l_b3 * l_b2**2

	def df3_dx( x_l, y_l, z_l  ):
		return a
	def df3_dy( x_l, y_l, z_l  ):
		return b
	def df3_dz( x_l, y_l, z_l ):
		return c

	##### Newton's method solver section
	
	err 	= 0.001		# acceptable error
	# first round
	u  		= init
	def F( v ):
		return np.array( [ f1(*v), f2(*v), f3(*v) ] )		# vector of functions
	def J( v ):
		return np.array( [ 
						[ df1_dx( *v ), df1_dy( *v ), df1_dz( *v ) ], 
						[ df2_dx( *v ), df2_dy( *v ), df2_dz( *v ) ], 
						[ df3_dx( *v ), df3_dy( *v ), df3_dz( *v ) ] ] )
	#print F(u), '\n'
	#print J(u)
	e 		= ( np.linalg.norm( F(u) ) )**0.5
	tries 	= 20
	while e > err and tries > 0:
		del_u = np.linalg.solve( J(u),-1*F(u) )
	#print del_u
		u += del_u
		e = ( np.linalg.norm( F(u) ) )**0.5
#		print e, err
		tries -= 1
	if tries == 0:
		print 'SOLVER DID NOT CONVERGE, DYING'
		sys.exit()
	else:
		return u, [ '%.3f'% n for n in u ] 		# return strong formatted array

 

	# initial guess should be pretty close



	return


# input: 3 consecuative 3 CA coords (ca_0,ca_1+,ca_2+) relative to point N (resi 0), an distance (ca_0,N), angle (ca_1+,ca_0,N), dihedral (ca_2+,ca_1+,ca_0,N), then initial coordinate guess
# return: coordinates to that free point N, in R3 when N is a terminal nitrogen
def calcCoordsCore( ca_arr, l_b3, ang, dihed, init ):

	# Set up known info/vectors for points i,j,k... solving for point l (x_l, y_l, z_l)
	i 		= ca_arr[2]
	j 		= ca_arr[1]
	k 		= ca_arr[0]

	b_1 		= j - i
	b_2 		= k - j
	l_b1	= np.linalg.norm( b_1 )
	l_b2	= np.linalg.norm( b_2 )

	n1 		= np.cross( b_1, b_2 )

	# Set initial guess

	## prepare vector valued functions, functions of (x_l, y_l, z_l)
	## And below them, list their partial derivatives for each: x_l, y_l, z_l

	# distance restraint aka kl distance
	def f1( x_l, y_l, z_l ):
		return (x_l - k[0])**2 + (y_l - k[1])**2 + (z_l - k[2])**2 - l_b3**2

	def df1_dx( x_l, y_l, z_l  ):
		return 2*( x_l - k[0] )
	def df1_dy( x_l, y_l, z_l  ):
		return 2*( y_l - k[1] )
	def df1_dz( x_l, y_l, z_l ):
		return 2*( z_l - k[2] )

	# Angle restraint aka j,k,l angle... b2, b3 vector/dot product
	def f2( x_l, y_l, z_l ):
		return b_2[0]*( x_l - k[0] ) + b_2[1]*( y_l - k[1] ) + b_2[2]*( z_l - k[2] ) + np.cos( ang ) * l_b2 * l_b3

	def df2_dx(x_l, y_l, z_l  ):
		return b_2[0]
	def df2_dy( x_l, y_l, z_l ):
		return b_2[1]
	def df2_dz( x_l, y_l, z_l ):
		return b_2[2]


	# dihedral restraint aka angle between planes b1 x b2 and b2 x b3, via dot product of plane's normal vectors
	#   coefficients for b3_vector, used in function 3 and it's partial derivatives...
	a = n1[1] * b_2[2] - n1[2] * b_2[1]
	b = n1[2] * b_2[0] - n1[0] * b_2[2]
	c = n1[0] * b_2[1] - n1[1] * b_2[0]
	#if nt:
	#	sign = -1
	#else:
	#	sign = 1

	def f3( x_l, y_l, z_l ):
		return a * ( x_l - k[0] ) + b * ( y_l - k[1] ) + c * ( z_l - k[2] )  +np.cos( dihed ) * l_b1 * l_b3 * l_b2**2

	def df3_dx( x_l, y_l, z_l  ):
		return a
	def df3_dy( x_l, y_l, z_l  ):
		return b
	def df3_dz( x_l, y_l, z_l ):
		return c

	##### Newton's method solver section
	
	err 	= 0.001		# acceptable error
	# first round
	u  		= init
	def F( v ):
		return np.array( [ f1(*v), f2(*v), f3(*v) ] )		# vector of functions
	def J( v ):
		return np.array( [ 
						[ df1_dx( *v ), df1_dy( *v ), df1_dz( *v ) ], 
						[ df2_dx( *v ), df2_dy( *v ), df2_dz( *v ) ], 
						[ df3_dx( *v ), df3_dy( *v ), df3_dz( *v ) ] ] )
	#print F(u), '\n'
	#print J(u)
	e 		= ( np.linalg.norm( F(u) ) )**0.5
	tries 	= 20
	while e > err and tries > 0:
		del_u = np.linalg.solve( J(u),-1*F(u) )
	#print del_u
		u += del_u
		e = ( np.linalg.norm( F(u) ) )**0.5
#		print e, err
		tries -= 1
	if tries == 0:
		print 'SOLVER DID NOT CONVERGE, DYING'
		sys.exit()
	else:
		return u, [ '%.3f'% n for n in u ] 		# return strong formatted array



########### END Helper functions ##################


########## MAIN ###################################

path2Funct 	= sys.argv[1]
paramsF 	= sys.argv[2]	

params 		= []

with open( paramsF ) as fin:
	for p in fin:
		if p[0] == '#': continue

		# parse parameters, format for octave function call
		par 		= p.rstrip().split(',')
		chN, chL 	= tuple( [ int( n ) for n in par[:2] ] )
		call 		= path2Funct + '(%s)' % p.rstrip()

		# call octave with given parameters, record output. Quit if function error found
		cmd 	= ['octave', '--eval', call]
		out 	= sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
		stdout, err 	= out.communicate()
		if len(err) > 0:
			print 'SOME ERROR FOUND... quitting\n', err
			sys.exit()


		# parse coordinates, save into chains
		coords 	= defaultdict(list)
		atm		= 0 
		ch 		= 65 		# Ascii value for chain to start as A
		for k in stdout.split('\n'):
			if len( k.rstrip() ) < 8: continue

			x, y, z = tuple( [ round( float( n ), 3) for n in k.rstrip().split() ] )
			#print atm, x,y,z, chr(ch)
			coords[ chr(ch) ].append( np.array( [x,y,z] ) )
			
			# make new chain if atom named is multiple of chain length
			atm += 1
			if atm%chL == 0:
				ch +=1

		# With coordinates, write full backbone atoms for chain

		for c in sorted( coords.keys() ):
			stp = 0
			print 'writing backbone coordinates for Chain', c
			
			# calculate N-terminal N and C, then O
			
			Ninit 	= coords[c][0]  + [0.0, 0.0, -1.45]
			Cinit 	= coords[c][0]  + [0.0, 0.0, 1.52]
			Nt 		= calcCoordsNt( coords[c][:3], Nt_dist, Nt_ang, Nt_dih, Ninit )
			Ct 		= calcCoordsNt( coords[c][:3], Ct_dist, Ct_ang, Ct_dih, Cinit )
			print '\n', coords[c][:3]
			print Nt
			print Ct
			print np.linalg.norm( Nt[0] - Ct[0] ), '\n'

			Ninit 	= coords[c][1]  + [0.0, 0.0, -1.45]
			Cinit 	= coords[c][1]  + [0.0, 0.0, 1.52]
			Nt 		= calcCoordsNt( coords[c][1:4], Nt_dist, Nt_ang, Nt_dih, Ninit )
			Ct 		= calcCoordsNt( coords[c][1:4], Ct_dist, Ct_ang, Ct_dih, Cinit )
			print '\n', coords[c][:3]
			print Nt
			print Ct
			print np.linalg.norm( Nt[0] - Ct[0] ), '\n'
			# calculate all middle N, C, O
			# calculate C-terminal N, C, O


			sys.exit()



######## EXTRA CRAP/NOTES BELOW
'''
	stp = 0
	for k in coordAr:
		#txt = 'ATOM    281  N   HIS A  39     %.2f %.2f  %.2f  1.00 20.00           N  ' % ( k[0], k[3] , k[2]  )
		lineInfo = ('ATOM', stp + 1, 'CA',' ', 'ALA', chain, stp + 1, ' ', k[0], k[1], k[2], 1.00, 20.00, 'C', '  ')
		txt = "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s" % ( lineInfo )
'''