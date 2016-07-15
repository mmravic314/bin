# Marco Mravic DeGrado lab July 2016
# Generate a parameters file 
# python ccc3_ParamsGen.py ~/tertBuilding/MS1/params.txt

import sys, os, numpy as np

outFile = sys.argv[1]

# (chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin)


### Variable Params

R_o 	= np.arange( 5.5, 7.625, 0.125 )
W_o 	= np.arange( -0.125, -3.125, -0.25 )
Phi_1	= np.arange( 75, 135, 5 )

###

### Fixed params

chains		= '3'
chL 		= '32'
#r_o		VARIABLE
r_1			= '2.26'
phi_o		= '[120, 240]'
topolopy 	= '[1 1]'
z_aa		= '[0.0 0.0]'

#N_alpha		= 3.6
#fix_w 		= 100.7

###

## Remaining parameters are calculated given the variable parameter

pSet = []
# (chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin)
txt = '# (chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin)\n'
for r_o in R_o:
	for phi1 in Phi_1:
		for w_o in W_o:
		

			## Caclulated parameters
			w_1 	= 100 - w_o
#			n_minor	= 360/w_1

			N_super = 360/w_o # can also calculate by inverse of major/minor twist differential: 1/3.6 - (w_1/360)

			Pitch 	= -1 * ( (N_super * 1.51)**2 - 4 * ( np.pi**2 ) * ( r_o**2 ) )**0.5
			alpha_Pang	= round( 180* np.arctan( 2*np.pi*r_o / Pitch) / np.pi, 3)

			vals = [ chains, chL, str( r_o ), str( r_1 ), str( w_o ), str( w_1 ), str(alpha_Pang), 
			'[%.1f, %.1f, %.1f]' % ( phi1, phi1, phi1 ), topolopy, phi_o, z_aa, '\'zoffaa\''

			]
			txt += ', '.join( vals ) + '\n'


		#	print 360*( (1/ 3.6) - ( 1/ 3.5) )
		
#		sys.exit()
outF = open( outFile, 'w' )
outF.write( txt )



