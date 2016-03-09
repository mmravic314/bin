# print list file of params for the grid search of single helices patterned on amyloid


import os, sys, numpy as np

## Global defs
hLen 	= 25.8	    # length of helix, based on 18 residue "idealized" alpha helix aka GCN4 
tLen	= 33		# longest distance in y direction	
##

## helper functions

def ang2rad( ang ):
	return ang*np.pi/180.0

def rnd( np_vector ):
	return [ round( x, 3 ) for x in np_vector ]

## 

# print to file as follows
# dBeta Theta N_t Z_n Z_c W_n
# But eval N_t last, since its range can depend on theta and phi

# left to right of these params, move through grid
# once previous value set, calculate elible range of remaining values for next parameter if limited by previous value
# example 

# Here the completely independents are spaced into vectors... modify to make grid search sparser or denser

dB 		= np.arange( -5,   6, 1  )
W_n		= np.arange(  0, 360, 20 ) 
Theta 	= np.arange( 50, 140, 10 )
Z_n		= np.arange(  3,   9,  1 )
Z_c		= np.arange(  3,  11,  1 )
# must calc N_t range from previous ranges

params  = []		# array to fill out, with arrays of these parameters

## Is a 5-layer nested loop the best way to do this?

## split it up... also skip some number that are already made (this is a tmp option), comment out lines with comment with # tmp to dump this

index = 0  ## tmp
for b in dB:
  for w in W_n:
    for t in Theta:
      for zn in Z_n:
        for zc in Z_c:
   		   # calc dependent param Phi from Z values of axis points
          phi = np.arcsin( ( zn - zc ) / hLen )	

          inc = 1
   		# calculate range of Nt given remaining params... make sure to convert Theta to radians for this
          N_t = np.arange( -0.5*(tLen - hLen * np.sin( ang2rad( t ) ) * np.cos(phi) ), inc + 0.5 * (tLen - hLen * np.sin(  ang2rad(t) ) * np.cos(phi) ), inc)
          for n in N_t:
            pSet = ' '.join( [str(x) for x in [ b, t, n, zn, zc, w ] ] )
            

            if index < 349160: 
                   index += 1
                   continue    #tmp


            params.append( pSet )
          


            if index == 500000:
                  txt1 = '#349160 dBeta Theta N_t Z_n Z_c W_n\n' + ' '.join( [ x + '\n' for x in params ] )
                  params =[]
                  outFile = open( sys.argv[1][:-4] + '1.txt'  , 'w' )
                  outFile.write( txt1 )

            if index == 650000:
                  txt2 = '#500000 dBeta Theta N_t Z_n Z_c W_n\n' + ' '.join( [ x + '\n' for x in params ] )
                  params =[]
                  outFile = open( sys.argv[1][:-4] + '2.txt'  , 'w' )
                  outFile.write( txt2 )

            if index == 800000:
                  txt3 = '#750000 dBeta Theta N_t Z_n Z_c W_n\n' + ' '.join( [ x + '\n' for x in params ] )
                  params =[]
                  outFile = open( sys.argv[1][:-4] + '3.txt'  , 'w' )
                  outFile.write( txt3 )

            if index == 970793:
                  txt4 = '#750000 dBeta Theta N_t Z_n Z_c W_n\n' + ' '.join( [ x + '\n' for x in params ] )
                  params =[]
                  outFile = open( sys.argv[1][:-4] + '4.txt'  , 'w' )
                  outFile.write( txt4 )



            index += 1 #

print index

print 'total param sets:', len( params )
txt = '# dBeta Theta N_t Z_n Z_c W_n\n' + ' '.join( [ x + '\n' for x in params ] )
#outFile = open( sys.argv[1], 'w' )
#outFile.write( txt )






 

