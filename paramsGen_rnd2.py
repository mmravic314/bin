# print list file of params for the grid search of single helices patterned on amyloid


# input 1: directory to print output params#_rnd2.txt


import os, sys, numpy as np, random

## Global defs
hLen  = 30.216     # length of helix, based on 21 residue "idealized" alpha helix   
tLen	= 35		     # longest distance in x direction	
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

dB 		= np.arange( -4.93, 5, 4.93/5  )
W_n		= np.arange(  15, 360, 20 ) 
Theta = np.arange( -102, -84, 2 ) 
Z_n		= np.arange(  4,  8,  0.5 )
Z_c		= np.arange(  3,  8,  0.5 )
# must calc N_t range from previous ranges

params  = []		# array to fill out, with arrays of these parameters



mN = 10
mxN= 10

for b in dB:
  for w in W_n:
    for t in Theta:
      for zn in Z_n:
        for zc in Z_c:
   		   # calc dependent param Phi from Z values of axis points
          phi = np.arcsin( ( zn - zc ) / hLen )	

          t_x = ang2rad( np.abs( t ) )

          ext = 0.5 *  ( tLen - hLen * np.sin( t_x ) * np.cos(phi) ) 

          inc = 3.5/5  # max N_t value for accessed PHI and THETA range

   		# calculate range of Nt given remaining params... make sure to convert Theta to radians for this
          N_t =  np.array( [ x for x in np.arange( -1 * 3.5,  3.5 + inc, inc )  if np.abs( x ) <= ext ] )

          for n in N_t:

            pSet = ' '.join( [str(  round( x, 2 ) ) for x in [ b, t, n, zn, zc, w ] ] )
            
            params.append( pSet )

i = len( params )
print 'total param sets:', len( params )



index = 0
step  = 20013
txt   = '%d# dBeta Theta N_t Z_n Z_c W_n\n' % (index/step)     # + ' '.join( [ x + '\n' for x in params ] 
#
# SYNTAX ####
# Goofy note: must end each xtxt string with '##' since it adds newline char to final params text line.
#             This allows the final line to be accessed by bash while loops and 'read' command
#             This is handled in main bash script that calls param file outputs of this script 


rng = random.sample( xrange( i ), i )
for k in rng:


  txt += params[k] + '\n'
  index +=1

  if index % step  == 0:

      outFile = open(  os.path.join( os.path.dirname( sys.argv[1] ), 'params%d_rnd2.txt' % ( index/step ) )    , 'w' )
      outFile.write( txt + '##')      
      outFile.close()
      txt = '%d# dBeta Theta N_t Z_n Z_c W_n\n' % (index/step)
 
# final print out of txt to file:
outFile = open(  os.path.join( os.path.dirname( sys.argv[1] ), 'params%d_rnd2.txt' % ( int( str( index/step ).split('.')[0] ) +1 ) )    , 'w' )
outFile.write( txt + '##')
outFile.close()
print 'wrote %d files! last file had %d lines' % ( index, len( txt.split('\n') ) )

#txt = '# dBeta Theta N_t Z_n Z_c W_n\n' + ' '.join( [ x + '\n' for x in params ] )
#outFile = open( sys.argv[1], 'w' )
#outFile.write( txt )






 

