# print list file of params for the grid search of single helices patterned on amyloid


import os, sys, numpy as np, random

## Global defs
hLen 	= 25.8	    # length of helix, based on 18 residue "idealized" alpha helix aka GCN4 
tLen	= 33		# longest distance in y direction	
##

## helper functions

def ang2rad( ang ):
	return ang*np.pi/180.0

def rnd( np_vector ):
	return [ round( x, 3 ) for x in np_vector ]

def shuffle( lst ):
  return



## 

# print to file as follows
# dBeta Theta N_t Z_n Z_c W_n
# But eval N_t last, since its range can depend on theta and phi

# left to right of these params, move through grid
# once previous value set, calculate elible range of remaining values for next parameter if limited by previous value
# example 

# Here the completely independents are spaced into vectors... modify to make grid search sparser or denser

dB 		= np.arange( -4.93, 5, 4.93/5  )
W_n		= np.arange(  0, 360, 30 ) 
Theta = np.append( np.arange( 50, 140, 15 ), np.arange( -130, -50, 10 ) )
Z_n		= np.arange(  3,  7,  1 )
Z_c		= np.arange(  3,  9,  1 )
# must calc N_t range from previous ranges

params  = []		# array to fill out, with arrays of these parameters



## split it up... by var 'step'
step = 50000

for b in dB:
  for w in W_n:
    for t in Theta:

      if t > 0: continue


      for zn in Z_n:
        for zc in Z_c:
   		   # calc dependent param Phi from Z values of axis points
          phi = np.arcsin( ( zn - zc ) / hLen )	
          if t > 180:
            t_x = ang2rad( t-180 )
          else:
            t_x = ang2rad( t )
          gap = tLen - hLen

          inc = 1
   		# calculate range of Nt given remaining params... make sure to convert Theta to radians for this
          N_t = np.arange( -0.5*( gap * np.sin( t_x ) * np.cos(phi) ), 1 + 0.5 * (gap * np.sin( t_x ) * np.cos(phi) ), inc)
          print N_t, t

          sys.exit()


          for n in N_t:

            pSet = ' '.join( [str(  round( x, 2 ) ) for x in [ b, t, n, zn, zc, w ] ] )
            
            params.append( pSet )

i = len( params )
print 'total param sets:', len( params )

sys.exit()
index = 0
step  = 50000
txt   = '#%d dBeta Theta N_t Z_n Z_c W_n\n' % (index/step)     # + ' '.join( [ x + '\n' for x in params ] 

rng = random.sample( xrange( i ), i )
for k in rng:


  txt += params[k] + '\n'
  index +=1

  if index % step  == 0:

      outFile = open(  os.path.join( os.path.dirname( sys.argv[1] ), 'params%d.txt' % ( index/step ) )    , 'w' )
      outFile.write( txt )
      outFile.close()
      txt = '%d# dBeta Theta N_t Z_n Z_c W_n\n' % (index/step)
 



#txt = '# dBeta Theta N_t Z_n Z_c W_n\n' + ' '.join( [ x + '\n' for x in params ] )
#outFile = open( sys.argv[1], 'w' )
#outFile.write( txt )






 

