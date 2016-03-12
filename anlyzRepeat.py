# print some stats about vector between repeating atoms in input file, given pairs of chains with unit cell for 2D wal paper symmetry


from prody import *
import sys, numpy as np, math

osk = parsePDB( sys.argv[1], subset ='bb' )

s1, s2 = osk.select( 'chain G' ).copy(), osk.select( 'chain H' )

x, y, z = [], [], []
for a1,a2 in zip( s1, s2 ):
	c = a2.getCoords() - a1.getCoords()

	x.append( c[0] )
	y.append( c[1] )
	z.append( c[2] )

print np.mean( x ), np.std( x )
print np.mean( y ), np.std( y )
print np.mean( z ), np.std( z )