# analyze geometry of CLuster 6 (SOL) helica pairs to DF1-L13A resi 27-48 [helix 2] 
#(that dont clash to helix 1 resi 1-24)

# Define center of mass to helix 1 and 2. Then find point on fit 'incident' helix  that crosses 
# center of mass in z direction. Find angle in x-y plane between three helical axes. 
# find/plot crossing angle & dZ to center of mass

from prody import * 
import sys, os, numpy as np, cPickle as pic, subprocess as sp
from numpy.linalg import norm
from numpy import dot


#### Helper functions

# Given coordinates of a helical axis, find projection on to XY plane using closest point on axis to z=0
# If helix comes 1.5 A < d < 3.2 A of XY plane, linearly interpolate using last three points.
# Return False if helical axis stops >= 3.2 A from XY plane and helica curvature < 100 degrees
def XYpoint( coordList, radius ):
	if len( coordList ) < 1:
		return 0

	closest = sorted( coordList, key = lambda x: np.fabs( x[-1] ) )[0]
	if closest[-1] <= 1.5:
		return ( closest[0], closest[1] , 0 )

	elif  closest[-1] < 3.2 and radius < 90:
		print 'XXXXXXXXXXXXXXXXX'
		return 0

	# If 'straight' helix is far away from XY-plane, take X-Y plane from line intersetion 
	elif radius >= 90:
		return ( closest[0], closest[1] , 0 )
	
	# Return FALSE is too far and too curved to easily determine XY plane interaction of helix axis
	else:
		return 0

# calc angle in degrees given three points via cosine rule
def angle( mobile, vertex, origin ):

	a = sum((mobile - vertex)**2)**0.5
	b = sum((origin - vertex)**2)**0.5
	c = sum((mobile - origin)**2)**0.5
	
	return np.arccos( ( a**2 + b**2 - c**2 ) / (2 * a * b) ) 
 


#### Function def's end  ########

# python ~/bin/analyzeSplayGeo.py ~/splayBundle/ df1L13G_targetPair001.pdb
# python ~/bin/analyzeSplayGeo.py ~/splayBundle/df1L13G-DIMER.pdb ~/splayBundle/pairs_Cluster-006 ~/splayBundle/mtchInfo-C6-v2.pkl ~/bin/PS
# Path to centered dimer, path to aligned PDBs, path to pickled hash of alignment info, path to PS binary

# import target helices B and B from different chains (resi 27-28), target from alignment (centered 0,0,0,)
target 		= parsePDB( sys.argv[1], subset = 'bb')
mtchDict 	= pic.load( open( sys.argv[3], 'rb' ) )


# Calculate angle of splayed axis from matched helix and 0,0,0 
# by finding points where helical axes of incident and target helices intersect XY-plane (z = 0)
cmdWTsplay 	= [
	 sys.argv[4],
	 '-i', sys.argv[1], 
	 '-f', '5', 
	 '-l', '20',
	 '-c', 'B',
	  ] 
cmdTarget	 = [ 
	 sys.argv[4],
	 '-i', sys.argv[1], 
	 '-f', '32', 
	 '-l', '45',
	 '-c', 'B',
	  ] 

print 
outWTsplay  = sp.Popen( cmdWTsplay  , stdout = sp.PIPE ).communicate()[0]
outTar 		= sp.Popen( cmdTarget  ,  stdout = sp.PIPE ).communicate()[0]
#print outTar 
#axis 	= open( './tmpAxis.pdb', 'w' )
coords 	= []
for i in outTar.split('\n'):
	if i[:6] == 'HETATM':
		coords.append( ( float( i[30:38]  ) , float( i[38:46] ) , float( i[46:54] ) ) )
		#print i 
		#axis.write( i + '\n' )
	if i[:15] == ' Sphere radius:':
		radCurve = int( round( float( i.split()[-2] ), 0 ) )
#axis.close()
targetInt = np.array( XYpoint( coords, radCurve ) )

coords 	= []
for i in outWTsplay.split('\n'):
	if i[:6] == 'HETATM':
		coords.append( ( float( i[30:38]  ) , float( i[38:46] ) , float( i[46:54] ) ) )
	if i[:15] == ' Sphere radius:':
		radCurve = int( round( float( i.split()[-2] ), 0 ) )

WTsplayInt = np.array( XYpoint( coords, radCurve ) )

wtAng = round( angle(  WTsplayInt, targetInt, np.zeros(3) ) * 180/np.pi , 3 )

# Determine if angle is positive or negative

#print calcCenter( target )
refAngle = round( angle(  WTsplayInt, targetInt, calcCenter( target ) ) * 180/np.pi, 3 )

if refAngle > 90:
	wtAng = 360 - wtAng
print 'WT splay angle:', wtAng

angleData = []
step = 0
#angleDataFile = open( 'c6-angleLog.txt', 'w' )
# find helical axis in fit helices
for f in os.listdir( sys.argv[2] ):
	break
	if f[-4:] != '.pdb':
		continue
	path 	= os.path.join( sys.argv[2], f )
	
	aligned = parsePDB( path )
	span 	= mtchDict[ f.split('_aligned')[0] ][-1].split('_') 
	span  	= np.arange( int( span[0] ), int( span[-1] ) + 1 )

	incident 	= aligned.select( 'not resnum %s' % ( ' '.join( [ str(x) for x in span] ) ) )
	inResiRng 	= ( incident.getResnums()[0], incident.getResnums()[-1] )
	chain 		= incident.getChids()[0]

	# find helical axis of incident helix & intersection to z = 0 plane
	cmd = [ sys.argv[4],
	 '-i', path, 
	 '-f', str( inResiRng[0] ), 
	 '-l', str( inResiRng[-1] ),
	 '-c', chain
	  ] 

	out = sp.Popen( cmd , stdout = sp.PIPE ).communicate()[0]

	coords 	= []
	for i in out.split('\n'):
		if i[:6] == 'HETATM':
			try:
				coords.append( ( float( i[30:38]  ) , float( i[38:46] ) , float( i[46:54] ) ) )
			except ValueError:
				break

		if i[:15] == ' Sphere radius:':
			radCurve = int( round( float( i[20:28].strip() ), 0 ) )

	if len(coords) < 6:	 continue

	intersect = XYpoint( coords, radCurve ) 

	if intersect:
		#if radCurve > 25:
			intersect = np.array( XYpoint( coords, radCurve )  )
	else:
		'failed', f
		continue

	# Calculate angle based on XY-plane projections, and note positive negative angles
	ang 	 = round( angle( intersect , targetInt, np.zeros( 3 ) ) * 180/np.pi, 3 )
	refAngle = round( angle( intersect, targetInt, calcCenter( target ) ) * 180/np.pi, 3 )

	if refAngle > 90:
		ang = 360 - ang

	print ang, refAngle, f, radCurve, intersect, step 
	step += 1
	angleData.append( ang )
	angleDataFile.write( str(ang) + '\t' + f + '\n' )
#angleData = np.array( angleData )
#angleDataFile.close()
with open( 'c6-angleLog.txt' ) as file:
	for i in file:
		angleData.append( float( i.split()[0] ) )

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

n, bins, patches = plt.hist( angleData , 50, facecolor='green', alpha=0.75)
plt.xlabel('Angle splayed from bundle')
plt.ylabel('Counts')
plt.title(r'$\mathrm{Histogram\ of\ Splay Angles:}$')
plt.axis([0, 360, 0, 10])
plt.grid(True)

plt.show()

