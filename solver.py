#solver.py

import numpy as np, sys

from math import *

def ang2rad( ang ):
	return ang*np.pi/180.0

def rad2ang( ang ):
	return ang*180/np.pi

def fun( r ):
	return sin(  0.5* ( ang2rad(120) - ( 2* asin( 4.7/r ) ) )     ) - 4.0/r

def fun2( r ):
	return 2*asin( 4.7/r )

print 65.32
print 65.32 + 


sys.exit()

trial = np.arange( 8.6, 8.8, 0.000001 )

minV = 10
minT = 10
for i in trial:
	v = np.fabs( fun(i) )
	print  rad2ang( fun2( i ) ), i,  v
	

	if v < minV:
		minV = v
		minT = i

print minV, minT, rad2ang( fun2( minT ) )

## SOLUTION
#   e = 1.00839172013e-08, Ro= 8.70938199992, phi0_B =  65.3193507677