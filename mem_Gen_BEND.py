## Marco Mravic DeGrado Lab April 2016
## 

# input 1: Input file,  
# input 2: Outfile file, pdb with DUM atoms as layers

##  python ~/CHAMP/bin/mem_Gen.py > ~/tertiaryBuilding/mem25.pdb


import sys, os, numpy as np
from itertools import product
# parse input

txt = ''

## length of 25
span = np.arange( -24.0, 26.0, 0.5 ) 

lines = ''

normal_in 	= [0.0, 0.0, 1.0]
normal_out 	= [0.0, 0.0, -1.0]
center 		= [0.0, 0.0, 0.0]
d_out		= 14.0
d_in		= 14.0

max_dZ 		= 10
y_cen	 	= 10.0

step = 4000
for k in product( span , repeat = 2 )  :
	


	# given x and y, solve for z with plane equation
	x, y = k[0], k[1]

	# for modification bending
	#if 0.0 <= x  or np.fabs( y ) < 10: continue
	if x >= 0.0: continue
	if x < -20: continue
	if x >= -12 and np.fabs( y ) < 10: continue

	if x < -12 and np.fabs( y ) < 10: 
		dZ = 0 

	else:
		dZ 		=  ( y / np.fabs( y ) ) * ( max_dZ / ( 1 + 0.5* np.fabs( x ) ) ) * np.exp( -0.2 * ( np.fabs( y ) - y_cen )  )
	#print x, y, dZ

	## Inner 
	z_o =   round( dZ + ( -1 * d_out - normal_out[0] * x - normal_out[1] * y ) / normal_out[2], 3 )
	z_i =   round( dZ + ( -1 * d_in - normal_in[0] * x - normal_in[1] * y ) / normal_in[2], 3 )

	## outer 

	x += center[0]
	y += center[1]
	z_o += center[2]
	z_i += center[2]

	txt += 'HETATM %4d  O   DUM  %4d    %8.3f%8.3f%8.3f \n' % ( step, step, x,y,z_o )
	step +=1
	txt += 'HETATM %4d  N   DUM  %4d    %8.3f%8.3f%8.3f \n' % ( step, step, x,y,z_i )
	step +=1 
#print

#sys.exit()
#print
print txt
print
sys.exit()
##Length of 30
span = np.arange( -25.0, 26.0, 2 ) 

lines = ''

normal_in 	= [0.0, 0.0, 1.0]
normal_out 	= [0.0, 0.0, -1.0]
center 		= [0.0, 0.0, 0.0]
d_out		= 15
d_in		= 15

step = 4000
for k in product( span , repeat = 2 )  :
	
	# given x and y, solve for z with plane equation
	x, y = k[0], k[1]

	z_o =   round( ( -1 * d_out - normal_out[0] * x - normal_out[1] * y ) / normal_out[2], 3 )
	z_i =   round( ( -1 * d_in - normal_in[0] * x - normal_in[1] * y ) / normal_in[2], 3 )


	x += center[0]
	y += center[1]
	z_o += center[2]
	z_i += center[2]

	txt += 'HETATM %4d  O   DUM  %4d    %8.3f%8.3f%8.3f \n' % ( step, step, x,y,z_o )
	step +=1
	txt += 'HETATM %4d  N   DUM  %4d    %8.3f%8.3f%8.3f \n' % ( step, step, x,y,z_i )
	step +=1 
print
print txt
print





sys.exit()




help = '''
HETATM 1104 THKN MEM Y  64      18.235  -0.341   0.267  1.00  0.00           X  
HETATM 1105 CNTR MEM Y  64       3.295   0.968   0.005  1.00  0.00           X  
HETATM 1106 NORM MEM Y  64       3.255   0.708   0.970  1.00  0.00           X  

HETATM 4093  N   DUM  4093      22.000   6.000 -17.400                          
HETATM 4093  O   DUM  4093      22.000   6.000  17.400 
'''

memDict = {}
with open( sys.argv[1] ) as fin:
	for i in fin:

		if i[:4] 	=='ATOM':
			txt += i
		elif i[:6] 	== 'REMARK':
			txt += i 
#		elif i[:3] 	== 'TER':
#			txt += i
		elif i[:6]  == 'HETATM':
			memDict[ i.split()[2] ] = np.array( [ float( i[30:38] ), float( i[38:46] ), float( i[46:54] ) ] )
		else:
			pass

#print txt

# calculate unit normal, 

center 		= memDict[ 'CNTR' ]
unitNorm 	= memDict[ 'NORM' ]	- center
thickness	= memDict[ 'THKN' ] - center

#print unitNorm, thickness, np.linalg.norm( unitNorm ), d

point_out 		= unitNorm * 15
normal_out		= point_out + unitNorm
d_out = np.dot( -1* point_out, normal_out )

point_in 		= unitNorm * -15
normal_in		= point_in + unitNorm
d_in = np.dot( -1* point_in, normal_in )

#print d_out, normal_out

span = np.arange( -30.0, 31.0, 2 ) 

lines = ''

step = 4000
for k in product( span , repeat = 2 )  :
	
	# given x and y, solve for z with plane equation
	x, y = k[0], k[1]

	z_o =   round( ( -1 * d_out - normal_out[0] * x - normal_out[1] * y ) / normal_out[2], 3 )
	z_i =   round( ( -1 * d_in - normal_in[0] * x - normal_in[1] * y ) / normal_in[2], 3 )


	x += center[0]
	y += center[1]
	z_o += center[2]
	z_i += center[2]

	txt += 'HETATM %4d  O   DUM  %4d    %8.3f%8.3f%8.3f \n' % ( step, step, x,y,z_o )
	step +=1
	txt += 'HETATM %4d  N   DUM  %4d    %8.3f%8.3f%8.3f \n' % ( step, step, x,y,z_i )
	step +=1 

print txt
print help
