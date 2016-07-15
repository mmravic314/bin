#! /bin/python 
######################
#
## input directory, parse all .par files, write a bash script to run a 'generateCrickBB.m' command
# Use input path to write output files
import sys,os,re

bashStr = '#!/bin/bash\n'

## octave --silent --eval  "fcrick('/Users/mmravic/splayBundle/1_extSearch/3BVX-006_007-0423_0445_A-0454_0469_A_aligned-SYM_24-39.pdb',4,'GENERAL',1,'/Users/mmravic/splayBundle/1_extSearch/3BVX-006_007-0423_0445_A-0454_0469_A_aligned-SYM_24-39',[],[],[],[])"


# function XYZ = generateCrickBB(chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin)
# hard coded to search for 4 chains, that is 3 parameters relative to chain A

for f in [ os.path.join( sys.argv[1], x ) for x in os.listdir( sys.argv[1] ) if x[-4:] == '.par' ]:
	rF =  open( f, 'rU').read() 

	print "Entering...", f, '\n'

	match = re.search( r'(R0.*)\n(R1.*)\n(w0.*)\n(w1.*)\n', rF)
	print match.group()

	print 
	print rF

	sys.exit()
