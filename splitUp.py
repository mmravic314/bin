# input a text files list of pathes, one per line
# Split into several listXr.txt files. in user specific limits

import sys,os

fCnt		= 1
pathCnt		= 0
txt 		= ''
with open( sys.argv[1] ) as file:
	for i in file:

		txt += i
		pathCnt +=1

		if pathCnt == 80000:
			outFile = open( 'list%dR.txt' % ( fCnt ) , 'w' )
			outFile.write( txt )
			outFile.close()	
			fCnt += 1
			pathCnt = 0
			txt = ''