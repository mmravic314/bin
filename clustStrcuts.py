## Marco Mravic UCSF Biophysics DeGrado Lab August 2015 ##
# Use master to write PDB files aligning MASTER searches to 

## NOt-sustainable way to read clustering, but saved std out from hClusterAndAnalysis.py
## Parse into a hash with list of pathes to each term indexed by the centroid path

import os, cPickle as pic, sys, subprocess as sp

dMatrix = pic.load( open( '/home/xray/tertBuilding/pD/distMatrixRedund.pkl', 'rb' ) )	
lookUp	= pic.load( open( '/home/xray/tertBuilding/pD/distMatrixRedund.pkl__lookUpHash', 'rb'	) 	)

clustering = {136: 11, 191: 7, 226: 6, 59: 5, 92: 5, 111: 5, 157: 5, 170: 5, 185: 5, 229: 5, 103: 4, 179: 4, 198: 4, 20: 3, 50: 3, 54: 3, 78: 3, 82: 3, 97: 3, 120: 3, 148: 3, 163: 3, 211: 3, 1: 2, 5: 2, 21: 2, 24: 2, 28: 2, 29: 2, 32: 2, 34: 2, 36: 2, 42: 2, 46: 2, 65: 2, 66: 2, 69: 2, 72: 2, 75: 2, 79: 2, 94: 2, 99: 2, 104: 2, 107: 2, 115: 2, 117: 2, 123: 2, 132: 2, 142: 2, 144: 2, 147: 2, 165: 2, 167: 2, 171: 2, 189: 2, 201: 2, 206: 2, 223: 2, 236: 2, 0: 1, 2: 1, 3: 1, 4: 1, 6: 1, 8: 1, 9: 1, 10: 1, 11: 1, 12: 1, 13: 1, 14: 1, 15: 1, 16: 1, 17: 1, 22: 1, 23: 1, 25: 1, 26: 1, 30: 1, 31: 1, 35: 1, 37: 1, 39: 1, 40: 1, 41: 1, 44: 1, 45: 1, 48: 1, 49: 1, 61: 1, 62: 1, 63: 1, 64: 1, 68: 1, 71: 1, 74: 1, 77: 1, 81: 1, 86: 1, 88: 1, 89: 1, 90: 1, 93: 1, 98: 1, 114: 1, 121: 1, 122: 1, 126: 1, 127: 1, 128: 1, 129: 1, 150: 1, 151: 1, 154: 1, 159: 1, 161: 1, 169: 1, 184: 1, 195: 1, 200: 1, 203: 1, 212: 1, 215: 1, 216: 1, 218: 1}


fclus = open( '/home/xray/tertBuilding/pD/clusterOutPut.txt', 'rU' )

centroidDict = {}
for i in fclus:
	#print i.rstrip() 
	if len( i.rstrip() )  == 0: continue

	if i.rstrip()[-1] == '/':
		clust = i.rstrip().split('/')[-2].split('_')[-1]
		
		continue
	if i[:7] == 'Cluster':
		centroidDict[ clust ] = int( i.split( )[1] )

masterPath = '/home/xray/termanal/master'
subdirPath = '/home/xray/tertBuilding/pD'
for f in sorted( [c for c in os.listdir(subdirPath) if 'Cluster_' in c ], key = lambda x: int( x.split('_')[-1] ) ) :
	membersList  	= os.listdir( os.path.join( subdirPath, f ) )

	if len( membersList ) > 1:
		clus 		=  f.split('_')[-1] 
		centroid 	= lookUp[  str( ( centroidDict[ clus ], centroidDict[ clus ] ) ) ][0] 

		centroidPdsPath = os.path.join( os.path.join( subdirPath, 'queries/'), centroid + '.pds' )
		print f, 'centroid:', centroid

		for p in membersList:
			if p[0] =='m': continue

			if p[:-4] != centroid:
				memberPdsPath 	= os.path.join( os.path.join( subdirPath, 'database/'), p[:-4] + '.pds' )

				outPath		= os.path.join( os.path.join( subdirPath, f), 'm' + p[:-4] + '.pdb' )
			
				cmd = [ masterPath,  '--query', centroidPdsPath, '--target', memberPdsPath, '--structOut' , os.path.join( subdirPath, f), 
						 '--outType', 'full', '--topN', '1', '--rmsdCut', '3.0']
				print cmd
				print 		 
				sp.call(cmd)
				rename = [ 'mv', os.path.join( os.path.join( subdirPath, f), 'full1.pdb' ),   outPath ]
				sp.call(rename)
			
		print 

