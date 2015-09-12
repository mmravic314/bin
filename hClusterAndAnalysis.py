# Marco Mravic UCSF DeGrado Lab August 2015 
# Cluster and analyze TERM fragments with pre-defined RMSD similarity matrix
# Plot & print info about clusters and members
## input: 1) path to distance matrix
## 2) path to indexing between matrix and strucutral pairs 
## 3) path to directory full of .pdb files that were clusters
## output filepath for cluster ID's (shuttles copy of pdb files into subdirectory specified for file)
# example cmd line 
# python ~/bin/hClusterAndAnalysis.py ~/tertBuilding/pA/distMatrixRedund.pkl ~/tertBuilding/pA/distMatrixRedund.pkl__lookUpHash ~/tertBuilding/pA/ ~/tertBuilding/pA/clustering.txt


from PDButil import *
import cPickle as pic, sys,os, shutil

dMatrix = pic.load( open( sys.argv[1], 'rb' ) 	)
lookUp	= pic.load( open( sys.argv[2], 'rb'	) 	)

# heirarchical clustering with flattening for clusters at rmsd cutoff
#clustering = heirarchy_cluster( dMatrix, 2.5 )

# kmedioids clustering with number of clusters chose as smallest with mean RMSD to medoid as threshold
kClustering = kmedoid_clustering( dMatrix, 2.5, lookUp)

from collections import defaultdict
cluster = defaultdict(list)
x = 0

for p in kClustering[0]:
	ID 			= lookUp[ str( ( p,x ) )  ]
	centroid 	= lookUp[ str( ( p,p ) )  ][0]
	print ID, centroid
	if ID[0] == ID[-1]:
		job = centroid
	else:
		job = ''.join( [ g for g in list(ID) if g != centroid ] )
	
	
	cluster[p].append( job )
	x += 1

cnum = 1
outFile = open( sys.argv[4], 'w') 

for k , v in cluster.items():
	path = os.path.join( os.path.dirname( sys.argv[3] ), 'Cluster_%s/' % ( str(cnum) ) )
	print path
	if not os.path.exists (path ):
		os.mkdir(path)
	for p in sorted(v):
		outFile.write( 'Cluster: %s %s\n' % (k, p) )
		print 'Cluster: %s %s\n' % (k, p)
		srcPath = os.path.join( sys.argv[3], '%s.pdb' % (p) )
		newpath = os.path.join( path, '%s.pdb' % (p) )
		print newpath
		shutil.copy( srcPath, newpath )
	cnum += 1

sys.exit()











while x < len(clustering):
	#clustering[x], loopUp[ str(( x, x )) ][0]
	cluster[clustering[x]].append( lookUp[ str(( x, x )) ][0] )

	x += 1 

outFile = open( sys.argv[4], 'w')
for k , v in cluster.items():
	path = os.path.join( os.path.dirname( sys.argv[3] ), 'Cluster_%s/' % (k) )
	print path
	if not os.path.exists (path ):
		os.mkdir(path)
	for p in sorted(v):
		outFile.write( 'Cluster: %s %s\n' % (k, p) )
		print 'Cluster: %s %s\n' % (k, p)
		srcPath = os.path.join( sys.argv[3], '%s.pdb' % (p) )
		newpath = os.path.join( path, '%s.pdb' % (p) )
		#print newpath
		shutil.copy( srcPath, newpath )


#from scipy.cluster import hierarchy as h
#import matplotlib.pyplot as plt


#matrix = np.random.rand(6,6)



#tree = h.linkage( matrix, method='ward', metric='euclidean')
#print matrix
#print 'linkage'
#print tree

#h.dendrogram( tree )
#plt.show()
#import scipy.spatial.distance as sp_dist
#import scipy.cluster.hierarchy as sp_clust
''' xray@dgl-xray:~/tertBuilding$ python ~/bin/hClusterAndAnalysis.py ~/tertBuilding/pD/distMatrixRedund.pkl ~/tertBuilding/pD/distMatrixRedund.pkl__lookUpHash ~/tertBuilding/pD/ ~/tertBuilding/pD/clustering.txt
397.188693047 elapsed for  20 clusters 1 identical trajectories from total 1142420 max RMSD of 3.5596
512.422659874 elapsed for  25 clusters 1 identical trajectories from total 1428025 max RMSD of 2.8576
629.151689053 elapsed for  30 clusters 1 identical trajectories from total 1713630 max RMSD of 2.8576
753.630767822 elapsed for  35 clusters 1 identical trajectories from total 1999235 max RMSD of 2.8576
24.18572402 elapsed for  50 clusters 1 identical trajectories from total 57121 max RMSD of 3.4672
26.345138073 elapsed for  75 clusters 1 identical trajectories from total 57121 max RMSD of 2.8032
28.7958381176 elapsed for  100 clusters 1 identical trajectories from total 57121 max RMSD of 2.7308

284.655477047 elapsed for  100 clusters 1 identical trajectory(s) from total 571210 max RMSD of 2.9772
298.134126902 elapsed for  120 clusters 1 identical trajectory(s) from total 571210 max RMSD of 2.816
302.209667206 elapsed for  125 clusters 1 identical trajectory(s) from total 571210 max RMSD of 2.6916
Counter({111: 8, 148: 7, 136: 6, 157: 6, 193: 6, 96: 5, 170: 5, 226: 5, 82: 4, 103: 4, 141: 4, 166: 4, 180: 4, 198: 4, 3: 3, 33: 3, 54: 3, 65: 3, 92: 3, 120: 3, 151: 3, 163: 3, 186: 3, 220: 3, 6: 2, 16: 2, 19: 2, 28: 2, 34: 2, 40: 2, 42: 2, 50: 2, 52: 2, 61: 2, 66: 2, 69: 2, 72: 2, 75: 2, 78: 2, 79: 2, 90: 2, 94: 2, 104: 2, 117: 2, 123: 2, 128: 2, 142: 2, 144: 2, 167: 2, 173: 2, 179: 2, 187: 2, 189: 2, 201: 2, 206: 2, 211: 2, 216: 2, 221: 2, 234: 2, 236: 2, 0: 1, 1: 1, 4: 1, 5: 1, 7: 1, 8: 1, 9: 1, 10: 1, 11: 1, 12: 1, 13: 1, 14: 1, 15: 1, 18: 1, 22: 1, 23: 1, 24: 1, 25: 1, 26: 1, 27: 1, 30: 1, 31: 1, 35: 1, 36: 1, 37: 1, 38: 1, 39: 1, 44: 1, 45: 1, 46: 1, 47: 1, 48: 1, 49: 1, 53: 1, 55: 1, 57: 1, 59: 1, 62: 1, 64: 1, 71: 1, 74: 1, 77: 1, 86: 1, 88: 1, 89: 1, 106: 1, 121: 1, 122: 1, 126: 1, 131: 1, 135: 1, 137: 1, 149: 1, 150: 1, 153: 1, 154: 1, 160: 1, 195: 1, 203: 1, 204: 1, 205: 1, 209: 1, 210: 1, 217: 1, 233: 1})
3995.31619501 elapsed for  125 clusters 1 identical trajectory(s) from total 7140125 max RMSD of 2.6392




313.125841856 elapsed for  150 clusters 1 identical trajectory(s) from total 571210 max RMSD of 2.816
310.919182062 elapsed for  150 clusters 1 identical trajectory(s) from total 571210 max RMSD of 2.3872

Counter({68: 7, 136: 7, 157: 7, 226: 6, 131: 5, 180: 4, 188: 4, 193: 4, 36: 3, 92: 3, 96: 3, 111: 3, 120: 3, 151: 3, 198: 3, 211: 3, 220: 3, 3: 2, 12: 2, 25: 2, 29: 2, 32: 2, 44: 2, 52: 2, 54: 2, 69: 2, 70: 2, 75: 2, 78: 2, 79: 2, 94: 2, 102: 2, 104: 2, 107: 2, 115: 2, 117: 2, 123: 2, 128: 2, 142: 2, 148: 2, 165: 2, 167: 2, 170: 2, 171: 2, 173: 2, 179: 2, 189: 2, 206: 2, 221: 2, 231: 2, 234: 2, 236: 2, 0: 1, 1: 1, 2: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, 11: 1, 13: 1, 15: 1, 16: 1, 17: 1, 18: 1, 19: 1, 20: 1, 22: 1, 23: 1, 24: 1, 27: 1, 28: 1, 30: 1, 31: 1, 34: 1, 35: 1, 38: 1, 39: 1, 40: 1, 41: 1, 42: 1, 43: 1, 46: 1, 47: 1, 48: 1, 49: 1, 50: 1, 51: 1, 53: 1, 55: 1, 56: 1, 57: 1, 59: 1, 61: 1, 62: 1, 63: 1, 64: 1, 71: 1, 72: 1, 73: 1, 74: 1, 77: 1, 81: 1, 82: 1, 83: 1, 84: 1, 85: 1, 86: 1, 87: 1, 88: 1, 89: 1, 90: 1, 93: 1, 98: 1, 100: 1, 101: 1, 106: 1, 108: 1, 109: 1, 114: 1, 121: 1, 122: 1, 125: 1, 127: 1, 133: 1, 146: 1, 147: 1, 149: 1, 150: 1, 153: 1, 154: 1, 155: 1, 156: 1, 159: 1, 174: 1, 182: 1, 195: 1, 196: 1, 203: 1, 205: 1, 212: 1, 215: 1, 216: 1, 218: 1, 219: 1, 230: 1, 233: 1})



xray@dgl-xray:~/tertBuilding$ python ~/bin/hClusterAndAnalysis.py ~/tertBuilding/pE/distMatrixRedund.pkl ~/tertBuilding/pE/distMatrixRedund.pkl__lookUpHash ~/tertBuilding/pE/ ~/tertBuilding/pE/clustering.txt
4332.21388912 elapsed for  125 clusters 1 identical trajectory(s) from total 7442000 max RMSD of 2.9788
(array([  0,   1,   2,   3,   4,   5,   5,   7,   8,   9,  10,  14,  14,
        13,  14,  15,  16,  17,  18,  19,  18,  19,  18,  23,  24,  25,
        26,  27,  27,  26,  30,  31,  32,  33,  34,  35,  36,  37,  38,
        39,  40,  40,  42,  43,  43,  45,  46,  47,  48,  49,  49,  51,
        52,  57,  54,  54,  57,  57,  54, 224,  60,  61,  62,  63,  64,
        64,  66,  67,  66,  69,  70, 202, 202,  73,  74,  75,  76,  77,
        79,  79,  79,  79,  82,  30,  84,  85,  86,  87,  88,  88,  90,
        61,  92,  93,  95,  95,  95,  97, 101, 101, 101, 101, 102, 103,
       102, 105, 107, 107, 108, 108, 107, 105, 105, 113, 113, 115, 116,
       118, 118, 119, 120, 120, 122, 123,  87, 125, 125, 127, 128, 133,
       130, 128, 128, 133, 128, 133, 133, 128, 128, 139, 139,  95, 143,
       143, 144, 145, 146,  76, 148, 149, 149, 161, 161, 153, 153, 153,
       153, 153, 153, 153, 153, 161, 153, 167, 164, 165, 166, 167, 167,
        61, 170,  31, 176, 176, 176, 176, 176, 176, 178, 128, 180, 180,
       180, 180, 176, 185, 185, 185, 185, 189, 189, 189, 167, 178, 196,
       196, 196, 196, 198,  61, 200, 201, 202, 202,  23, 206, 206, 165,
         1, 209, 209, 118,  39, 213, 213, 213, 213,  93, 220, 220, 220,
       123, 222, 224, 224,   3, 228, 228, 228, 228, 228, 231, 231, 224,
       206, 143, 143, 237, 237, 239, 239, 241, 241, 243], dtype=int32), 10.485799999999994, 1)

'''
'''
Counter({136: 11, 191: 7, 226: 6, 59: 5, 92: 5, 111: 5, 157: 5, 170: 5, 185: 5, 229: 5, 103: 4, 179: 4, 198: 4, 20: 3, 50: 3, 54: 3, 78: 3, 82: 3, 97: 3, 120: 3, 148: 3, 163: 3, 211: 3, 1: 2, 5: 2, 21: 2, 24: 2, 28: 2, 29: 2, 32: 2, 34: 2, 36: 2, 42: 2, 46: 2, 65: 2, 66: 2, 69: 2, 72: 2, 75: 2, 79: 2, 94: 2, 99: 2, 104: 2, 107: 2, 115: 2, 117: 2, 123: 2, 132: 2, 142: 2, 144: 2, 147: 2, 165: 2, 167: 2, 171: 2, 189: 2, 201: 2, 206: 2, 223: 2, 236: 2, 0: 1, 2: 1, 3: 1, 4: 1, 6: 1, 8: 1, 9: 1, 10: 1, 11: 1, 12: 1, 13: 1, 14: 1, 15: 1, 16: 1, 17: 1, 22: 1, 23: 1, 25: 1, 26: 1, 30: 1, 31: 1, 35: 1, 37: 1, 39: 1, 40: 1, 41: 1, 44: 1, 45: 1, 48: 1, 49: 1, 61: 1, 62: 1, 63: 1, 64: 1, 68: 1, 71: 1, 74: 1, 77: 1, 81: 1, 86: 1, 88: 1, 89: 1, 90: 1, 93: 1, 98: 1, 114: 1, 121: 1, 122: 1, 126: 1, 127: 1, 128: 1, 129: 1, 150: 1, 151: 1, 154: 1, 159: 1, 161: 1, 169: 1, 184: 1, 195: 1, 200: 1, 203: 1, 212: 1, 215: 1, 216: 1, 218: 1})

'''