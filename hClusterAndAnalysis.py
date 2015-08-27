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
loopUp	= pic.load( open( sys.argv[2], 'rb'	) 	)

clustering = heirarchy_cluster( dMatrix, 2.5 )

x = 0
from collections import defaultdict
cluster = defaultdict(list)
while x < len(clustering):
	#clustering[x], loopUp[ str(( x, x )) ][0]
	cluster[clustering[x]].append( loopUp[ str(( x, x )) ][0] )

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