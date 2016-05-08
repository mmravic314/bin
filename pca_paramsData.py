# Marco Mravic

# cluster successful models by parameters to grab representitive

import sys, os, numpy as np, cPickle as pic

# input paramsSummary file
# input path to pickle of extracted params scores

# example command line
# python ~/bin/pca_paramsData.py ~/peptideAmyloid/rnd2_Summary.txt ~/peptideAmyloid/rnd2_Summary_GOOD.pkl ~/peptideAmyloid/rnd2_clusters.txt ~/peptideAmyloid/rnd2_OsakaModels/

# load in data, save in Pkl

def load_Data( inF ):
	pDict = {}

	with open( inF ) as fin:
		for l in fin:

			if l[0] == 'm': continue

			k 		= l.split()
			scores 	= [ float(n) for n in k[9:13] ]

			if min( scores ) < 120: continue


			model 	= k[0]
			params 	= ' '.join( k[2:8] )
			#print model, params, scores

			pDict[ params ] = ( model, scores )


	return pDict

##### MAIN #####

# load input data, or make if first time in script
if not os.path.exists( sys.argv[2] ):
	pDict = load_Data ( sys.argv[1] )
	pic.dump( pDict, open( sys.argv[2], 'wb' ) )
else:
	pDict = pic.load( open( sys.argv[2], 'rb' ) )


## Run PCA to visualize clusters and detect variance

from sklearn.decomposition import PCA
from sklearn.preprocessing import scale, StandardScaler
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans, DBSCAN,MeanShift, estimate_bandwidth
#from itertools import cycle
from collections import Counter

# convert params data to matrix, store array of models indexed the same as the params data matrix
data 		= []
modelLookup = []
stp 		= 0
for p in sorted( pDict.keys() ):
	data.append( [ float( n ) for n in p.split() ] )
	modelLookup.append( ( pDict[p], p ) )


data 	= np.array( data )
ss_data = StandardScaler().fit_transform( data )
n_samples, n_features = ss_data.shape
print '\nsamples:', n_samples, 'features:', n_features


## Clustering
## Parameter scan for band width
# Mean Shift

pca 			= PCA(n_components=n_features)
X 				= pca.fit_transform( ss_data )
#X 				= ss_data
quant 			= 0.041
bandwidth 		= estimate_bandwidth( X , quantile=quant )
ms 				= MeanShift(bandwidth=bandwidth, bin_seeding=True)
ms.fit(X)
labels 			= ms.labels_
cluster_centers = ms.cluster_centers_
labels_unique 	= np.unique(labels)
n_clusters_	 	= len(labels_unique)


def rnd( lst ):
	return [ round( n, 5 ) for n in lst ]


# define cluster centers, and their model ID
# this is done by slow crappy comparison since the cluster centroid indices are not given
ind 	= 0
recall 	= {}
smples 	= X[:]


print 'clusters found', n_clusters_
print '\nhere are the centroids'
print
print cluster_centers

cnt = Counter( labels )
print 


print '\nhere are the points closest to each centroid\n'

recall 		= []
ind 		= 0
outTXT		= ''

for cen in cluster_centers:

	#print cen
	dif 	= [ np.linalg.norm( cen - k  ) for k in X ]
	minV 	= min(dif)
	index	= dif.index( minV )
	print 'cluster', ind, 'closest to sample', index, 'model', modelLookup[index]
	path 	= os.path.join( sys.argv[4], 'model_%s' % ( modelLookup[index][0][0] ) )

	outTXT	+= '%s cluster-%d members-%d %s\n' % ( path, ind, cnt[ind], modelLookup[index][1] )

	print

	ind += 1

print 
print outTXT
print
#for k, v in sorted( recall.items() ):
#	print k
#	print v
#	print

#print len( recall.keys() ), n_clusters_

sys.exit()

# plot
colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(X[my_members, 0], X[my_members, 1], col + '.')
    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
plt.title('Mean Shift estimated Clusters: %d, Samples: %d, Features: %d' % ( n_clusters_, n_samples, n_features) )
plt.xlabel( 'PC1' )
plt.ylabel( 'PC2' )


plt.show()

sys.exit()

## Parameter scan for band width
qSet = np.arange( 0.01, 0.05, 0.001 )

for q in qSet:
	print q,
	bandwidth 		= estimate_bandwidth( ss_data , quantile=q )

	ms 				= MeanShift(bandwidth=bandwidth, bin_seeding=True)
	ms.fit( ss_data )
	labels 			= ms.labels_
	cluster_centers = ms.cluster_centers_
	labels_unique 	= np.unique(labels)
	n_clusters_	 	= len(labels_unique)

	print n_clusters_


sys.exit()

## PCR dimensionality reduction noot useful, eigen vectors all about 20-10%

pca = PCA(n_components=n_features)
data_reduced = pca.fit_transform( ss_data )

print 'eigen values', pca.explained_variance_ratio_
print
print 



fig = plt.figure(1, figsize=(8, 6))
ax = Axes3D(fig, elev=-150, azim=110)
X_reduced = PCA(n_components=6).fit_transform( sc_data )
ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2],
           cmap=plt.cm.Paired)
ax.set_title("First three PCA directions")
ax.set_xlabel("1st eigenvector")
ax.w_xaxis.set_ticklabels([])
ax.set_ylabel("2nd eigenvector")
ax.w_yaxis.set_ticklabels([])
ax.set_zlabel("3rd eigenvector")
ax.w_zaxis.set_ticklabels([])

plt.show()





