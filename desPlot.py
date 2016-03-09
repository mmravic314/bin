#!/bin/python

#python ~/bin/desPlot.py scoreDict.pkl

# suite to analyze and visualize designability data in certain ways
import os, sys, cPickle as pic, numpy as np

import matplotlib.mlab as mlab
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

sDict = pic.load( open( sys.argv[1], 'rb') )


# sorted 

scores 	= []
numRes	= []

params 	= []
goodPs	= []

Beta 	= []
Theta	= []
N_t		= []
Z_n		= []
Z_c		= []
W_n		= []

gBeta 	= []
gTheta	= []
gN_t	= []
gZ_n	= []
gZ_c	= []
gW_n	= []

for num, dat in sorted( sDict.items(), key=lambda n: (  n[1][0][1], n[1][0][0], n[1][0][2]  ) , reverse=True):
	
	score, numR, duds = dat[0][0], dat[0][1], dat[0][2]

## rnd1 design filter
#	if numR < 8 or score < 2300 or duds != 0: continue
	#scores.append( dat[0][0] )
	#numRes.append( dat[0][1] )

	#print num, score, numR, duds
	#x 	= np.array( scores )
	#y 	= np.array( numRes )

## rnd 1 end ##

## plot params of all models, compared to 'okay to good' models
	
	beta, theta, n_t, z_n, z_c, w_n = tuple( dat[1] ) 
	Beta.append(beta)
	Theta.append(theta)
	N_t.append(n_t)
	Z_n.append(z_n)
	Z_c.append(z_c)
	W_n.append(w_n)

	params.append( dat[1] )

	if numR < 6 or score < 1300 or duds > 1: continue
	goodPs.append( dat[1] )

	# grab params for each... all 1D 

	gBeta.append(beta)
	gTheta.append(theta)
	gN_t.append(n_t)
	gZ_n.append(z_n)
	gZ_c.append(z_c)
	gW_n.append(w_n)



## params analysis end 


print '\ntrials:', len(params), ';  good:', len( goodPs )


Beta 	= np.array(Beta)
Theta	= np.array(Theta)
N_t		= np.array(N_t)
Z_n		= np.array( Z_n )
Z_c		= np.array( Z_c )
W_n		= np.array( W_n )
params 	= [ Beta, Theta, N_t, Z_n, Z_c, W_n ] 


gBeta 	= np.array( gBeta )
gTheta	= np.array( gTheta )
gN_t	= np.array( gN_t )
gZ_n	= np.array( gZ_n )
gZ_c	= np.array( gZ_c )
gW_n	= np.array( gW_n )
gparams 	= [ gBeta, gTheta, gN_t, gZ_n, gZ_c, gW_n ]



pSet = [ np.array( [ t, g ] ) for t,g in zip( params, gparams )  ]


#print [np.random.randn(n) for n in [10000, 5000, 2000]]

#sys.exit()

## make overlapping histograms for each set of distrbutions for that param set

fig, axes 	= plt.subplots(nrows=2, ncols=3)
ax0, ax1, ax2, ax3, ax4, ax5 = axes.flat
setAX 		= [ax0, ax1, ax2, ax3, ax4, ax5]
label		= ['total', 'select']
titles		= [ 'Beta', 'Theta', 'N_t', 'Z_n', 'Z_c', 'W_n' ]

for t, g, ax, title in zip( params, gparams, setAX, titles ):

	ax.hist( [t, g], 8, histtype='step', label=label)
	
	ax.set_title( title )


fig.suptitle( 'Rnd1 params for TORCH models; total: blue (n=13194), select: green (n=1727)' , fontsize=20 )

plt.show()



sys.exit()

####### 2D histogram of designability & number of residues ######

# define the colormap
cmap = plt.cm.jet
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
cmaplist[0] = (.5,.5,.5,1.0)
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)


fig, ax = plt.subplots()
H = ax.hist2d(x, y, bins=12, cmap = cmap)
fig.colorbar(H[3], ax=ax)

plt.title( 'TORCH Designability via TERMS (n=13194)', fontsize =24 )
plt.ylabel('Number of TERMS found', fontsize=14)
plt.xlabel('Design-ability score', fontsize=14)

plt.show()
