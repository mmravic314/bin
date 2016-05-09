## Main function has option to do the following (has & unhash to do or not)
## Anlayze model designibility scores from peptide design

# input 1: summary txt file of scores
# input 2: directory to put outputs in

## python ~/bin/scores_analyzePlot.py ~/peptideAmyloid/negPaSummary.txt ~/peptideAmyloid/rnd1_analysis/negTheta/

import sys, os, numpy as np, cPickle as pic
from collections import defaultdict
import matplotlib.pyplot as plt

# holds info of each evaluation
class dataPoint:

	def __init__( self, model, params, scores, cls ):
		self.model 	= model
		self.params	= params
		self.scores	= scores
		self.cls 	= cls

	def __repr__( self ):
		return self.model




## print sorted best designibility, return pickle or path to pickle
def sortBest( fPath, outDir ):

	total 	= []
	good 	= {}
	pklPath = os.path.splitext( os.path.basename( fPath ) )[0]
	print
	ind 	= 0 
	with open( fPath ) as fin:
		for i in fin:
			if i[0] == 'm': continue		# skip first line

			#print i.split()[-1]
			params 	= [ float(j) 	for j in i.split()[2:8]  ]
			scores 	= [ int(j) 		for j in i.split()[9:13] ]
			cls 	= int( i.split()[-1] )
			model 	= i.split()[0]

			ind +=1
			if int( i.split()[-1] ):
				if min( scores ) > 140:
					good[model] = dataPoint( model, params, scores, cls )
					print model, ''.join( params ), scores

	print 
	print ind, len( good )
	for k,v in sorted( good.items(), key=lambda x: min( x[1].scores ), reverse=True ):
		pass


	return

#sortBest( sys.argv[1], sys.argv[2] )
#sys.exit()

## plot 1D histograms for each parameter given its classification
## return overall normalized frequencies (obs) each class

tru 	= [] 
fals 	= []

# gather data in each class, to select bins
objs = defaultdict(list)
indDic = { 0:'Bx', 1:'Theta', 2:'Ny', 3:'Zn', 4:'Phi', 5:'W'  }		# probably unnecessary hashing of param labels to parameter array

with open( sys.argv[1] ) as fin:
		for i in fin:
			if i[0] == 'm': continue		# skip first line

			params 	= [ float(j) 	for j in i.split()[2:8]  ]
			ind = 0
			for k in params:
				try:
					if k not in objs[ind]:
						objs[ind].append( k )
				except KeyError:
						objs[ind].append( k )
				ind +=1


			if int( i.split()[-1] ):
				tru.append( params )
			#	print i

			else:
				fals.append( params )


tru 	= np.array( tru )
fals	= np.array( fals )

## make frequency plots of data in panels and overall bar



# plot for True class (good models)
fig, axes 	= plt.subplots(nrows=2, ncols=3)
ax0, ax1, ax2, ax3, ax4, ax5 = axes.flat
setAX 		= [ax0, ax1, ax2, ax3, ax4, ax5]
titles		= [ 'B_x', 'Theta', 'N_y', 'Z_n', 'Phi', 'W_n' ]

for ax, title, k in zip( setAX, titles, [0,1,2,3,4,5] ):

	bin = sorted( objs[k] ) 
	data = tru[:,k]

	ax.hist(  data, bins=bin , normed=1, facecolor='green', alpha=0.75)
	
	ax.set_title( title )


fig.suptitle( 'Rnd2 params for TORCH model hits: class=1 (n=%d)' % ( len(data) ), fontsize=20 )

plt.show()

print 

# plot for False class (good models)
fig, axes 	= plt.subplots(nrows=2, ncols=3)
ax0, ax1, ax2, ax3, ax4, ax5 = axes.flat
setAX 		= [ax0, ax1, ax2, ax3, ax4, ax5]
titles		= [ 'B_x', 'Theta', 'N_y', 'Z_n', 'Phi', 'W_n' ]

for ax, title, k in zip( setAX, titles, [0,1,2,3,4,5] ):

	bin = sorted( objs[k] ) 
	data = fals[:,k]

	ax.hist(  data, bins=bin, normed=1, facecolor='red', alpha=0.75)
	
	ax.set_title( title )


fig.suptitle( 'Rnd2 params for TORCH model misses: class=0 (n=%d)' % ( len(data) ), fontsize=20 )

plt.show()


## pi graph for absolute class variables
# The slices will be ordered and plotted counter-clockwise.
truRate =  float( len( tru ) ) / ( len( tru ) + len( fals ) )
falsR	= 1 - truRate 

labels = 'Pass\nClass=1', 'Fail\nClass=0'
sizes = [truRate, falsR]
colors = ['yellowgreen', 'red']



plt.pie(sizes, labels=labels, colors=colors,
        autopct='%1.1f%%', shadow=True, startangle=20)
# Set aspect ratio to be equal so that pie is drawn as a circle.

plt.title( 'Observed RND2 TORCH Model \n\'Design-ibility\' Classification' , fontsize=20)

plt.axis('equal')
plt.show()



