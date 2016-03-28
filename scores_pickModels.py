## Main function has option to do the following (has & unhash to do or not)
## Anlayze model designibility scores from peptide design

# input 1: summary txt file of scores
# input 2: directory to put outputs in

## python ~/bin/scores_analyzePlot.py ~/peptideAmyloid/negPaSummary.txt ~/peptideAmyloid/rnd1_analysis/negTheta/

import sys, os, numpy as np, cPickle as pic

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
				if min( scores ) > 100:
					good[model] = dataPoint( model, params, scores, cls )

	print '#model_ID min_frag_hits params_Bx-Ny_Theta_Zn_Phi_Wn' 
#	print ind, len( good )
	for k,v in sorted( good.items(), key=lambda x: min( x[1].scores ), reverse=True ):
		print v, min( v.scores ), ' '.join( [ str( b ) for b in v.params ] )

	print
	return

## plot 1D histograms for each parameter given its classification
## return overall normalized frequencies (obs) each class

# 

##

########## Main #########


## in the shell command, pipe the output to a design model 

# python ~/bin/scores_pickModels.py ~/peptideAmyloid/negPaSummary.txt ~/peptideAmyloid/rnd1_analysis/negTheta/ > models2designRND1.txt

sortBest( sys.argv[1], sys.argv[2] )

