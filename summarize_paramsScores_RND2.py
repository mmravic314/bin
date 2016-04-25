import sys, os, numpy as np, gzip

#hLen 	= 25.71 ## Global parameter, describing helix used, calculation for PHI ... round 1
hLen	= 30.216


#
#scl enable python27 'python ~/bin/summarize_paramsScores.py ~/peptideAmyloid/rnd1_scores/ ~/peptideAmyloid/rnd1OsakaModels/  rnd2_Summary.txt '
#
# scl enable python27 'python ~/bin/summarize_paramsScores.py ~/peptideAmyloid/tmp_PA/ ~/peptideAmyloid/tmp_PA/ rnd2_Summary.txt '

# open score files and model files to read parameters, write into a summary file (one of two actually, depending on sign of Theta parameter)

# open the two files to write info into
sumF = open( sys.argv[3], 'w'  )

#negF = open( sys.argv[4], 'w'  )
txt = 'modelID |params| dB_x N_y Theta Z_n Phi W_n |scores| HH SH1 SH2 SH3 | class\n'

sumF.write( txt )

# search score files
for i in os.listdir( sys.argv[1] ):
	
	if i[-2:] != 'sc':
		print 'odd file found, skipping', i
		print 
		continue

	scPath = os.path.join( sys.argv[1], i )

	#print scPath
	# fetch designibility scores
	with  open(scPath) as fin:
		for l in fin:
			if l[1] == '#':	continue
			scores = l.split()
 	
	

	scoreV = [ int( x ) for x in scores ]
	#print scores
	scoreSTR = ' '.join( scores )
	
	# fetch params
	modelP = os.path.join( sys.argv[ 2 ], i[:-3] + '.pdb.gz' )
	if not os.path.exists( modelP ):
		print 'missing', modelP
		continue
	with gzip.open( modelP, 'rb' ) as fin:
		for k in fin:
			#paramSTR = ' '.join( k.split()[2:]  )
			#print k, k.split()[3]
			theta = float( k.split()[3]  )
			params = [float(n) for n in k.split()[2:]]
			phi = round(  np.arcsin( ( params[3] - params[4] ) / hLen ), 3 )
			paramSTR = params[:4]
			paramSTR.append( phi )
			paramSTR.append( params[5] )
			paramSTR = ' '.join([str(m) for m in  paramSTR] )

	#		print params, phi, paramSTR
	
			break

	#print theta, paramSTR
	modID = i[6:-3]
	
	# initial classification: all > 50, class = 1
	cls = '0' 
	if np.mean( scoreV ) >= 100 or min( scoreV ) >= 80  :
		cls = '1'
	
	txt = '%s | %s | %s | %s\n' %  ( modID, paramSTR, scoreSTR, cls) 
#	if theta > 0:
#		posF.write( txt  )	
#	else:
	sumF.write( txt  )	
#	print txt
	


