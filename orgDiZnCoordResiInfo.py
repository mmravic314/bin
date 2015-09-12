# python org

import sys, cPickle as pic

metalDict 	= pic.load( open( sys.argv[1], 'rb') )
validZn 	= [ x.rstrip() for x in open(sys.argv[2], 'rU').readlines() if len(x.rstrip()) ==4  ]

for k, v in metalDict.items():
	if k not in validZn: 
		continue
	print k 
	for p in v:
		#print p, '\t', ' '.join( [ r[2:5] for r in p.contacts] )
		print p, p.contacts	


	print 
