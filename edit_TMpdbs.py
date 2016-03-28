
## following pairIntructions.txt, with 


# parse TM definintions file and rewrite fastas

import os, sys, re, cPickle as pic
from prody import *
from collections import defaultdict, Counter
from PDButil import UnNatAA


# xray@dgl-xray:~/jan2016TMpdb$ python ~/bin/breakup_TMpdbs.py pdbDimerList.txt extendedDefinitions.txt selectedOPM/ TM_redundant/ TM_redundant_Fasta


def get_exp_info( listFile, pklPath, headDir ):
	
	infoDir = {}
	total  	= 0
	good 	= 0
	with open( listFile ) as fin:
		for i in fin:
			print i.rstrip()
			headerP = os.path.join( headDir, '%s.pdb' % i.rstrip() )

			## skip if no header found
			if not os.path.exists( headerP ):
				continue

			res = 4
			dat = False
			# find resolution and exp data
			with open( headerP ) as fin2:
				for k in fin2:
					
					if k[:6] == 'EXPDTA':
						if k[6:].strip() == 'X-RAY DIFFRACTION': 

							dat = True
						else: 
							print k[6:].strip()
							break
					
					if k[:3] != 'REM': continue

					if 'REMARK   2 RESOLUTION.' in k: 
						
						res = float( k.split()[3].rstrip() )

			if dat and res <= 3.2:
				good +=1
				infoDir[i.rstrip()] = res

			total += 1

			print

	print good, total	

	pic.dump(  infoDir, open( pklPath, 'wb') )
	return

#get_exp_info( sys.argv[1], sys.argv[2], sys.argv[3] )
infoDir = pic.load( open( sys.argv[6], 'rb' ) )

sys.exit()

# extendeddefinintions.txt
storeDict 		= defaultdict(list)
storeDictFull 	= defaultdict(list)
with open( sys.argv[1] ) as file:
	for i in file:

		# skip ineligibles
		try: 
			infoDir = i.rstrip() 
				
		#write prody selection string if segment > 8, hash by pdbID
		if int( i[18:23].strip() ) - int( i[13:18].strip() ) > 8:

			strk = i[5], i[13:18].strip(), i[18:23].strip()
			storeDict[ i[:4] ].append( 'chain %s resnum %s to %s'% ( i[5], i[13:18].strip(), i[18:23].strip() ) )

		storeDictFull[ i[:4] ].append( 'chain %s resnum %s to %s'% ( i[5], i[13:18].strip(), i[18:23].strip() ) )

motif 	= defaultdict(list)
tm_info = {} 

# access each tertiary pdb, make TM only regions
pdbList = []
with i in ...
		print "Entering", i.rstrip(), '\n'

		pdbPath = os.path.join( sys.argv[3], i.rstrip() + '.pdb' )
		inPDB 	= parsePDB(  pdbPath )

		oPdb	= os.path.join( sys.argv[4], i.rstrip() + '.pdb' )
		oFasta	= os.path.join( sys.argv[5], i.rstrip() + '.fasta' )

		tmPDB 	= False

		chASCII	= 97
		prvCh	= ''
		txt 	= ''

		# skip beta barrels
		if len( storeDictFull[ i.rstrip() ] ) == 0:
			continue

		for sel in storeDictFull[ i.rstrip() ]:

			# grab each
			obj 	= inPDB.select( sel ).copy()
			#print obj
			fasta 	= ''.join( [ UnNatAA[x.getResname()] for x in inPDB.select( sel ).copy().iterResidues() if x.getResname() in UnNatAA.keys()] )
			ch 		= obj[0].getChid()

			# reset the chain to this pseudo chain per TM pass 
			if ch 	== prvCh:
				chASCII += 1
			else:
				chASCII = 97

			newCh 	= ch + chr( chASCII )

			obj.setChids( [ newCh for j in obj.iterAtoms()  ] )
			txt += '>%s%s\n%s\n' % ( i.rstrip(), newCh, fasta )

			if not tmPDB:
				tmPDB = obj

			else: tmPDB += obj

			prvCh 	= ch

			#find right-parallel GAS-x3-gas-x3
			#match = re.search( r'[GAS]\w\w\w[GAS]\w\w\w[GAS]', fasta  )
			#if match:
				# remove lo complexity 
			#	if 'GGG' in match.group() or 'AAA' in match.group() or 'G' not in match.group(): continue 

			#	print match.group(), obj[0].getChid(), i.rstrip()


			#	motif[ i.rstrip() ].append( match.group() )

		# write PDB 
		writePDB( oPdb, tmPDB )
		oF = open( oFasta, 'w' )
		print txt
		oF.write( txt )
		#	print fasta, obj[0].getChid()

		#print 

		#sys.exit()		

for k, v in motif.items():

	if len(v) == 0: continue
	c = Counter( v )


	# find single pass multimers only
	if len( c.keys() ) > 1 : continue
	if c.most_common()[0][1] < 2: continue 

	print k, c.most_common()[0][0], c.most_common()[0][1]

	print 





