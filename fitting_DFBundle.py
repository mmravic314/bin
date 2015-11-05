
import sys, os, subprocess as sp, cPickle as pic, numpy as np, time
from prody import *
from PDButil import *

# input target PDB, path to list .txt file with one path per PDB to align to target, path to runCE.sh (CE must be in original directory)

#  python ~/bin/fitting_DFBundle.py ~/splayBundle/df1L13G_targetPair001CENz.pdb ~/splayBundle/Cluster-006 ~/bin/protein-comparison-tool_20130326/runCE.sh 

target 	 		= parsePDB( sys.argv[1], subset = 'bb', chain = 'B' )
#targetFull  	= parsePDB( sys.argv[1], subset = 'bb')


#Centering with coiled-coil axis aligned to Z axis... only do this once
#moveAtoms( targetFull.select( 'resnum %s' % (' '.join( [str(x) for x in np.arange(27, 49)] ) ) ) , to = np.zeros(3), ag = True )
#writePDB( sys.argv[1][:-4] + 'X.pdb' , targetFull )

#sys.exit()

targPath = sys.argv[1][:-4] + 'chB.pdb'

writePDB( targPath, target )


#targetPath 	= os.path.join( sys.argv[2], os.path.basename( sys.argv[1] ) )
#targetID	= os.path.basename( sys.argv[1] )[:-4]
matchDir 	= os.path.join( os.path.dirname( sys.argv[1] ), 'matches_%s/' % ( sys.argv[2].split('/')[-1]  ) )
fragDir 	= os.path.join( os.path.dirname( sys.argv[1] ), 'pairs_%s/'   % ( sys.argv[2].split('/')[-1] ) )

if not os.path.exists( matchDir ):
	os.mkdir( matchDir )

if not os.path.exists( fragDir ):
	os.mkdir( fragDir )


### commented out section used to write mobile PDB paths to list
listF = open( 'list_%s.txt' % ( sys.argv[2].split('/')[-1] ) , 'w' )
for f in os.listdir( sys.argv[2] ):
		mobilePath 	= os.path.join( os.path.abspath( sys.argv[2] ), f )
		print mobilePath
		listF.write( mobilePath + '\n')

		# Renumber/rewrite all the helical pairs in the cluster from 1 to n
		# instead of starting from 1 on each chain
		pdb 	= parsePDB( mobilePath )
		length 	= 1 
		for r in pdb.iterResidues():
			r.setResnum( length )
			length += 1
		writePDB( mobilePath, pdb )
#sys.exit()
time.sleep(2)


## Do CE backbone alignments to DF helix B. Write PDB files in match directory if length is > 10 & RMSD < 2
## Remember matched region in hash
import re
#logf 	= open( '/Users/mmravic/splayBundle/matchLog.txt', 'a' )
out  	= ''
cutCnt 	= 0
mtchDict = {}  #		= pic.load( open( '/Users/mmravic/splayBundle/mtchInfo-C6.pkl' ,'rb' ) )

with open( 'list_%s.txt' % ( sys.argv[2].split('/')[-1] ) ) as file:
	for f in file:

		print f.rstrip()
		print 

		#File to write PDB of alignment to DF chain B
		outF = os.path.join( matchDir , os.path.basename( f.rstrip() )[:-4] + '_dfBmatch.pdb' )

		# Full command to output a PDB file and shortened file to just read CE output text
		cmd  	= ['bash', sys.argv[3], '-file1', targPath, '-file2', f.rstrip(), '-outFile', outF,  '-outputPDB', '-printCE']
		#miniCmd = ['bash', sys.argv[3], '-file1', targPath, '-file2', f.rstrip(), '-printCE']


		# Run in directory with CEAlign Java files
		os.chdir( os.path.dirname( sys.argv[3] ) )
		out = sp.Popen( cmd  , stdout = sp.PIPE ).communicate()[0]			#  full command
		#out = sp.Popen( miniCmd  , stdout = sp.PIPE ).communicate()[0]		#  abbreviated command

		match = re.search( r'Alignment\slength\s=\s(\d+)\sRmsd\s=\s(\d\.\d\d)A' , out)
		if not match:
			AlignLength, RMSD = 0.0, 0.0
			print "Alignment error... Cut!!!!!!!!!!!!!"
			continue
		else:
			AlignLength, RMSD = int( match.group(1) ) , float( match.group(2) )
		print AlignLength, RMSD

		# Remove poor matches
		if AlignLength < 10 or RMSD > 2:
			try:
				os.remove( outF )
			except OSError:
				pass
			cutCnt += 1
			print "Cut!!!!!!!!!!!!!"
			continue

		# Store alignment information	
		match2 = re.search( r'Chain\s1:\s+(\d+)\s.+\nChain\s2:\s+(\d+)\s', out )
		#print out
		print match2.group(1), match2.group(2)

		# Record which regions matched ( format of original PDB file )
		key 	= os.path.basename( f.rstrip() )[:-4]
		ch1		= '%s_%s' % ( match2.group(1), str( int(match2.group(1)) + AlignLength - 1 ) ) 
		ch2		= '%s_%s' % ( match2.group(2), str( int(match2.group(2)) + AlignLength - 1 ) ) 
		matchVal= ( ch1, ch2 )			# store tuple of residue range for each chain e.g. ( 2_22, 29_49 )
		mtchDict[ key ] = matchVal		# Indexed by mobile's pdb name, without the file extension


pic.dump( mtchDict, open( '/Users/mmravic/splayBundle/mtchInfo-C%s.pkl' % ( sys.argv[2][-1] ) ,'wb' ) )


# After writing good matches to PDB (See below), clear any that clash with the DF A helix 
targetPair		= parsePDB(  os.path.join( '/Users/mmravic/splayBundle', sys.argv[1]) , subset = 'bb', chain = 'A' )		# Df1 helix A, not the matched template
mtchDict		= pic.load( open( '/Users/mmravic/splayBundle/mtchInfo-C%s.pkl' % ( sys.argv[2][-1] ),'rb' ) )
fastaMatch 		= open( sys.argv[1][:-4] + '.fasta', 'a' )



# Write aligned fasta files and aligned PDB files containing only relevant regions
for f in os.listdir( matchDir ):

	print f,
	path 	= os.path.join( matchDir, f )
	mobile 	= parsePDB( path, subset = 'bb', model = 2 )
	target 	= parsePDB( path, subset = 'bb', model = 1 )

	# Recall residues aligned for this match
	segTargX, segMobileX = mtchDict[ f.split( '_dfBmatch.pdb' )[0] ]
	
	# renumber residues to match numbering of pdb file
	segTarg 	= np.arange( int(segTargX.split('_')[0] ) - 1, int(segTargX.split('_')[1] )  ) + target.getResnums()[0]
	
	# Get fasta from match and atom selection string of subset of aligned atoms/residues
	segMobile 	= np.arange( int(segMobileX.split('_')[0] ) - 1, int(segMobileX.split('_')[1] )  ) + mobile.getResnums()[0]
	mobSelStr 	= 'resnum '


	fasta		= []
	step = 1
	chainMobile = []
	for x in mobile.iterResidues():
		if step in segMobile:
			fasta.append( natAA[ x.getResname() ] )
			mobSelStr += '%s ' % ( x.getResnum() )
			if x.getChid() not in chainMobile:
				chainMobile.append( x.getChid() )
		step += 1
	mobSelStr += 'chain %s' % ( chainMobile[0] ) 
	fasta = ''.join(fasta)

	# get object containing just aligned residues from mobile and target selection strings
	targSelStr 	= 'resnum %s' % ( ' '.join( [ str( x ) for x in segTarg ] ) )

	targMtch 	= target.select( targSelStr )
	mobMtch		= mobile.select( mobSelStr  )

	mobileLabels = []

	# From mobile, include all residues with backbone atoms within 15 Angstroms from any atom in aligned helix
	otherChainMobile = ''
	if chainMobile[0] == 'A':
		otherChainMobile = 'B'
	else:
		otherChainMobile = 'A'
	contact = mobile.select( 'backbone chain %s and within 15 of %s' % ( otherChainMobile, mobSelStr ) )

	## Disclude any residues with backbone atoms 7 Angstroms from DF helix A
	clashing= mobile.select( 'backbone chain %s and within 7 of targetPair' % ( otherChainMobile), targetPair = targetPair   )
	try:
		clashing = [ ( a.getResnum(), a.getChid() ) for a in clashing.iterAtoms() ]

	except AttributeError:
		clashing = []


	# Record all resdiues from aligned helix and 'close' residues on incident helix form mobile, to rewrite PDB
	try:
		for a in contact.iterAtoms():
			ID = ( a.getResnum(), a.getChid() )
			if ID in clashing:
				continue
			if ID not in mobileLabels:
				mobileLabels.append( ID )
	except AttributeError:
		print
		continue

	# if interacting portion of indicent helix length <8 residues, forget alignment and match 
	if len( mobileLabels ) < 12: 
		print 
		continue

	# If eligible residues from incident helix not continuous (check previous Res number), forget alignment and match
	breakFlg 	= False
	first		= 0
	for label in mobileLabels:
		if first ==0:
			first 		+= 1	
			prvResnum 	= int( label[0] )
			continue

		if int( label[0] ) - 1 != prvResnum:
			print '\tDiscontinuous!'
			breakFlg = True
			break 

		prvResnum 	= int( label[0] )



	if breakFlg:
		continue


	# Else, get remainder of ID for helix pair fragment to be written in new PDB
	for a in mobile.select( mobSelStr ).iterAtoms():
		ID = ( a.getResnum(), a.getChid() )
		if ID not in mobileLabels:
			mobileLabels.append( ID )

	# Rewrite PDB from aligned match file including only 'interacting' region of incident helix & DF1-aligned region
	#print mobileLabels
	inFile 	= open( path, 'rU' )
	outFile	= open( os.path.join( fragDir, f.split( '_dfBmatch.pdb' )[0] + '_aligned.pdb' ) , 'w' )
	modelFlg= 0
	with open( path ) as file:
		for i in file:

			if i[:5] == 'MODEL':
				if i.rsplit()[-1] == '2':
					modelFlg += 1

			if modelFlg == 0:
				continue

			try:
				if ( int(i[22:26].strip()), i[21:22].strip() ) in mobileLabels:
					outFile.write( i )
			except ValueError:
				continue


	# fasta of mobile sequence aligned to target; initialize as '...' length of target and replace aligned region
	# Write to file
	alignedFasta= ''
	index = 0
	for a in target.iterResidues():
		if a.getResnum() in segTarg:
			alignedFasta += fasta[index]
			index += 1
		else: 
			alignedFasta += '.'
	alignedFasta = '%s \t%s' % ( f.split( '_dfBmatch.pdb' )[0] , alignedFasta )
	fastaMatch.write( alignedFasta + '\n')


	print mtchDict[ f.split( '_dfBmatch.pdb' )[0] ], f.split( '_dfBmatch.pdb' )[0]
	print alignedFasta

	print '\n'
	
fastaMatch.close()
sys.exit()



