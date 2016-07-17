# Repack input for variable cccp


import sys, os, subprocess as sp, time
from collections import defaultdict

##
#  python ~/bin/ms1_repack.py ~/tertBuilding/MS1/ms1_trial_1259/ ~/rosetta/ ~/tertBuilding/MS1/Xms1_Spread28.span ~/tertBuilding/MS1/helix_Relax-MS1mini.xml ~/tertBuilding/MS1/resfile5 ~/tertBuilding/MS1/resfile6
##
### I/O

workingDir	= sys.argv[1]

rosiBase 	= sys.argv[2]
rosiDB 		= os.path.join( rosiBase, 'database/')
rosiScrps	= os.path.join( rosiBase, 'source/bin/rosetta_scripts.linuxgccrelease' )
rosiSpanGen = os.path.join( rosiBase, 'source/bin/spanfile_from_pdb.linuxgccrelease')

spanFile 	= sys.argv[3]
designXML 	= sys.argv[4]
resfile5 	= sys.argv[5]
resfile6 	= sys.argv[6]

input_5 	= os.path.join( workingDir, os.path.basename( os.path.dirname( workingDir ) ) + '_5.pdb' ) 
input_6 	= os.path.join( workingDir, os.path.basename( os.path.dirname( workingDir ) ) + '_6.pdb' ) 
input_8 	= os.path.join( workingDir, os.path.basename( os.path.dirname( workingDir ) ) + '_8.pdb' ) 
input_9 	= os.path.join( workingDir, os.path.basename( os.path.dirname( workingDir ) ) + '_9.pdb' ) 

synonym = {}
synonym[5] 	= (resfile5, input_5)
synonym[6] 	= (resfile6, input_6)
synonym[8] 	= (resfile6, input_8)
synonym[9] 	= (resfile6, input_9)

FNULL 		= open( os.devnull, 'w' )	

### HARD CODED
baseSC 		= {}
baseSC['6'] 	= -185.311
baseSC['8'] 	= -185.311
baseSC['9'] 	= -185.311
baseSC['5'] 	= -168.387

os.chdir(workingDir)



### 50 repacks of 9
n = '1'
for k, v in synonym.items():

	cmd = [ rosiScrps, 
'-parser:protocol', designXML, 					# Path to Rosetta script (see above)
'-in:file:s', v[1],						# Input PDB structure
'-nstruct', n, 								# Generate 1000 models
'-mp:setup:spanfiles', spanFile,				# Input spanfile
'-mp:scoring:hbond', 						# Turn on membrane hydrogen bonding
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:prefix',workingDir,
'-packing:resfile', v[0],
'-out:overwrite',
'-packing:pack_missing_sidechains', '0' ]

	start = time.time()
	sp.call( cmd , stdout=FNULL, stderr=sp.STDOUT )
	print cmd
	print 'elapsed:', time.time() - start


	#sys.exit()

## Summarize scores
sumF 	= os.path.join( workingDir, 'score_summary.txt' ) 
scD 	= defaultdict( list )
#for f in [x for x in os.listdir() if len( x.split('_') ) == 4 ] :
for f in os.listdir(workingDir):

	info = f.split('_')
	if len( info ) != 5: continue

	model 	= info[2]
	thrID 	= info[3][0]
	path 	= os.path.join( workingDir, f )
	print path, thrID
	with open( path) as fin:
		for i in fin:

			if i[:5] == 'score':
					scD[thrID].append( float( i.split()[1] ) ) 

outFile = os.path.join( workingDir , '%s_scoreSummary.txt' % model )
txt = ''
for k, v in scD.items():
	txt += '%s %f\n' % (k, min( v ) - baseSC[k])
	print '%s %f' % (k, min( v ) - baseSC[k])
oF = open( outFile, 'w' )
oF.write( txt )



