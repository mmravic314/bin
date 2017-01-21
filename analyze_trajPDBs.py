# Do analysis on 100 output models from peptide amyloid tranjectories
# plot those properties

# input 1: PDB file to align everything to, equilibrated final frame
# input 2: stack of MD files, each as a model

# in ACMD, make sure all waters
# 
# 

import sys, os, numpy as np
from prody import * 
from collections import defaultdict, Counter


### thing to analyze and plot
#
# 1. per helix CA distance plot v time
# 2. per helix backbone/CA RMSD v time
# 3. per lysine contcants 		v time 
# 4. SASA per strand hydrophobics 
# 5. Energy plots 
# 
### these are done in separate functions 

# notes
# Helices P13, P14, P15.... these are chains  -, -, - respectively. 1376-2361

#python ~/bin/analyze_trajPDBs.py EQ_reference.pdb traj_frames.dcd 

# input a model PDB file, return an array 

plotz = os.path.abspath( 'plots' )
if not os.path.exists( plotz ):
	os.mkdir( plotz )

# Align trajectories to the middle strands calphas
refPDB		= parsePDB( sys.argv[1] )
refStrands 	= refPDB.select( 'segment P3 P5 P7 calpha' )
refHelices	= [ refPDB.select( 'segment P13 ca' ), refPDB.select( 'segment P14 ca' ), refPDB.select( 'segment P15 ca' )]

# for some reason, his has different names in different runs, so note this
hisName 	= refPDB.select('segment P1 calpha').copy()[4].getResname()


ensemble 	= parseDCD( sys.argv[2] )
ensemble.setAtoms( refStrands )
ensemble.superpose()

### work through each model and do RMSD calculations or function to do analysis of each feature

## Based on strand-aligned helices, per helix RMSD
helixRMSD = []
for h in refHelices:
	ensemble.setAtoms( h )
	helixRMSD.append( ensemble.getRMSDs() ) 



## In each frame, calc Lysine and histidine contacts & SASA of apolars in strand
traj = Trajectory( sys.argv[2] )
traj.link( refPDB )
traj.setCoords( refPDB )

# package selections for important his and lys atoms into a hash
std_groups	= {}

LysHyd		= [ refPDB.select( 'segment %s name HZ1 HZ2 HZ3' % a ) for a in [ 'P3', 'P5', 'P7', 'P9'] ] 
for k, v in zip( [ 'S3', 'S5', 'S7', 'S9'], LysHyd):
	std_groups[ k + ' Lys_amine' ] = v

HisH		= [ refPDB.select( 'segment %s name HD1 and (resname %s) and (resnum 14)' % ( a, hisName) ) for a in [ 'P3', 'P5', 'P7', 'P9'] ]
if hisName == 'HSE':
	HisH		= [ refPDB.select( 'segment %s name HE1 and (resname %s) and (resnum 14)' % ( a, hisName) ) for a in [ 'P3', 'P5', 'P7', 'P9'] ]
for k, v in zip( [ 'S3', 'S5', 'S7', 'S9'], HisH):
	std_groups[ k + ' His_NH' ] = v


HisN		= [ refPDB.select( 'segment %s name NE2 and (resname %s) and (resnum 14)' % ( a, hisName) ) for a in [ 'P3', 'P5', 'P7', 'P9'] ]
if hisName == 'HSE':
	HisN		= [ refPDB.select( 'segment %s name ND1 and (resname %s) and (resnum 14)' % ( a, hisName) ) for a in [ 'P3', 'P5', 'P7', 'P9'] ]
for k, v in zip( [ 'S3', 'S5', 'S7', 'S9'], HisN):
	std_groups[ k + ' His_N' ] = v

#


# Package selections for VAL to grab SASA of each over time
#apolr_grps	= {}

#val12 		= [ 'segment %s resname VAL resnum 12' % a for a in [ 'P3', 'P5', 'P7', 'P9'] ]
#for k, v in zip( [ 'S3', 'S5', 'S7', 'S9'], val12):
#	apolr_grps[ k + ' Val_12' ] = v

#val18		= [ 'segment %s resname VAL resnum 18' % a for a in [ 'P3', 'P5', 'P7', 'P9'] ]
#for k, v in zip( [ 'S3', 'S5', 'S7', 'S9'], val18):
#	apolr_grps[ k + ' Val_18' ] = v
#

SASA_calc	= defaultdict( list )
cnt_Counter = defaultdict( list )
tmpPath 	= os.path.join( os.getcwd() , 'wtvr.pdb' )
dssPath		= os.path.join( os.getcwd() , 'wtvr.dssp' )

# hard coded hash to link dssp numbering to apolar valine groups of interest
apolr_grps	= {}
apolr_grps['15'] = 'S3 Val_12'; apolr_grps['21'] = 'S3 Val_18';
apolr_grps['27'] = 'S5 Val_12'; apolr_grps['33'] = 'S5 Val_18';
apolr_grps['39'] = 'S7 Val_12'; apolr_grps['45'] = 'S7 Val_18';
apolr_grps['51'] = 'S9 Val_12'; apolr_grps['57'] = 'S9 Val_18';

name = os.path.basename(  os.getcwd() )


for i, frame in enumerate(traj):

	## calc SASA for each strand, writing a DSSP to a tmp file and pulling the data
	writePDB( tmpPath, refPDB )
	execDSSP( tmpPath )
	
	startFlg = 0
	with open( dssPath ) as fin:
		for i in fin:
			num = i.split()[0]
			if num == '#':
				startFlg = 1
				continue
			if startFlg:
				try:
					SASA_calc[ apolr_grps[ num ] ].append( int( i[35:38].strip() ) )
				except KeyError:
					continue


	## for each contact of interest record number of atoms within range... maxed out at 1 or 3, for his and lys respectively
	for k, v in std_groups.items():

		hYes 		= 0	
		# collect ID info
		if 'His' in k:
			badRes 	= hisName
			maxV	= 1 
			# look for hydrogen only contacts ... for NE2. default = FALSE 
			if 'His_NH' not in k:
				hYes= 1
		else:
			badRes 	= 'LYS'
			maxV 	= 3


		# Get raw contacts to heavy atoms
		if hYes:
			atoms 	= refPDB.select( '(within 2.7 of pdb) and (not resname %s) and hydrogen' % badRes, pdb=v )
		else:
			atoms 	= refPDB.select( '(within 2.7 of pdb) and (not resname %s) and heavy' % badRes, pdb=v )


		# Store valid counts to array
		cnts 	= []
		if not atoms:
			cnt_Counter[ k ].append( 0 )
		#	print k, 0, '\n'
			continue

		# ignore atoms from same residue/water and carbons/sulfurs
		for a in atoms:
		
			if 'C' in a.getName() or 'S' in a.getName(): continue
			if a.getResnum() not in cnts and len(cnts) < maxV:
				cnts.append( a.getResnum() )

		
		#print k, len( cnts )
		# add to hash
		cnt_Counter[ k ].append( len( cnts ) )
		#print 
	##


import matplotlib.pyplot as plt

# plotting section
print 'helix RMSD'
num_frames 	= len( helixRMSD[0] )
time 		= np.arange( 1, num_frames + 1 )
stp = 1
color = 'obgr'

for h in helixRMSD:
	plt.plot(time, h, '.%s-' % color[stp], label='H-%d'% stp )
	#plt.show()
	stp +=1
plt.xlabel('Frame Number, (100 ps)')
plt.ylabel('RMSD (Angstrom)') 
plt.title( 'TORCH Helices C-alpha RMSD, 10 ns\n%s' %  name  )
plt.axis( [0, 100, 0, 4.5] )
plt.legend( ncol=3, mode="expand", borderaxespad=0.)
path = os.path.join( plotz, 'helix_rmsd.pdf' )
plt.savefig( path, bbox_inches='tight' ) 
plt.close()
print

print 'CONTACTS'

stp 	= 1
color 	= { 'S3':'b', 'S5':'g', 'S7':'r', 'S9':'k' }

for k, v in sorted( cnt_Counter.items() ):

	clr = color[ k.split()[0] ]
	ax 	= plt.subplot(4,3, stp)
	plt.plot(time, v, '.%s-' % clr) 

	if stp in [10,11,12]:
		plt.xlabel('Frame Number, (100 ps)')
	else:
		plt.setp(ax.get_xticklabels(), visible=False)
	if stp in [1,4,7,10]:
		plt.ylabel('# Contacts')
		plt.yticks( [0.0, 1.0, 2.0, 3.0] )
	else:
		plt.setp(ax.get_yticklabels(), visible=False)

	plt.title( k, fontsize=12 )
	#print k, stp
	stp += 1
	plt.axis( [0, 100, 0, 3.1] )
	
plt.suptitle( 'Contacts to fiber chemical groups, 10 ns\n%s' % name, fontsize=16  )
path2 = os.path.join( plotz, 'contacts.pdf' )
plt.tight_layout()
plt.subplots_adjust(top=0.85)
plt.savefig( path2, bbox_inches='tight' ) 
plt.close()

print 'SASA'
stp = 1
for k, v in sorted( SASA_calc.items() ):

	clr = color[ k.split()[0] ]
	ax 	= plt.subplot(4,2, stp)
	plt.plot(time, v, '.%s-' % clr) 

	if stp in [7,8]:
		plt.xlabel('Frame Number, (100 ps)')
	else:
		plt.setp(ax.get_xticklabels(), visible=False)
	if stp in [1,3,5,7]:
		plt.ylabel('SASA ($\AA^2$)')
		plt.yticks( [0.0, 40.0, 80.0, 120.0] )
	else:
		plt.setp(ax.get_yticklabels(), visible=False)

	plt.title( k, fontsize=12 )
	#print k, stp
	stp += 1
	plt.axis( [0, 100, 0, 121] )
	
plt.suptitle( 'apolar SASA of fiber VAL, 10 ns\n%s' % name, fontsize=16  )
path3 = os.path.join( plotz, 'apolar_SASA.pdf' )
plt.tight_layout()
plt.subplots_adjust(top=0.85)
plt.savefig( path3, bbox_inches='tight' ) 
plt.close()
print


