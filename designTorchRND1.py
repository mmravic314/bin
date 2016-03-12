#!/usr/bin/python

#$ -S /usr/bin/python
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G
#$ -l h_rt=00:29:00
#$ -cwd
#$ -j y
#$ -o /netapp/home/mmravic/tertBuilding/CMP_bobo/redesign_inputs/logs
#$ -t 1-300

#############
#############


### local design
# python ~/bin/designTorchRND1.py ~/binLocal/Rosetta/main/ ~/peptideAmyloid/RND1_designTrj030515/model_643016/model_643016.pdb ~/peptideAmyloid/RND1_designTrj030515/designProtocol_RND1.xml ~/peptideAmyloid/RND1_designTrj030515/resfile 
#

# 1) path 2 rosetta main
# 2) path to input structure file
# 3) path to XML script for protocol
# 4) resfile

import sys, os, subprocess as sp, re, time

start = time.time()

################## MAIN #######################
# Non-variable args
rosetta_database_path   = os.path.join( sys.argv[1] , 'database/' )
rosetta_scriptsEXE_path = os.path.join( sys.argv[1], 'source/bin/rosetta_scripts.linuxgccrelease' )
design_script_path      = sys.argv[3]
struc_path 				= sys.argv[2]
resfile_path			= sys.argv[4]

# Variable args
index = os.path.basename( struc_path[:-4].split('_')[-1] ) 

output_prefix			= os.path.dirname( sys.argv[2] ) + '/'
try:
	output_suffix 			= '_out%s' % (str(  os.environ["SGE_TASK_ID"]) )
except KeyError:
	output_suffix                      = '_out%s' % ( 'Local' )

cmd = [
		rosetta_scriptsEXE_path,
		'-database', rosetta_database_path,
		'-parser:protocol', design_script_path,
		'-in:file:s', struc_path,
		'-out:prefix', output_prefix,   
		'-out:suffix', output_suffix,                               
		'-out:no_nstruct_label',
		'-out:overwrite',
        '-packing:resfile', resfile_path,
]

print
print cmd
print

sp.call( cmd )

print "\n\nWall clock time used:", time.time() - start