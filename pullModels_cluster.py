import sys, os, subprocess as sp


## given list of got torch models, write their full remote path and tell bash to download to cetrain location

# python ~/bin/pullModels_cluster.py ~/peptideAmyloid/models2designRND1.txt ~/peptideAmyloid/rnd1_design/ mmravic@pass2.compbio.ucsf.edu/netapp/home/mmravic/peptideAmyloid/rnd1OsakaModels/


with open( sys.argv[1] ) as fin:
	for i in fin:

		if len(i) < 10: continue

		SSHpath = os.path.join( sys.argv[3], 'model_%s.pdb.gz' % i.split()[0] )
		locPath =  sys.argv[2]

		cmd = [ 'scp', SSHpath, locPath ]

		print
		print cmd
		print

		sp.call( cmd )
