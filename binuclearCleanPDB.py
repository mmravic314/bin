## Marco Mravic UCSF Biophysics DeGrado Lab August 2015 ##
# Input list of PDB's with PDB chain mappings to UniProt entries
# Download all chains once, comment out that section
# Determine which pdb chains have zinc binding capabilities in nature (despite crystallization metals) from uniprot
# Create hash table out of chain mappings to compare coordination ligands and UniProt descriptions
# Output list of PDBs with natural zinc2+ (and Mg2+?) to cull for seq ID

import PDButil as *  ## this is a homemade module
import prody as *	 ## http://prody.csb.pitt.edu/
import sys, nump as np, cPickle as pic

#Input file format 
# 4BBJ [4BBJa, Q5ZWR1]
# 4ORB [4ORBa, P63328, 4ORBb, Q63810]
# 2HZH []
# 3F79 [3F79a, Q9I045, 3F79b, Q9I045, 3F79c, Q9I045, 3F79d, Q9I045, 3F79e, Q9I045, 3F79f, Q9I045]
