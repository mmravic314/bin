#!/usr/bin/python

import os, sys, subprocess as sp, time
from prody import * 
## Marco Mravic DeGrado lab Dec 2015 UCSF 
## input 1: FULL/ABSOLUTE PATH to directory containing backbone template files
## input 2: path to directory with both createPDS and master binaries
## input 3: path to text file list of pathes to TM-helix-only pds (database)
## input 4: path to output text list: path & number of matches below 1 A to each fragment/interface

## SAMPLE
## python ~/~/designInterfaceExample/support.default//breakUpInterfaces.py ~/designInterfaceExample/ ~/designInterfaceExample/support.default/ ~/designInterfaceExample/support.default/TM_Database/list_tm.txt


## Take interface files
## convert to PDS and submit each to master search of transmembrane data base
## record number of matches under 1 A for each fragment in .match files


step = 0
for filename in os.listdir( sys.argv[1] ):
	print filename






