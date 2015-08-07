#!bin/bash

#DeGradoLab rotation Marco Mravic Fall 2014 Biophysics UCSf
#Input a list
#downloads pdbs into working directory

while read line 
do
    wget http://www.rcsb.org/pdb/files/$line.pdb

done < $1
