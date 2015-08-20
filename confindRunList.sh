#!/bin/bash 
## Run many confind runs from a list of PDB files

# Prior to running: >> mkdir zN_rOUT/

for line in `cat $1`
do
    #echo ~/bin/confind --p $line --o ./zN_rOUT/${line:30:-4}.cmap --rLib ~/bin/RR2000.rotlib --rout ./zN_rOUT/${line:30:-4}.rout
    
    echo $line
    ~/termanal/confind --p $line --o ./zN_rOUT/${line:30:-4}.cmap --rLib ~/bin/RR2000.rotlib --rout ./zN_rOUT/${line:30:-4}.rout 
    echo 
done