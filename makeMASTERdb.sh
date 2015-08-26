# Making a MASTER database (http://grigoryanlab.org/master/)
# Make a 'target'  type of each file as well as 'query' type
# just run inside the working directory (e.g. /home/xray/tertBuilding/pA/)


ls ./ > list
mkdir database
mkdir queries
while read  p
do
	/home/xray/termanal/createPDS --type target --pdb $p --pds ./database/${p:0:-4}.pds --unnaturalAA
	/home/xray/termanal/createPDS --type query --pdb $p --pds ./queries/${p:0:-4}.pds   --unnaturalAA


done < list
