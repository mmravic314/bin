#!/bin/bash

#$ -S /bin/bash
#$ -l mem_free=16G
#$ -l arch=linux-x64
#$ -l netapp=1G
#$ -l h_rt=00:25:00
#$ -cwd
#$ -j y
#$ -o /netapp/home/mmravic/peptideAmyloid/logfiles/

### Know task from number of length of list file or number of pdb's in directory... hardcode in
#$ -t 1-20013

# For each task, read the input file and stop on the line defined by the SGE task ID (on cluster)
# if parameter file has a skip number (e.g. params2.txt starts with 1#), then all model id's should increase by ( skip number * 50000 ) via parameter subset size, from paramsGen.py


### Example command line for submission
#
# qsub ~/bin/torch_evalQsub.sh ~/peptideAmyloid/parameterization/params5.txt ~/peptideAmyloid/OsakaModels/ ~/peptideAmyloid/scores_rnd1/ ~/peptideAmyloid/parameterization/helix_prep/18_Ideal_allALA.pdb ~/peptideAmyloid/parameterization/FullInterfaceCENTERED.pdb ~/termanal/support.default/152901_bc_30-scPDB/list.txt ~/termanal/
#
####

## actual
# qsub ~/bin/torch_evalQsub.sh ~/peptideAmyloid/parameterization/params12.txt ~/peptideAmyloid/OsakaModels/ ~/peptideAmyloid/scores_rnd1/ ~/peptideAmyloid/parameterization/helix_prep/18_Ideal_allALA.pdb ~/peptideAmyloid/parameterization/FullInterfaceCENTERED.pdb ~/binLocal/termanal/support.default/162901_bc_30-scPDB/list.txt ~/binLocal/termanal/
##

## 
# global params for i/o

paramsF=$1
modelDir="$2"
scoresDir="$3"
helixF=$4
fibrilF=$5
dbListF=$6
termDir=$7

echo
echo recieved input files/directories:
echo params file: $paramsF
echo model directory $modelDir
echo \'design-ibility\' scores directory $scoresDir
echo \'ideal\' helix file: $helixF
echo target fibril file: $fibrilF
echo MASTER database list: $dbListF
echo directory with master/CreatePDS: $termDir
echo 

echo
echo Begin process...
echo
echo
echo

echo Reading params \& generating model
echo
##



header="$(head -n 1 < $1)"

hsh=$((` expr index "$header" \#` -1 ))		# the index of the hash character in the first line of the params file
skipNum=${header::$hsh}


# subtracting accounts for the first line that is skipped
inx=$(($skipNum*20013 + 1)) 
start=$inx

# 

# tasks start at 1, so {1,2,...50000}, {50001,...100000}
taskID=$SGE_TASK_ID 			#CLUSTER
#taskID=5						# LOCAL 
stop=$(( $start + $taskID - 1 ))

basename="model_$stop.pdb"
base2="tmp_model_$stop"
modelPath="$modelDir$basename"
fragDir="$modelDir$base2"/



while read -r line
	do 
		# skip first line, by only operating when second place in string of $line is length=0 or not equal to '#'
		if [[ "${line:$hsh}" != '# dBeta Theta N_t Z_n Z_c W_n' ]]
		then
			# stop and execute command when on the proper parameter set
			if [ $inx -eq $stop ]
			then

				params=$line
				break

			fi
			inx=$(( inx + 1))
		fi


done < $1



echo params $params 

#### make model
echo
echo Making $modelPath
## local
#python ~/bin/params2coords.py "$params" $helixF $modelPath $fibrilF


## Cluster


scl enable python27 "python ~/bin/params2coords_rnd2.py $params $helixF $modelPath $fibrilF"


####

echo 
echo "attempting fragment generation from model..."
echo
##### try to generate fragments

## local
#python ~/bin/torch_desginibility.py $modelPath $fragDir $scoresDir


## Cluster
scl enable python27 "python ~/bin/torch_desginibility_rnd2.py $modelPath $fragDir $scoresDir"

#####
echo
echo "attempting to find 'design-ibility' score via searching frags (if present)..."
echo
###### try to send MASTER jobs to assess fragment design-ibility


## local
#python ~/bin/torch_master.py $fragDir $scoresDir $termDir $dbListF

## Cluster
scl enable python27 "python ~/bin/torch_master.py $fragDir $scoresDir $termDir $dbListF"


######

gzip -f $modelPath

echo
echo "gzipped/compressed model file"
echo
echo
echo
echo
echo ...process Terminated
echo


