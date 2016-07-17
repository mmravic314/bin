#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /netapp/home/mmravic/tertBuilding/MS1/logs                       #-- output directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=1G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=1G,scratch=1G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=0:29:59                #-- runtime limit (see above; this requests 24 hours)
##$ -t 1-2448                        #-- remove first '#' to specify the number of
                                   #-- tasks if desired (see Tips section)


inx=1

#taskID=$SGE_TASK_ID
taskID=2448

for dir in /home/xray/tertBuilding/MS1/ms1_trial_*/
	
	do 
	
	# stop and execute command when on the proper parameter set
	if [ $inx -eq $taskID ]
			then

				python ~/bin/ms1_repack.py $dir ~/rosetta/ ~/tertBuilding/MS1/Xms1_Spread28.span ~/tertBuilding/MS1/helix_Relax-MS1mini.xml ~/tertBuilding/MS1/resfile5 ~/tertBuilding/MS1/resfile6
				break

			fi
			inx=$(( inx + 1))
	
	#python ~/bin/ms1_repack.py $dir ~/rosetta/ ~/tertBuilding/MS1/Xms1_Spread28.span ~/tertBuilding/MS1/helix_Relax-MS1mini.xml ~/tertBuilding/MS1/resfile5 ~/tertBuilding/MS1/resfile6
	done
