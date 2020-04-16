#!/bin/bash
#
#SBATCH -p shared	 				# Partition to submit to (comma separated)
#SBATCH -J prospyall					# Job name
#SBATCH -n 1						# Number of cores
#SBATCH -N 1						# Ensure that all cores are on one machine
#SBATCH -t 0-4:00					# Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 2000					# Memory in MB
#SBATCH -o Logs/.%j.out					# file for standard out writing
#SBATCH -e Logs/.%j.err					# file for standard err writing
#SBATCH --mail-type=ALL					# Type of email notification: BEGIN,END,FAIL,ALL
#SBATCh --mail-user=locke.patton@cfa.harvard.edu	# email for notifications
#SBATCH --array=0-31					# iterations of same script

i=$SLURM_ARRAY_TASK_ID

my_output_directory=/n/holystore01/LABS/berger_lab/Lab/lockepatton

cd $my_output_directory

echo Starting postcreatefast.py ...
echo $SLURM_ARRAY_TASK_ID

python ./prospectTools/postcreatefast.py $SLURM_ARRAY_TASK_ID
