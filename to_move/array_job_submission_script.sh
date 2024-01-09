#!/bin/bash
#SBATCH --array 0-98
#SBATCH --output dsq-job_list-%A_%2a-%N.out
#SBATCH --job-name dsq-job_list
#SBATCH --partition psych_day --nodes 1 --ntasks 1 --mem-per-cpu 16g --cpus-per-task 10 -t 04:00:00 --mail-type ALL

# DO NOT EDIT LINE BELOW
/gpfs/milgram/apps/hpc.rhel7/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/milgram/pi/chang/pg496/repositories/eisg/scripts/job_list.txt --status-dir /gpfs/milgram/pi/chang/pg496/repositories/eisg/scripts

