module load dSQ
dsq --job-file job_list.txt --batch-file array_job_submission_script.sh --partition psych_day --nodes 1 --ntasks 1 --mem-per-cpu 16g --cpus-per-task 10 -t 04:00:00 --mail-type ALL
sbatch array_job_submission_script.sh
