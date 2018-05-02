#!/bin/bash -l
#SBATCH --time=0-0:30:00
#SBATCH --nodes=1
#SBATCH --array=60-126
#SBATCH --job-name=mf
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --mail-user=samuelczipper@gmail.com
#SBATCH --mail-type=END
#SBATCH --output=mf%a-%j.out

# change to directory
cd /home/zipper/scratch/TNC-PilotProject/modflow/HTC/Navarro/SteadyState/RIV/mfnwt/mf$SLURM_ARRAY_TASK_ID

# launch script in directory
./launch_thisRun.sh

