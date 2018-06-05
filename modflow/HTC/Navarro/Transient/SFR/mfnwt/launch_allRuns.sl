#!/bin/bash -l
#SBATCH --time=0-05:59:59
#SBATCH --nodes=1
#SBATCH --array=0,1-776:25
#SBATCH --job-name=mf
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1200M
#SBATCH --mail-user=samuelczipper@gmail.com
#SBATCH --mail-type=END
#SBATCH --output=mf%a-%j.out

# change to directory
cd /home/zipper/scratch/TNC-PilotProject/modflow/HTC/Navarro/Transient/SFR/mfnwt/mf$SLURM_ARRAY_TASK_ID

# launch script in directory
./launch_thisRun.sh

