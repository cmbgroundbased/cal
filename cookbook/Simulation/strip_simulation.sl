#!/bin/bash -l
#
#SBATCH --job-name=strip_atm
#SBATCH --account=strip
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu 1000
#SBATCH --time=2-00:00:00
#SBATCH --exclusive
#SBATCH --error=job_log/job.%J.err 
#SBATCH --output=job_log/job.%J.out

module load Mlnx/hpcx

cd $HOME/STRIP_PYCAL/cookbook/Simulation

pwd

export OMP_NUM_THREADS=32
singularity exec pycal.SIF python3 strip_pipeline.py @strip_file.par
ls

seff -d $SLURM_JOBID
  
	
