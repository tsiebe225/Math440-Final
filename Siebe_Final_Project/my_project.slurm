#!/bin/bash -x
#SBATCH --job-name="TS_O8X64"
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks=64
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --account=math540s18
#SBATCH --time=80:00:00
#SBATCH -o conNOO8X64.out
#SBATCH -e conNOO8X64.err

srun  ./mainOpen8_exe
