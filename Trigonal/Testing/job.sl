#!/usr/bin/tcsh
#SBATCH --job-name=test    # Job name
#SBATCH  --ntasks 1 --cpus-per-task 1
#SBATCH --mem=100mb                     # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec

cd ~/MonteCarlo/Trigonal/Testing/

./make


