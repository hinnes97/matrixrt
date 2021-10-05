#!/bin/bash
#SBATCH --output=timestep_test.out
#SBATCH --output=timestep_test.err
#SBATCH --partition=priority-rp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=1000

main.exe
