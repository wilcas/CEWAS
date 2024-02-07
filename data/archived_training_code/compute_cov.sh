#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
source activate stat
python -u make_methy_covar.py
