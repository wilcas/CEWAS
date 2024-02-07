#!/bin/bash
#PBS -q small
#PBS -l walltime=2:00:00,mem=22gb,nodes=1:ppn=4
#PBS -m beaf
#PBS -M william.casazza@stat.ubc.ca

set -e -x -o pipefail

Rscript train_cmc_enet.R

