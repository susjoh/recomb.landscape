#!/bin/sh

#$ -cwd
#$ -l h_rt=01:00:00
#$ -V
#$ -l h_vmem=5200M

. /etc/profile.d/modules.sh
module load R/3.1.1
R CMD BATCH 2_generate_files.R


