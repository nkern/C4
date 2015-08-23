#!/bin/sh
#PBS -S /bin/sh
#PBS -A christoq_flux
#PBS -q flux
#PBS -N C4_Clusters
#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=5:00:00
#PBS -V
#PBS -o /nfs/christoq_ls/nkern/C4/
#PBS -j oe
#

# Go to Working Directory
cd /nfs/christoq_ls/nkern/C4/

# Run Code 
python c4_stack.py 0

