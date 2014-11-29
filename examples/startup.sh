#!/bin/bash -l

#PBS -N referenceImplementation
#PBS -l nodes=4:ppn=32
#PBS -l walltime=0:05:00
#PBS -q siwir
#PBS -M sebastian.eibl@fau.de -m abe
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err

. /etc/profile.d/modules.sh
module load openmpi/1.6.5-ib
module load gcc/4.8.2

cd ~/SiWiR/siwir_ex3/examples

mpirun -np 128 ./cg 10000 10000 100 -1