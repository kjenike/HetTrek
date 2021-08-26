#!/bin/bash
# specify BASH shell
#$ -S /bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
#  specify that the job requires X GB of memory per thread (25G for kmc, 1G for smoothing)
#$ -l m_mem_free=2G
# specify the number of threads to be used (5 for kmc, 20 for smoothing)
#$ -pe threads 24


#/seq/schatz/kjenike/software/KMCs/active/KMC-multithreading/bin/smooth 

#/grid/schatz/home/ranallo/software/KMC/bin/smooth


time /seq/schatz/kjenike/software/KMCs/active/KMC-multithreading/bin/smooth -i SIMULATION/simulatedreads_template1_g100000_cov30_rl1000_err1.0_het1.0_indel1.fasta  -j KMC/g100k_k35 -o OUT_100k/ -k 35 -l 12 -t 24





