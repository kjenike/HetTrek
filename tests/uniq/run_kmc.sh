#!/bin/bash
# specify BASH shell
#$ -S /bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
#  specify that the job requires X GB of memory per thread (25G for kmc, 1G for smoothing)
#$ -l m_mem_free=5G
# specify the number of threads to be used (5 for kmc, 20 for smoothing)




#folder=ERR1.0_HET1.0_INDEL1
#mkdir "$folder"
mkdir SIMULATION
mkdir KMC
mkdir KMC/tmp

#python /seq/schatz/tbenavi/software/HetSmoother/simulator.py -o SIMULATION -g 100000 -c 30 -r 1000 --het_rate 1 -e 1 -i 1 -k 35 -t t2t_chr22_26_100kb.fasta

/seq/schatz/kjenike/software/KMCs/k21/KMC-smooth/bin/kmc -k35 -t1 -m50 -sm -ci1 -cs10000000 -fm SIMULATION/simulatedreads_template1_g5000_cov40_rl1000_err1.0_het1.0_indel1.fasta KMC/g5k_k35 KMC/tmp/

/seq/schatz/kjenike/software/KMCs/k21/KMC-smooth/bin/kmc_tools transform KMC/g5k_k35 histogram KMC/g5k_k35.hist -cx10000000

