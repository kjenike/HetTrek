#!/usr/bin/env python
#Makes simulated reads 
#Arg 1: genome size (ex. 2000000)
#Arg 2: Coverage (ex. 30 for 30x)
#Arg 3: Read size (ex. 10000)
#Arg 4: Percentage of heterozygosity, must be between 0-100 (ex. 1)
#Arg 5: Whether the het should be randomly distributed or not (1 == randomly distributed)
#Arg 6: Percent GC, must be between 0-100 (ex. 50) 
#Arg 7: Error rate, must be between 0-100 (ex. 1 for 1% error in the reads)
#Arg 8: Repeat rate, must be between 0-100 (ex. 2 for 2% repeats)
#Arg 9: Path to a template genome (must be single line fasta file with no "N" or "n"). Set as "0" if you don't want to provide a template genome.

import csv
import sys
import random
import statistics
import numpy as np
import scipy.stats as stats
from math import floor
from collections import defaultdict

genome_size = int(sys.argv[1]) #1000000
coverage = int(sys.argv[2]) #40
read_size = int(sys.argv[3]) #10000
het = float(sys.argv[4]) #0 #This needs to be a percentage (so 1 would mean 1%)
random_het = int(sys.argv[5])
rc = "rc" 
gc = float(sys.argv[6])
err_rate = float(sys.argv[7]) #This also needs to be a percentage (so "1" would mean 1%)
spread=1
repeats = float(sys.argv[8]) #This means there will be X % of the genome that is repetative 
copy_cat_genome = sys.argv[9]
out =     open("simulation.random.gs" + sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".fasta", "w")
out_mat = open("mat_genome_random.gs" + sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".txt"  , "w")
out_pat = open("pat_genome_random.gs" + sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".txt"  , "w")
out_het = open("het_loc_random.gs"    + sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".txt"  , "w")
out_max_kmers = open("max_kmers.gs"   + sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".txt"  , "w")
out_pairs    = open("true_pairs.gs"+    sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".tsv"  , "w") 
out_pairs_m1 = open("true_pairs_m1.gs"+ sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".tsv"  , "w") 
out_pairs_a1 = open("true_pairs_a1.gs"+ sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".tsv"  , "w")
out_pairs_m2 = open("true_pairs_m2.gs"+ sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".tsv"  , "w")
out_pairs_a2 = open("true_pairs_a2.gs"+ sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".tsv"  , "w")

print("Genome size: " + str(genome_size))
print("Coverage: " + str(coverage))
print("read size: " + str(read_size))
print("Random het: " + str(random_het))
print("Het: " + str(het))
print("GC: " + str(gc))
print("Error rate: " + str(err_rate))

def find(ch1, ch2, string1):
	pos = []
	for i in range(len(string1)):
		if ch1 == string1[i] or ch2 == string1[i]:
			pos.append(i)
	return pos

def copycat(copy_cat_genome, gs, gc):
	#Use a provided genome as a template 
	print("Using the template")
	copy_cat = ""
	with open(copy_cat_genome, encoding="utf8", errors='ignore') as f: 
		print("Opening template genome")
		Lines = f.readlines()
		for line in Lines:
			line = line.strip()
			if (line.count(">") == 0):
				copy_cat = copy_cat + line
	#Take the first X base pairs from the template to use 	
	mat_genome = copy_cat[:gs]

	#Option to change the amount of GC 
	#The main issue with this is it could break up repeats within the template genome. 
	if int(gc) != -1:
		print("changing gc") 
		#What is the current gc content? 
		gc_content = (mat_genome.count("C")+mat_genome.count("G"))/(len(mat_genome))

		if (gc_content*100) > gc: #The mat genome has too many gcs 
			while ((gc_content*100) > gc):
				num_bp_to_chng = int(((((gc_content*100) - gc)/100)*len(mat_genome))+1)
				nums = random.sample(find("G", "C", mat_genome), num_bp_to_chng)
				counter = 0 #Just for our progress report. 
				for i in nums:
					counter = counter + 1
					if (i%2) == 0:
						mat_genome = mat_genome[:(i)] + "A" + mat_genome[(i+1):]
					else:
						mat_genome = mat_genome[:(i)] + "T" + mat_genome[(i+1):]
					if counter%100000 == 0:
						gc_content = (mat_genome.count("C") + mat_genome.count("G"))/(len(mat_genome))
						print(str((gc_content*100)))
				gc_content = (mat_genome.count("C")+mat_genome.count("G"))/(len(mat_genome))

		elif (gc_content*100) < gc: #Here, our mat genome doesn't have enough gc! 		
			while ((gc_content*100) < gc):
				num_bp_to_chng = int((((gc - (gc_content*100))/100)*len(mat_genome))+1)
				nums = random.sample(find("A", "T", mat_genome), num_bp_to_chng)
				counter = 0 #Just for our progress report. 
				for i in nums:
					counter = counter + 1
					if (i%2) == 0:
						mat_genome = mat_genome[:(i)] + "G" + mat_genome[(i+1):]
					else:
						mat_genome = mat_genome[:(i)] + "C" + mat_genome[(i+1):]
					if counter%100000 == 0:
						gc_content = (mat_genome.count("C") + mat_genome.count("G"))/(len(mat_genome))
						print(str((gc_content*100)))
				gc_content = (mat_genome.count("C")+mat_genome.count("G"))/(len(mat_genome))

	print("Done with GC addition")
	return mat_genome

def add_errors(string, err_rate): #Will want to add errors to the reads. This will be on a per-read basis
	read_w_errors=list(string)
	#read_w_errors=string
	#nums = random.sample(range(0,len(string)), int((err_rate/100)*len(string)))
	nums = get_random_sample(len(string), err_rate/100)
        #nums.sort()
	read_err_loc = ""
	
	for i in nums:
		read_w_errors[i] = random_string_exclude(read_w_errors[i])
		#read_w_errors = read_w_errors[:(i)] + random_string_exclude(read_w_errors[i]) + read_w_errors[i+1:]
		read_err_loc = read_err_loc + str(i) + ","
	#return read_w_errors, read_err_loc
	return "".join(read_w_errors), read_err_loc

def find_het(het_is_here, read_length, starting_position, rc_or_c):
	#There totally must be an easier way to do this. But this is sufficient for now 
	het_loc_for_read = ""
	y = [x for x in het_is_here if starting_position <= x < (starting_position+read_length)]
	
	if rc_or_c == 1:
		z = [(read_length-1 - (x-starting_position)) for x in y]
	else:
		z = [x-starting_position for x in y]
	
	for x in z:
		het_loc_for_read = het_loc_for_read  + str(x) + ","
	
	return het_loc_for_read

def random_string(length):
	#Self explanatory
	string = ""
	for i in range(length):
		nuc = random.choice("ACTG")
		string += nuc
	return string

def weighted_string(length, gc):
	#Like random_string, but will randomly choose from a string that contains a specified amount
	#of Gs and Cs
	#ref = "C"*int(gc*100) + "G"*int(gc*100) + "A"*(100-int(gc*100)) + "T"*(100-int(gc*100))
	ref = "C"*int(gc) + "G"*int(gc) + "A"*int(100-gc) + "T"*int(100-gc)
	string = ""
	for i in range(length):
		nuc = random.choice(ref)
		string += nuc
	return string

def random_string_exclude(char):
	#When we want to pick a random letter (ACTG) but it can't be the same character as provided
	new_char = char
	while new_char == char:
		new_char = random.choice("ACTG")
	return new_char

def reverse_complement(string):
  """
  Assumes characters of string are upper case.
  """
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
  return "".join(complement[base] for base in reversed(string))

def prioritize_snp(i, kmer1, kmer2):
	window = -1
	while True:
		window += 1
		flipped = False
		A = kmer1[i-window:i+window+1]
		B = kmer2[i-window:i+window+1]
		if len(A) != 2*window+1:
			return
		if A > B:
			flipped = True
			A, B = B, A
		A_rc = reverse_complement(A)
		B_rc = reverse_complement(B)
		if (A < B < B_rc < A_rc) or (A < B_rc < B < A_rc):
			if flipped:
				return 1
			else:
				return 0
		if (B_rc < A < A_rc < B) or (B_rc < A_rc < A < B):
			if flipped:
				return 0
			else:
				return 1
		if (A == B_rc < B == A_rc):
			continue
		if (B_rc < A < B < A_rc):
			if flipped:
				return 1
			else:
				return 0
		if flipped:
			return 0
		else:
			return 1

#this function is for the genome level
def get_kmer_pair_locations(het_locations, k, genome_size):
	"""m1 for middle one away
	a1 for any one away
	m2 for middle two away
	a2 for any two away"""
	m1 = []
	a1 = []
	m2 = []
	a2 = []
	last_position_checked = -1
	set_het_locations = set(het_locations)
	for het_location in het_locations:
		#for each kmer overlapping the snp
		for i in range(max(0, het_location-k+1), min(genome_size-k+1, het_location+1)):
			#if we haven't already checked this position
			if i > last_position_checked:
				last_position_checked = i
				overlapping_snps = [position for position in range(i, i+k) if position in set_het_locations]
				#print(overlapping_snps)
				if len(overlapping_snps) == 1:
					a1.append(i)
					a2.append(i)
					overlapping_snp = overlapping_snps[0]
					if overlapping_snp == i+floor(k/2):
						m1.append(i)
						m2.append(i)
				if len(overlapping_snps) == 2:
					a2.append(i)
					overlapping_snp1, overlapping_snp2 = overlapping_snps
					if overlapping_snp1 == i+floor(k/2) or overlapping_snp2 == i+floor(k/2):
						m2.append(i)
					#for g in range(1, k):
					#	if (overlapping_snp1 == i+floor(k/2)-floor((g+1)/2)) and (overlapping_snp2 == overlapping_snp1+g):
					#		m2.append(i)
	return m1, a1, m2, a2

def get_middle_kmer_pair_locations(het_locations, k, genome_size):
    m = []
    for het_location in het_locations:
        #if the kmer doesn't go past the end of the genomes
        start = het_location-int((k-1)/2)
        end = start + (k-1)
        if start >= 0 and end <= genome_size - 1:
            m.append(start)
    return m

def get_kmer_pairs(kmer_pair_locations, mat_genome, pat_genome, k):
    kmer_pairs = []
    for kmer_pair_location in kmer_pair_locations:
        mat_kmer = mat_genome[kmer_pair_location:kmer_pair_location+k]
        pat_kmer = pat_genome[kmer_pair_location:kmer_pair_location+k]
        kmer_pairs.append((mat_kmer, pat_kmer))
    return kmer_pairs

def snp_is_editable(mode, overlapping_snps, i, k):
	if mode == "m1":
		if (len(overlapping_snps) == 1) and (overlapping_snps[0] == i+floor(k/2)):
			return True
		else:
			return False
	if mode == "a1":
		if len(overlapping_snps) == 1:
			return True
		else:
			return False
	if mode == "m2":
		#if ((len(overlapping_snps) == 1) and (overlapping_snps[0] == i+floor(k/2))) or ((len(overlapping_snps) == 2) and (any([(overlapping_snps[0]==i+floor(k/2)-floor((g+1)/2)) and (overlapping_snps[1]==overlapping_snps[0]+g) for g in range(1,k)]))):
		if ((len(overlapping_snps) == 1) and (overlapping_snps[0] == i+floor(k/2))) or ((len(overlapping_snps) == 2) and ((overlapping_snps[0] == i+floor(k/2)) or (overlapping_snps[1] == i+floor(k/2)))):
			return True
		else:
			return False
	if mode == "a2":
		if (len(overlapping_snps)==1) or (len(overlapping_snps)==2):
			return True
		else:
			return False

#this function is for the genome level
def get_editable_snp_genome_locations(het_locations, k, genome_size):
	"""m1 for middle one away
	a1 for any one away
	m2 for middle two away
	a2 for any two away"""
	modes = ["m1", "a1", "m2", "a2"]
	results = defaultdict(list)
	set_het_locations = set(het_locations)
	for mode in modes:
		for het_location in het_locations:
			#for each kmer overlapping the snp
			for i in range(max(0, het_location-k+1), min(genome_size-k+1, het_location+1)):
				overlapping_snps = [position for position in range(i, i+k) if position in set_het_locations]
				#print(overlapping_snps)
				if snp_is_editable(mode, overlapping_snps, i, k):
					results[mode].append(het_location)
					break
	m1 = results["m1"]
	a1 = results["a1"]
	m2 = results["m2"]
	a2 = results["a2"]
	return m1, a1, m2, a2

#this function is for the read level
def get_editable_snp_locations(read_start, read_het_locations, read_err_locations, k, read_size, genome_het_location_to_radius, read_sequence, rep, mode):
	results = []
	uneditable_by_error = []
	uneditable_by_density = []
	uneditable_by_ambiguity = []
	uneditable_by_pair = []
	#a1 = []
	#m2 = []
	#a2 = []
	for read_het_location in read_het_locations:
		#m1_flag = False
		#a1_flag = False
		#m2_flag = False
		#a2_flag = False
		snp_editable = False
		uneditable_by_error_flag = False
		uneditable_by_density_flag = False
		uneditable_by_ambiguity_flag = False
		uneditable_by_pair_flag = False
		for i in range(max(0, read_het_location-k+1), min(read_size-k+1, read_het_location+1)):
			overlapping_errors = [position for position in range(i, i+k) if position in read_err_locations]
			#If there are any overlapping errors, go to next window
			if overlapping_errors:
				uneditable_by_error_flag = True
				continue
			#Note: there will always be at least one overlapping_snp in overlapping_snps below
			overlapping_snps = [position for position in range(i, i+k) if position in read_het_locations]
			#If the snp density is uneditable under the given mode, go to next window
			if (mode == "m1") and ((len(overlapping_snps)!=1) or (overlapping_snps[0]!=i+floor(k/2))):
				uneditable_by_density_flag = True
				continue
			if (mode == "a1") and (len(overlapping_snps)!=1):
				uneditable_by_density_flag = True
				continue
			if (mode == "m2") and [[(len(overlapping_snps)!=1) or (overlapping_snps[0]!=i+floor(k/2))] and [(len(overlapping_snps)!=2) or all([(overlapping_snps[0]!=i+floor(k/2)-floor((g+1)/2)) or (overlapping_snps[1]!=overlapping_snps[0]+g) for g in range(1,k)])]]:
				uneditable_by_density_flag = True
				continue
			if (mode == "a2") and ((len(overlapping_snps)!=1) and (len(overlapping_snps)!=2)):
				uneditable_by_density_flag = True
				continue
			genome_het_location = read_start + read_het_location
			radius = genome_het_location_to_radius[genome_het_location]
			#If the snp is ambiguous, go to next window
			if (read_het_location-radius < i) or (read_het_location+radius > i+k-1):
				uneditable_by_ambiguity_flag = True
				continue
			#If the kmer is not in the kmer pair, go to next window
			if read_sequence[i:i+k] not in rep:
				uneditable_by_pair_flag = True
				continue
			#If we get this far, then the snp is editable
			snp_editable = True
			results.append(read_het_location)
			break
		if not snp_editable:
			if uneditable_by_error_flag:
				uneditable_by_error.append(read_het_location)
			if uneditable_by_density_flag:
				uneditable_by_density.append(read_het_location)
			if uneditable_by_ambiguity_flag:
				uneditable_by_ambiguity.append(read_het_location)
			if uneditable_by_pair_flag:
				uneditable_by_pair.append(read_het_location)
	return results, uneditable_by_error, uneditable_by_density, uneditable_by_ambiguity, uneditable_by_pair 

def get_random_sample(genome_size, het):
    population = range(genome_size)
    return [x for x in population if random.random() <= het]

def add_het(first_genome, het_level): 
	#Adds the het randomly  
	#nums = random.sample(range(0,len(first_genome)), int((het_level/100)*len(first_genome)))
	nums = get_random_sample(len(first_genome), het_level/100)
        #nums.sort()
	#with open('track1.bed', 'w') as f:
	#	for num in nums:
	#		f.write(f"chr20\t{num}\t{num+1}\n")
	mat_nums = []
	pat_nums = []
	new_genome = list(first_genome)
	x = 0
	for i in nums:
		x = x + 1
		if x%10000==0:
			print(str(x/len(nums)*100) + "%")
		#####new_genome = new_genome + first_genome[i]
		new_genome[i] = random_string_exclude(new_genome[i]) 
		#new_genome = new_genome[:i] + random_string_exclude(new_genome[i]) + new_genome[i+1:]
		#kmer1 = first_genome[max(0,i-10):i+11]
		kmer1 = "".join(first_genome[max(0,i-10):i+11])
		#kmer2 = new_genome[max(0,i-10):i+11]
		kmer2 = "".join(new_genome[max(0,i-10):i+11])
		priority = prioritize_snp(i-max(0, i-10), kmer1, kmer2)
		if priority == 0:
			mat_nums.append(i)
		elif priority == 1:
			pat_nums.append(i)
		else:
			print("het change at position " + str(i) + " is ambiguous.")
			f = open('ambiguous_positions.txt', 'a')
			f.write(str(i)+'\n')
			f.close()
	print("added het")
	#return new_genome, mat_nums, pat_nums
	return "".join(new_genome), mat_nums, pat_nums
	
def add_het_evenly(first_genome, het_level):
	#This time we aren't going to distribute the het randomly. This time it will be even. 
	new_genome = first_genome
	nums = list(range(0,(len(first_genome))))[(int(100 / het_level)-1)::int(100 / het_level)]
	for i in nums:
		new_genome = new_genome[:i] + random_string_exclude(new_genome[i]) + new_genome[i+1:]
	return new_genome, nums

def add_repeats(sequence, repeat_rate, spread):
	#There must be a better name than "spread", but basically this is how many groups 
	#of repeats to add 
	spread = int(len(sequence)/1000) 
	if spread == 0:
		spread = 1
	print(spread)
	repeat_seq = sequence
	#First, calculate how many bases will be effected.
	effected_bases = len(sequence) * (repeat_rate / 100) 
	#Next calculate the size of the repeat "chunks" 
	chunk_size = effected_bases / spread
	#Then determine where in the genome we want to place the repeated regions
	nums = random.sample(range(0,len(sequence)), spread) 
	#Now, let's go through and replace them
	for i in nums:
		mer_len=1 #Only adding homopolymers for now 
		# Make the mer
		mer = repeat_seq[i:(i+mer_len)]
		# Make the entire repeat seq
		add_this = ""	
		while len(add_this) <= chunk_size :
			add_this = add_this + mer
		add_this = add_this[:int(chunk_size)]
		# Add the repeats to the genome 
		repeat_seq = repeat_seq[:i] + add_this + repeat_seq[(i+len(add_this)):]
	print("Added repeats")
	return repeat_seq[:len(sequence)] 

#####################################
#Main chunk
if sys.argv[9] == "0": #We aren't given a fasta file to use as a template
	mat_genome = weighted_string(genome_size,gc)
else: #We are given a copycat genome! 
	print("Using this as a template: "+copy_cat_genome)
	mat_genome = copycat(copy_cat_genome, genome_size, gc)

#Add repeats
if repeats != 0.0:
	mat_genome = add_repeats(mat_genome, repeats, spread)
mat_genome = mat_genome[:genome_size]

#Add het
if (het != 0.0) and (random_het == 1): #We want to randomly dist the het
	pat_genome, mat_het_is_here, pat_het_is_here = add_het(mat_genome, het)
	het_is_here = sorted(mat_het_is_here+pat_het_is_here)
	#m1, a1, m2, a2 = get_editable_snp_genome_locations(het_is_here, 21, genome_size) #assume k = 21
	#with open('track2.m1.bed', 'w') as f:
	#	for num in m1:
	#		f.write(f"chr20\t{num}\t{num+1}\n")
	#with open('track2.a1.bed', 'w') as f:
	#	for num in a1:
	#		f.write(f"chr20\t{num}\t{num+1}\n")
	#with open('track2.m2.bed', 'w') as f:
	#	for num in m2:
	#		f.write(f"chr20\t{num}\t{num+1}\n")
	#with open('track2.a2.bed', 'w') as f:
	#	for num in a2:
	#		f.write(f"chr20\t{num}\t{num+1}\n")
	m1, a1, m2, a2 = get_kmer_pair_locations(het_is_here, 21, genome_size) #assume k = 21
	m = get_middle_kmer_pair_locations(het_is_here, 21, genome_size)
	m1_pairs = get_kmer_pairs(m1, mat_genome, pat_genome, 21)
	a1_pairs = get_kmer_pairs(a1, mat_genome, pat_genome, 21)
	m2_pairs = get_kmer_pairs(m2, mat_genome, pat_genome, 21)
	a2_pairs = get_kmer_pairs(a2, mat_genome, pat_genome, 21)
	m_pairs  = get_kmer_pairs(m,  mat_genome, pat_genome, 21)
	out_max_kmers.write(str(len(m1))+'\n')
	out_max_kmers.write(str(len(a1))+'\n')
	out_max_kmers.write(str(len(m2))+'\n')
	out_max_kmers.write(str(len(a2))+'\n')
	#with open('trackp.m1.bed', 'w') as f:
	#	for num in m1:
	#		f.write(f"chr20\t{num}\t{num+1}\n")
	#with open('trackp.a1.bed', 'w') as f:
	#	for num in a1:
	#		f.write(f"chr20\t{num}\t{num+1}\n")
	#with open('trackp.m2.bed', 'w') as f:
	#	for num in m2:
	#		f.write(f"chr20\t{num}\t{num+1}\n")
	#with open('trackp.a2.bed', 'w') as f:
	#	for num in a2:
	#		f.write(f"chr20\t{num}\t{num+1}\n")
	for m1_pair in m1_pairs:
		kmer1, kmer2 = m1_pair
		out_pairs_m1.write(f"{kmer1}\t{kmer2}\n")
	for a1_pair in a1_pairs:
		kmer1, kmer2 = a1_pair
		out_pairs_a1.write(f"{kmer1}\t{kmer2}\n")
	for m2_pair in m2_pairs:
		kmer1, kmer2 = m2_pair
		out_pairs_m2.write(f"{kmer1}\t{kmer2}\n")
	for a2_pair in a2_pairs:
		kmer1, kmer2 = a2_pair
		out_pairs_a2.write(f"{kmer1}\t{kmer2}\n")
	for m_pair in m_pairs:
		kmer1, kmer2 = m_pair
		out_pairs.write(f"{kmer1}\t{kmer2}\n")
	#print(mat_het_is_here)
	#print(pat_het_is_here)
elif (het != 0.0) and (random_het == 0): #We want to evenly dist the het
	pat_genome, het_is_here = add_het_evenly(mat_genome, het)
else: 
	pat_genome = mat_genome
	het_is_here = []
	mat_het_is_here = []
	pat_het_is_here = []

#We want to have a copy of the mat and pat actual genome (as a reference)
out_mat.write(mat_genome)
out_pat.write(pat_genome)

#Calculate the number of reads that we'll need to make 
number_reads = (genome_size * coverage) / read_size

l = 0
for i in range(0,int(number_reads)): #This is where we make the reads. 
	#Pick a random place to start 
	if i%1000 == 0:
		print(str(i/int(number_reads)*100) + "%")
	random_number = random.randint(0,genome_size)
	#read_length = 
	read_het_locations = "" #Keeping track of where we add the het. 
	if (i%2) == 0: #Use mat 
		this_read = mat_genome[random_number:(random_number + read_size)]
		which_genome="mat"
	else: # Use pat 
		this_read = pat_genome[random_number:(random_number + read_size)]
		which_genome = "pat"
	
	#Decide if we will use the reverse complement or not 
	rc_or_c = random.randint(0,1)
	rc="original"
	if rc_or_c == 1:
		rc = "reverse_complement"
		this_read = reverse_complement(this_read)
	
	#Find where the het was added. There is totally an easier way, but for now this works well enough.
	if which_genome == "mat":
		read_het_locations = find_het(mat_het_is_here, len(this_read), random_number, rc_or_c)
	else:
		read_het_locations = find_het(pat_het_is_here, len(this_read), random_number, rc_or_c)

	#Finally, add error to each read. 
	this_read, read_err_loc = add_errors(this_read, err_rate) 
	if len(this_read) != read_size:
		print("Found one!")
		print("read_"+ str(i))
		print("random_number: " + str(random_number))
		print("read length: " + str(len(this_read)))
		print("This read: " + this_read)
	#write this in a fasta file
	#The name for the read will include the read #, the forward starting location of the read (or forward ending position of the reverse complement read, 
	#whether it is from the mat or pat,
	#Whether it is the original or reverse complement strand, the het locations, the error locations
	#and the read length
	out.write(">read_"+str(i) + "|"+ str(random_number) + "|"  + which_genome + "|" + rc +  "|" + str(read_het_locations) + "|" + str(read_err_loc) + "|" + str(len(this_read)) + "\n")
	out.write(this_read + "\n")

for i in mat_het_is_here:
	out_het.write(str(i) + ",")
out_het.write("\n")
for i in pat_het_is_here:
	out_het.write(str(i) + ",")
out_het.write("\n")

out_het.close()
out.close()
out_mat.close()
out_pat.close()
out_max_kmers.close()
out_pairs_m1.close()
out_pairs_a1.close()
out_pairs_m2.close()
out_pairs_a2.close()
