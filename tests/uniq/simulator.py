#!/usr/bin/env python

import argparse
import math
import os
import random

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outdir", help="Path for the output directory", default=".", dest="outdir")
parser.add_argument("-g", "--genome_length", help="Genome length (bp) for maternal genome", default=1000000, type=int, dest="mat_genome_length")
parser.add_argument("-c", "--coverage", help="Read coverage for sequencing data", default=30, type=int, dest="coverage")
parser.add_argument("-r", "--read_length", help="Read length (bp)", default=10000, type=int, dest="read_length")
parser.add_argument("--het_rate", help="Heterozygosity rate percent (e.g. 1 for 1 percent het between haplotypes)", default=1.0, type=float, dest="het_rate")
parser.add_argument("-e", "--err_rate", help="Error rate percent (e.g. 1 for 1 percent err in the reads)", default=1.0, type=float, dest="err_rate")
parser.add_argument("-i", "--indel_length", help="Optional parameter to set the indel length allowed in heterozygosity and error", default=0, type=int, dest="indel_length")
parser.add_argument("-t", "--template_genome", help="Optional parameter for path to a template genome (must be a single line fasta file with no 'N' or 'n' and all caps ('A', 'C', 'G', 'T')).", default="", dest="template_genome")
parser.add_argument("-k", "--kmer_length", help="Kmer length (necessary to calculate expected heterozygous and error blocks)", default=21, type=int, dest="k")

args = parser.parse_args()
outdir = args.outdir
mat_genome_length = args.mat_genome_length
coverage = args.coverage
read_length = args.read_length
het_rate = args.het_rate
err_rate = args.err_rate
indel_length = args.indel_length
template_genome = args.template_genome
k = args.k

if template_genome:
	template_bool = 1
else:
	template_bool = 0

#if output folder doesn't exist create it
if not os.path.exists(outdir):
	os.makedirs(outdir)

out_reads     = open(f'{outdir}/simulatedreads_template{template_bool}_g{mat_genome_length}_cov{coverage}_rl{read_length}_err{err_rate}_het{het_rate}_indel{indel_length}.fasta', 'w')
#these are the reads before adding any errors
out_purereads = open(f'{outdir}/purereads_template{template_bool}_g{mat_genome_length}_cov{coverage}_rl{read_length}_err{err_rate}_het{het_rate}_indel{indel_length}.fasta', 'w')
out_mat       = open(f'{outdir}/matgenome_template{template_bool}_g{mat_genome_length}_cov{coverage}_rl{read_length}_err{err_rate}_het{het_rate}_indel{indel_length}.fasta', 'w')
out_pat       = open(f'{outdir}/patgenome_template{template_bool}_g{mat_genome_length}_cov{coverage}_rl{read_length}_err{err_rate}_het{het_rate}_indel{indel_length}.fasta', 'w')
out_het       = open(f'{outdir}/hetloc_template{template_bool}_g{mat_genome_length}_cov{coverage}_rl{read_length}_err{err_rate}_het{het_rate}_indel{indel_length}.txt', 'w')
out_genomehetblocks = open(f'{outdir}/genomehetblocks_template{template_bool}_g{mat_genome_length}_cov{coverage}_rl{read_length}_err{err_rate}_het{het_rate}_indel{indel_length}.fasta', 'w')
out_purereadhetblocks = open(f'{outdir}/purereadhetblocks_template{template_bool}_g{mat_genome_length}_cov{coverage}_rl{read_length}_err{err_rate}_het{het_rate}_indel{indel_length}.fasta', 'w')
out_readerrblocks = open(f'{outdir}/readerrblocks_template{template_bool}_g{mat_genome_length}_cov{coverage}_rl{read_length}_err{err_rate}_het{het_rate}_indel{indel_length}.fasta', 'w')

mutations = {"A":"CTG", "C":"ATG", "G":"ACT", "T":"ACG"}

def reverse_complement(string):
  """
  Assumes characters of string are upper case.
  """
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
  return "".join(complement[base] for base in reversed(string))

def add_het(mat_genome, het_rate, indel_length):
	#currently assuming that the indel length is constant for all indels
	pat_genome = ""
	mat_to_pat = []
	pat_to_mat = []
	#het_locations are a list of three sets
	#the first set has the positions of substititions
	#the second set has the positions of insertions
	#the third set has the positions of deletions
	#the mat locations are with respect to mat_genome (i.e. mat -> pat)
	#the pat locations are with respect to pat_genome (i.e. pat -> mat)
	mat_het_locations = [set(), set(), set()]
	pat_het_locations = [set(), set(), set()]
	#i will keep track of the current position in mat_genome (not yet written to)
	#j will keep track of the current position in pat_genome (not yet written to)
	i = 0
	j = 0
	#iterate over all positions of the mat_genome to determine whether there is het
	while (i < len(mat_genome)):
		#If this position is a het
		if random.random() <= het_rate/100:
			#what kind of het is it?
			#if indels are allowed, then het can be sub, ins, or del
			if indel_length > 0:
				#0 = substition, 1 = insertion, 2 = deletion
				het_type = random.randint(0,2)
			#if indels are not allowed, then het is sub
			else:
				het_type = 0
			#substition (occurs at position)
			if het_type == 0:
				ref = mat_genome[i]
				alt = random.choice(mutations[ref])
				pat_genome += alt
				mat_het_locations[0].add(i)
				pat_het_locations[0].add(j)
				mat_to_pat.append(j)
				pat_to_mat.append(i)
				i += 1
				j += 1
			#insertion (occurs after position)
			if het_type == 1:
				ref = mat_genome[i]
				insertion = "".join(random.choices("ACTG", k=indel_length))
				pat_genome += ref + insertion
				mat_het_locations[1].add(i)
				pat_het_locations[2].add(j+1)
				mat_to_pat.append(j)
				pat_to_mat.append(i)
				pat_to_mat += [-1]*indel_length
				i += 1
				j += 1 + indel_length
			#deletion (occurs at position)
			if het_type == 2: #TO FIGURE OUT, WHAT HAPPENS FOR DELETION AT END OF GENOME, THE LENGTH MAY NOT MATCH INDEL_LENGTH
				mat_het_locations[2].add(i)
				pat_het_locations[1].add(j-1)
				mat_to_pat += [-1]*indel_length
				i += indel_length
				j += 0 
		#else this position is not a het
		else:
			pat_genome += mat_genome[i]
			mat_to_pat.append(j)
			pat_to_mat.append(i)
			i += 1
			j += 1
	return pat_genome, mat_het_locations, pat_het_locations, mat_to_pat, pat_to_mat

def add_err(pureread, err_rate, indel_length):
	#currently assuming that the indel length is constant for all indels
	read = ""
	pureread_to_read = []
	read_to_pureread = []
	read_start_diff = 0
	check_read_start_diff = True
	#err_locations are a list of three sets
	#the first set has the positions of substititions
	#the second set has the positions of insertions
	#the third set has the positions of deletions
	#the pureread locations are with respect to pureread (i.e. pureread -> read)
	#the read locations are with respect to read (i.e. read -> pureread)
	pureread_err_locations = [set(), set(), set()]
	read_err_locations = [set(), set(), set()]
	#i will keep track of the current position in pureread (not yet written to)
	#j will keep track of the current position in read (not yet written to)
	i = 0
	j = 0
	#iterate over all positions of the pureread to determine whether there is err
	while (i < len(pureread)):
		#If this position is an err
		if random.random() <= err_rate/100:
			#what kind of err is it?
			#if indels are allowed, then err can be sub, ins, or del
			if indel_length > 0:
				#0 = substition, 1 = insertion, 2 = deletion
				err_type = random.randint(0,2)
			#if indels are not allowed, then err is sub
			else:
				err_type = 0
			#substition (occurs at position)
			if err_type == 0:
				ref = pureread[i]
				alt = random.choice(mutations[ref])
				read += alt
				pureread_err_locations[0].add(i)
				read_err_locations[0].add(j)
				pureread_to_read.append(j)
				read_to_pureread.append(i)
				i += 1
				j += 1
				check_read_start_diff = False
			#insertion (occurs after position)
			if err_type == 1:
				ref = pureread[i]
				insertion = "".join(random.choices("ACTG", k=indel_length))
				read += ref + insertion
				pureread_err_locations[1].add(i)
				read_err_locations[2].add(j+1)
				pureread_to_read.append(j)
				read_to_pureread.append(i)
				read_to_pureread += [-1]*indel_length
				i += 1
				j += 1 + indel_length
				check_read_start_diff = False
			#deletion (occurs at position)
			if err_type == 2: #TO FIGURE OUT, WHAT HAPPENS FOR DELETION AT END OF READ, THE LENGTH MAY NOT MATCH INDEL_LENGTH
				pureread_err_locations[2].add(i)
				read_err_locations[1].add(j-1)
				pureread_to_read += [-1]*indel_length
				i += indel_length
				j += 0
				if check_read_start_diff:
					read_start_diff += indel_length
		#else this position is not an err
		else:
			read += pureread[i]
			pureread_to_read.append(j)
			read_to_pureread.append(i)
			i += 1
			j += 1
			check_read_start_diff = False
	return read, pureread_err_locations, read_err_locations, pureread_to_read, read_to_pureread, read_start_diff

#####################################
#Main chunk

#If we are given a template genome
if template_bool:
	with open(template_genome, "r") as input_template_file:
		for line in input_template_file:
			if not line.startswith(">"):
				mat_genome = line.strip()
	mat_genome_length = len(mat_genome) #this overwrites any parameter set by -g, may need to change later
#else we are not given a template genome
else:
	mat_genome = "".join(random.choices("ACTG", k=mat_genome_length))
mat_rev_genome = reverse_complement(mat_genome)
out_mat.write(f">matgenome_template{template_bool}_g{mat_genome_length}_cov{coverage}_rl{read_length}_err{err_rate}_het{het_rate}_indel{indel_length}\n")
out_mat.write(mat_genome+"\n")
out_mat.close()

#Add het to pat_genome where each nucleotide
#independently has a het_rate chance of being a het
#if indel_length > 0 then each het site
#is 1/3 substitution, 1/3 insertion, 1/3 deletion
if (het_rate != 0.0):
	pat_genome, mat_het_locations, pat_het_locations, mat_to_pat, pat_to_mat = add_het(mat_genome, het_rate, indel_length)
else: 
	pat_genome = mat_genome
	mat_het_locations = [set(), set(), set()]
	pat_het_locations = [set(), set(), set()]
	mat_to_pat = list(range(mat_genome_length))
	pat_to_mat = list(range(mat_genome_length))
pat_genome_length = len(pat_genome)
mat_rev_het_locations = []
pat_rev_het_locations = []
mat_rev_het_locations.append({mat_genome_length - 1 - x for x in mat_het_locations[0]})
mat_rev_het_locations.append({mat_genome_length - 1 - (x+1) for x in mat_het_locations[1]})
mat_rev_het_locations.append({mat_genome_length - 1 - x for x in mat_het_locations[2]})
pat_rev_het_locations.append({pat_genome_length - 1 - x for x in pat_het_locations[0]})
pat_rev_het_locations.append({pat_genome_length - 1 - (x+1) for x in pat_het_locations[1]})
pat_rev_het_locations.append({pat_genome_length - 1 - x for x in pat_het_locations[2]})
pat_rev_genome = reverse_complement(pat_genome)
out_pat.write(f">patgenome_template{template_bool}_g{mat_genome_length}_cov{coverage}_rl{read_length}_err{err_rate}_het{het_rate}_indel{indel_length}\n")
out_pat.write(pat_genome+"\n")
out_pat.close()

#hetloc will have 12 lines
#maternal forward (sub, ins, del)
#paternal forward (sub, ins, del)
#maternal reverse (sub, ins, del)
#paternal reverse (sub, ins, del)
for het_locations in mat_het_locations:
	out_het.write(",".join(str(x) for x in sorted(het_locations)) + "\n")
for het_locations in pat_het_locations:
	out_het.write(",".join(str(x) for x in sorted(het_locations)) + "\n")
for het_locations in mat_rev_het_locations:
	out_het.write(",".join(str(x) for x in sorted(het_locations)) + "\n")
for het_locations in pat_rev_het_locations:
	out_het.write(",".join(str(x) for x in sorted(het_locations)) + "\n")
out_het.close()

#Calculate the genome het blocks
#(these will be approximate since by chance may insert the same next base in the sequence)
#(or they may delete to a location with the same next base in the sequence)
#traverse mat_to_pat to find any consecutive sequences (with no substitutions) of length >= k
previous_type = -1
mat_hetblock_locations = []
mat_rev_hetblock_locations = []
pat_hetblock_locations = []
pat_rev_hetblock_locations = []
for i in range(mat_genome_length):
	#if this continues the previous consecutive sequence with no substitutions
	if (i > 0 and mat_to_pat[i] > 0 and (i-1 not in mat_het_locations[0]) and (i not in mat_het_locations[0]) and mat_to_pat[i] == mat_to_pat[i-1] + 1):
		count += 1
		if count >= k:
			#we have left the beginning of the genome
			if previous_type == -1:
				pass
			#het to hom, we have left the het block
			if previous_type == 0:
				mat_lastnonhom_idx = i-k
				pat_lastnonhom_idx = mat_to_pat[i-k+1]-1
				mat_rev_firstnonhom_idx = mat_genome_length - 1 - (mat_lastnonhom_idx + k - 1)
				mat_rev_lastnonhom_idx = mat_genome_length - 1 - (mat_firstnonhom_idx + k - 1)
				pat_rev_firstnonhom_idx = pat_genome_length - 1 - (pat_lastnonhom_idx + k - 1)
				pat_rev_lastnonhom_idx = pat_genome_length - 1 - (pat_firstnonhom_idx + k - 1)
				mat_hetblock_locations.append((mat_firstnonhom_idx, mat_lastnonhom_idx))
				mat_rev_hetblock_locations.append((mat_rev_firstnonhom_idx, mat_rev_lastnonhom_idx))
				pat_hetblock_locations.append((pat_firstnonhom_idx, pat_lastnonhom_idx))
				pat_rev_hetblock_locations.append((pat_rev_firstnonhom_idx, pat_rev_lastnonhom_idx))
				#write het block to out_genomehetblocks file
				mat_header = f">mat_forward_firstnonhomkmer{mat_firstnonhom_idx}_lastnonhomkmer{mat_lastnonhom_idx}"
				mat_hetblock = mat_genome[mat_firstnonhom_idx-1:mat_firstnonhom_idx-1+k] + " " + mat_genome[mat_firstnonhom_idx-1+k:mat_lastnonhom_idx+1] + " " + mat_genome[mat_lastnonhom_idx+1:mat_lastnonhom_idx+k+1]
				pat_header = f">pat_forward_firstnonhomkmer{pat_firstnonhom_idx}_lastnonhomkmer{pat_lastnonhom_idx}"
				pat_hetblock = pat_genome[pat_firstnonhom_idx-1:pat_firstnonhom_idx-1+k] + " " + pat_genome[pat_firstnonhom_idx-1+k:pat_lastnonhom_idx+1] + " " + pat_genome[pat_lastnonhom_idx+1:pat_lastnonhom_idx+k+1]
				mat_rev_header = f">mat_reverse_firstnonhomkmer{mat_rev_firstnonhom_idx}_lastnonhomkmer{mat_rev_lastnonhom_idx}"
				mat_rev_hetblock = mat_rev_genome[mat_rev_firstnonhom_idx-1:mat_rev_firstnonhom_idx-1+k] + " " + mat_rev_genome[mat_rev_firstnonhom_idx-1+k:mat_rev_lastnonhom_idx+1] + " " + mat_rev_genome[mat_rev_lastnonhom_idx+1:mat_rev_lastnonhom_idx+k+1]
				pat_rev_header = f">pat_reverse_firstnonhomkmer{pat_rev_firstnonhom_idx}_lastnonhomkmer{pat_rev_lastnonhom_idx}"
				pat_rev_hetblock = pat_rev_genome[pat_rev_firstnonhom_idx-1:pat_rev_firstnonhom_idx-1+k] + " " + pat_rev_genome[pat_rev_firstnonhom_idx-1+k:pat_rev_lastnonhom_idx+1] + " " + pat_rev_genome[pat_rev_lastnonhom_idx+1:pat_rev_lastnonhom_idx+k+1]
				out_genomehetblocks.write(mat_header + "\n")
				out_genomehetblocks.write(mat_hetblock + "\n")
				out_genomehetblocks.write(pat_header + "\n")
				out_genomehetblocks.write(pat_hetblock + "\n")
				out_genomehetblocks.write(mat_rev_header + "\n")
				out_genomehetblocks.write(mat_rev_hetblock + "\n")
				out_genomehetblocks.write(pat_rev_header + "\n")
				out_genomehetblocks.write(pat_rev_hetblock + "\n")
			#hom to hom, continuing hom block
			if previous_type == 1:
				pass
			previous_type = 1 #homozygous
		#if this is the start of the genome, then previous_type remains -1
		if previous_type == -1:
			continue
	#if this does not continue the previous consecutive sequence
	else:
		count = 1
		#if this is the start of the genome, then previous_type remains -1
		if previous_type == -1:
			continue
		#het to het, continuing het block
		if previous_type == 0:
			pass
		#hom to het, we have entered the het block
		if previous_type == 1:
			mat_firstnonhom_idx = i-k+1
			pat_firstnonhom_idx = mat_to_pat[i-k]+1
		previous_type = 0 #heterozygous
#We are at the end of the genome, but are we in the middle of a heterozygous block?
#If so, then we need to output the final het block NEED TO WRITE
#if previous_type == 0:
#	mat_lastnonhom_idx = mat_genome_length - 1
#	pat_lastnonhom_idx = pat_genome_length - 1
#	mat_rev_firstnonhom_idx = mat_genome_length - 1 - (mat_lastnonhom_idx + k - 1) #1-k, so we need to change this
#	mat_rev_lastnonhom_idx = mat_genome_length - 1 - (mat_firstnonhom_idx + k - 1)
#	pat_rev_firstnonhom_idx = pat_genome_length - 1 - (pat_lastnonhom_idx + k - 1)
#	pat_rev_lastnonhom_idx = pat_genome_length - 1 - (pat_firstnonhom_idx + k - 1)
#	mat_hetblock_locations.append((mat_firstnonhom_idx, mat_lastnonhom_idx))
#	mat_rev_hetblock_locations.append((mat_rev_firstnonhom_idx, mat_rev_lastnonhom_idx))
#	pat_hetblock_locations.append((pat_firstnonhom_idx, pat_lastnonhom_idx))
#	pat_rev_hetblock_locations.append((pat_rev_firstnonhom_idx, pat_rev_lastnonhom_idx))
	
out_genomehetblocks.close()
mat_rev_hetblock_locations.reverse()
pat_rev_hetblock_locations.reverse()

#Calculate the number of reads that we'll need to make
number_reads = math.ceil((mat_genome_length * coverage) / read_length)

#Simulate the reads
for i in range(number_reads):
	#Decide which genome (mat or pat) and which strand (forward or reverse) to simulate from
	genome_bool = random.randint(0,1)
	strand_bool = random.randint(0,1)
	#Use mat_genome
	if genome_bool == 0:
		genome_text = "mat"
		#Pick a random starting location for the read
		#randint(a,b): return a random integer N such that a <= N <= b
		pureread_start = random.randint(0,mat_genome_length-1)
		#forward strand
		if strand_bool == 0:
			strand_text = "forward"
			pureread = mat_genome[pureread_start:(pureread_start + read_length)]
			strand_het_locations = mat_het_locations
			strand_hetblock_locations = mat_hetblock_locations
		#reverse strand
		else:
			strand_text = "reverse"
			pureread = mat_rev_genome[pureread_start:(pureread_start + read_length)]
			strand_het_locations = mat_rev_het_locations
			strand_hetblock_locations = mat_rev_hetblock_locations
	#Use pat_genome
	else:
		genome_text = "pat"
		pureread_start = random.randint(0,pat_genome_length-1)
		#forward strand
		if strand_bool == 0:
			strand_text = "forward"
			pureread = pat_genome[pureread_start:(pureread_start + read_length)]
			strand_het_locations = pat_het_locations
			strand_hetblock_locations = pat_hetblock_locations
		#reverse strand
		else:
			strand_text = "reverse"
			pureread = pat_rev_genome[pureread_start:(pureread_start + read_length)]
			strand_het_locations = pat_rev_het_locations
			strand_hetblock_locations = pat_rev_hetblock_locations
	
	#The length of the pure read may be less than intended if read start occurs near end of strand
	pureread_actual_length = len(pureread)
	pureread_end = pureread_start + pureread_actual_length - 1 #inclusive end
	
	#Check which hetblocks overlap the pureread
	for first, last in strand_hetblock_locations:
		#We check whether each anchor (+ 1 nucleotide) overlaps the pureread (strand coordinates)
		rightanchor_overlaps = ((pureread_start <= last) & ((last+k) <= pureread_end))
		leftanchor_overlaps = ((pureread_start <= (first-1)) & ((first-1+k) <= pureread_end))
		#Convert to pureread coordinates
		first -= pureread_start
		last -= pureread_start
		if leftanchor_overlaps and rightanchor_overlaps:
			pureread_hetblock = pureread[first-1:first-1+k] + " " + pureread[first-1+k:last+1] + " " + pureread[last+1:last+1+k]
		if rightanchor_overlaps and not leftanchor_overlaps:
			if 0 < first-1+k:
				pureread_hetblock = pureread[:first-1+k] + " " + pureread[first-1+k:last+1] + " " + pureread[last+1:last+1+k]
			else:
				pureread_hetblock = "" + " " + pureread[:last+1] + " " + pureread[last+1:last+1+k]
		if leftanchor_overlaps and not rightanchor_overlaps:
			if last < pureread_actual_length - 1:
				pureread_hetblock = pureread[first-1:first-1+k] + " " + pureread[first-1+k:last+1] + " " + pureread[last+1:]
			else:
				pureread_hetblock = pureread[first-1:first-1+k] + " " + pureread[first-1+k:] + " " + ""
		#write to out_purereadhetblocks if there was an overlap (i.e. one of the three cases above)
		if leftanchor_overlaps or rightanchor_overlaps:
			pureread_hetblock_header = f">pureread{i}_firstnonhomkmer{first}_lastnonhomkmer{last}"
			out_purereadhetblocks.write(pureread_hetblock_header + "\n")
			out_purereadhetblocks.write(pureread_hetblock + "\n")
	
	#get het locations on pure read (with respect to pure read index)
	pureread_het_locations = []
	#initially we get the indices with respect to the strand het locations
	pureread_het_locations.append(strand_het_locations[0].intersection(range(pureread_start, pureread_start + pureread_actual_length)))
	pureread_het_locations.append(strand_het_locations[1].intersection(range(pureread_start, pureread_start + pureread_actual_length)))
	pureread_het_locations.append(strand_het_locations[2].intersection(range(pureread_start, pureread_start + pureread_actual_length)))
	#now we correct to get indices with respect to the pure read
	pureread_het_locations[0] = {x-pureread_start for x in pureread_het_locations[0]}
	pureread_het_locations[1] = {x-pureread_start for x in pureread_het_locations[1]}
	pureread_het_locations[2] = {x-pureread_start for x in pureread_het_locations[2]}
	pureread_het_locations_text = "[" + ",".join(str(x) for x in sorted(pureread_het_locations[0])) + "][" + ",".join(str(x) for x in sorted(pureread_het_locations[1])) + "][" + ",".join(str(x) for x in sorted(pureread_het_locations[2]))+"]"
	
	#Now we add errors to the read, in the same way that we added het to the genomes
	if (err_rate != 0.0):
		read, pureread_err_locations, read_err_locations, pureread_to_read, read_to_pureread, read_start_diff = add_err(pureread, err_rate, indel_length)
	else:
		read = pureread
		pureread_err_locations = [set(), set(), set()]
		read_err_locations = [set(), set(), set()]
		pureread_to_read = list(range(pureread_actual_length))
		read_to_pureread = list(range(pureread_actual_length))
		read_start_diff = 0
	
	#The length of the read may change based on any insertions/deletions added
	read_actual_length = len(read)
	
	#calculate the error blocks (read -> pureread)
	#(these will be approximate since by chance may insert the same next base in the sequence)
	#(or they may delete to a location with the same next base in the sequence)
	#traverse read_to_pureread to find any consecutive sequences (with no substitutions) of length >= k
	previous_type = -1
	for l in range(read_actual_length):
		#if this continues the previous consecutive sequence with no substitutions
		if (l > 0 and read_to_pureread[l] > 0 and (l-1 not in read_err_locations[0]) and (l not in read_err_locations[0]) and read_to_pureread[l] == read_to_pureread[l-1] + 1):
			count += 1
			if count >= k:
				#we have left the beginning of the read
				if previous_type == -1:
					pass
				#err to nonerr, we have left the err block
				if previous_type == 0:
					read_lasterr_idx = l-k
					pureread_lasterr_idx = read_to_pureread[l-k+1]-1
					#write err block to out_readerrblocks
					read_header = f">read{i}_firsterrorkmer{read_firsterr_idx}_lasterrorkmer{read_lasterr_idx}"
					read_hetblock = read[read_firsterr_idx-1:read_firsterr_idx-1+k] + " " + read[read_firsterr_idx-1+k:read_lasterr_idx+1] + " " + read[read_lasterr_idx+1:read_lasterr_idx+k+1]
					pureread_header = f">pureread{i}_firsterrorkmer{pureread_firsterr_idx}_lasterrorkmer{pureread_lasterr_idx}"
					pureread_hetblock = pureread[pureread_firsterr_idx-1:pureread_firsterr_idx-1+k] + " " + pureread[pureread_firsterr_idx-1+k:pureread_lasterr_idx+1] + " " + pureread[pureread_lasterr_idx+1:pureread_lasterr_idx+k+1]
					out_readerrblocks.write(read_header + "\n")
					out_readerrblocks.write(read_hetblock + "\n")
					out_readerrblocks.write(pureread_header + "\n")
					out_readerrblocks.write(pureread_hetblock + "\n")
				#nonerr to nonerr, continuing nonerr block
				if previous_type == 1:
					pass
				previous_type = 1 #non error
			#if this is the start of the read, then previous_type remains -1
			if previous_type == -1:
				continue
		#if this does not continue the previous consecutive sequence
		else:
			count = 1
			#if this is the start of the read, then previous_type remains -1
			if previous_type == -1:
				continue
			#err to err, continuing err block
			if previous_type == 0:
				pass
			#nonerr to err, we have entered the err block
			if previous_type == 1:
				read_firsterr_idx = l-k+1
				pureread_firsterr_idx = read_to_pureread[l-k]+1
			previous_type = 0 #error
	
	read_start = pureread_start + read_start_diff
	read_err_locations_text = "[" + ",".join(str(x) for x in sorted(read_err_locations[0])) + "][" + ",".join(str(x) for x in sorted(read_err_locations[1])) + "][" + ",".join(str(x) for x in sorted(read_err_locations[2]))+"]"
	
	#Write pureread and read to respective output files
	pureheader = f">pureread{i}|{genome_text}|{strand_text}|{pureread_start}|{pureread_het_locations_text}|{pureread_actual_length}"
	header = f">read{i}|{genome_text}|{strand_text}|{read_start}|{read_err_locations_text}|{read_actual_length}" #eventually will add read_het_locations_text
	out_purereads.write(pureheader + "\n")
	out_purereads.write(pureread + "\n")
	out_reads.write(header + "\n")
	out_reads.write(read + "\n")

out_purereadhetblocks.close()
out_readerrblocks.close()
out_purereads.close()
out_reads.close()
