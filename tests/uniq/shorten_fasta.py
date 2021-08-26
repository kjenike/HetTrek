import sys

input_fasta = "t2t_chr22_26_27.fasta"

size=100000

with open(input_fasta, "r") as f:
    line = f.readline()
    while line:
        if line.count(">") > 1:
            print(line.strip())
        else:
            print(line[:size])
        line = f.readline()






