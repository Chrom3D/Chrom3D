import pybedtools
import sys
import re

raw_fname=str(sys.argv[1])
first_chr=str(sys.argv[2])
sec_chr=str(sys.argv[3])

with open(raw_fname) as rf:
  raw_file=rf.readlines()

raw_file=[x.strip() for x in raw_file]

hg19=pybedtools.helpers.get_chromsizes_from_ucsc("hg19")
first_chr_size = hg19[first_chr][1]
sec_chr_size = hg19[sec_chr][1]

bedpe = []
for line in raw_file:
  line=line.split()
  end1 = int(line[0]) + 1000000
  end2 = int(line[1]) + 1000000
  bedpe.append([first_chr,line[0],end1,sec_chr,line[1],end2,int(float(line[2]))])

for n in bedpe:
  if n[2] > first_chr_size:
    n[2] = first_chr_size
  if n[5] > sec_chr_size:
    n[5] = sec_chr_size
  n=map(str,n)
  print '\t'.join(n)
