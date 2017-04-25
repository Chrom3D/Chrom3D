import numpy as np
import csv
import sys 

fname=str(sys.argv[1])
ofname=str(sys.argv[2])
chr=str(sys.argv[3])
bin_size=int(sys.argv[4])

chr_size = {'chr1':249250621,
	'chr2':243199373,
	'chr3':198022430,
	'chr4':191154276,
	'chr5':180915260,
	'chr6':171115067,
	'chr7':159138663,
	'chr8':146364022,
	'chr9':141213431,
	'chr10':135534747,
	'chr11':135006516,
	'chr12':133851895,
	'chr13':115169878,
	'chr14':107349540,
	'chr15':102531392,
	'chr16':90354753,
	'chr17':81195210,
	'chr18':78077248,
	'chr19':59128983,
	'chr20':63025520,
	'chr21':48129895,
	'chr22':51304566,
	'chrX':155270560}


with open(fname) as f:
  file=f.readlines()

file=[x.strip() for x in file]

d = {}

for line in file:
  x,y,n = line.split()
  d[x,y] = n

get_max_coord = []
for key, value in d.iteritems():
  get_max_coord.append(int(key[1]))

chr_max_value = max(get_max_coord) + bin_size 

query = list(np.arange(0,chr_max_value,bin_size))

list_of_list = []
for i in query:
  list1 = []
  col3 = i+bin_size
  list1.extend([chr,i,col3])
  for j in query:
    keytuple = (str(i),str(j)) 
    if keytuple in d:
      if(i==j):
        list1.append(0)
      else:
        list1.append(d[keytuple])
    else:
      reverse_keytuple = (str(j),str(i))
      if reverse_keytuple in d:
        if(i==j):
          list1.append(0)
        else:
          list1.append(d[reverse_keytuple])
      else:
        list1.append(0)
  list_of_list.append(list1)

list_of_list[-1][2] = chr_size[chr]

with open(ofname,"wb") as csvfile:
  writer = csv.writer(csvfile,delimiter="\t")
  writer.writerows(list_of_list)
