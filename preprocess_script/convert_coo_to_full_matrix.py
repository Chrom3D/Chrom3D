import numpy as np
import csv
import sys 

fname=str(sys.argv[1])
ofname=str(sys.argv[2])
chr=str(sys.argv[3])
bin_size=int(sys.argv[4])

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
  list1.extend([chr,i,i+bin_size])
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

with open(ofname,"wb") as csvfile:
  writer = csv.writer(csvfile,delimiter="\t")
  writer.writerows(list_of_list)
