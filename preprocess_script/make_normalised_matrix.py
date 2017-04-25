import sys



fname=str(sys.argv[1]) #"chr19_50kb.RAWobserved"
norm_file=str(sys.argv[2]) #"chr19_50kb.KRnorm"

with open(fname) as f:
  file=f.readlines()

file=[x.strip() for x in file]

d = {}

for line in file:
  x,y,n = line.split()
  d[x,y] = n

#print d["700000","750000"]

with open(norm_file) as norm_fh:
  norm_vector=norm_fh.readlines()

norm_vector = [x.strip() for x in norm_vector]

norm_dict = {}
bin_size = 50000
for key,value in d.iteritems():
  #key=('750000','750000')
  factor1_idx = (int(key[0])/bin_size) 
  factor2_idx = (int(key[1])/bin_size)
  factor1 = norm_vector[factor1_idx]
  factor2 = norm_vector[factor2_idx]
  #print key[0], factor1_idx, factor1, key[1], factor2_idx, factor2
  if(factor1 == 'NaN' or factor2 == 'NaN'):
    norm_dict[key] = 0
  else:
    norm_dict[key] = round(float(d[key])/(float(factor1)*float(factor2)),2)
#print norm_dict

for norm_key, norm_value in norm_dict.iteritems():
  print norm_key[0] + "\t" + norm_key[1] + "\t" + str(norm_dict[norm_key])
