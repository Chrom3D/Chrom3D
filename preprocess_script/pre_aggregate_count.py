import sys

domain_fname=str(sys.argv[1])
matrix_fname=str(sys.argv[2])

with open(domain_fname) as df:
  domain_file=df.readlines()

with open(matrix_fname)  as mf:
  matrix_file=mf.readlines()

domain_dict = {}

for dl in domain_file:
  chrs,start,end,cat = dl.split()
  domain_dict[start,end] = chrs + ":" + str(start) + "-" + str(end)

for ml in matrix_file:
  bin1,bin2,raw_count = ml.split()
  mod_bin1 = ''
  mod_bin2 = ''
  for key in domain_dict.iterkeys():
    if int(key[0]) <= int(bin1) < int(key[1]):
      mod_bin1 = domain_dict[key]
    if int(key[0]) <= int(bin2) < int(key[1]):
      mod_bin2 = domain_dict[key]
  print mod_bin1, mod_bin2, raw_count 
