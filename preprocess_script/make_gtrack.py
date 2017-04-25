import sys
import re

sig_fname=str(sys.argv[1])
domain_fname=str(sys.argv[2])

with open(sig_fname) as f:
  sig_file=f.readlines()

with open(domain_fname) as dom:
  domain_file=dom.readlines()

sig_file=[x.strip() for x in sig_file]
domain_file=[x.strip() for x in domain_file]

sig_list = []

sig_dict = {}
 
for line in sig_file:
  chr_1,start1,end1,chr_2,start2,end2,pval = line.split()
  key=chr_1 + ":" + start1 + "-" + end1
  val=chr_2 + ":" + start2 + "-" + end2
  sig_dict.setdefault(key,[]).append(val)


combined_dict = {}

for line in domain_file:
  char,start,end,cat = line.split()
  domain=char + ":" + start + "-" + end
  
  if domain in sig_dict:
    if domain in sig_dict.keys():
       #print domain, sig_dict[domain]
       combined_dict[domain] = sig_dict[domain]
    for key, value in sig_dict.iteritems():
       if domain in value:
         #print domain, key
         if domain in combined_dict.keys():
            combined_dict[domain].append(key)
  else:
    #print domain + "\t."
    combined_dict[domain] = "."


#print "##gtrack version: 1.0"
#print "##track type: linked segments"
#print "###seqid start	end	id	radius	edges"
for key, value in combined_dict.iteritems():
  temp_bed = re.split('[:-]',key)
  temp_bed = "\t".join(temp_bed)
  temp_edge = ";".join(value) 
  print temp_bed + "\t" +  key + "\t" + "1\t"+ temp_edge
