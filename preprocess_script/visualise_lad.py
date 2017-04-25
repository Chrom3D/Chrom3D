import pybedtools
import sys
import tempfile
import re

cmm_fname=str(sys.argv[1])
lad_fname=str(sys.argv[2])

with open(cmm_fname) as cf:
  cmm_file=cf.readlines()

cmm_dict = {}

for line in cmm_file:
  if line.startswith("<marker "):
    line=line.rstrip()
    line=line.split()
    key=line[-1].replace("beadID=","").replace("\"","").replace("/","").replace(">","")
    cmm_dict[key] = line

tf = tempfile.NamedTemporaryFile(delete=False)
tf_name = tf.name
for key in cmm_dict.iterkeys():
  if("A" in key):
    temp_bed = re.split('[:-]',key)
    temp_bed = "\t".join(temp_bed)
    temp_bed = re.sub(r'(_[AB])',r'',temp_bed)
    print >> tf, temp_bed
tf.close()
  
beads = pybedtools.BedTool(tf_name)
sorted_beads = beads.sort()

lad_bed = pybedtools.BedTool(lad_fname)
sorted_lad_bed = lad_bed.sort()

beads_inter_lad = beads.intersect(sorted_lad_bed,c=True)

for l in beads_inter_lad:
  line=str(l).split()
  if(int(line[-1]) > 0):
    line1 = line[0] + "_A:" + line[1] + "-" + line[2]
    line2 = line[0] + "_B:" + line[1] + "-" + line[2]
    cmm_dict[line1][6] = "r=\"0.864\""
    cmm_dict[line1][7] = "g=\"0.000\""
    cmm_dict[line1][8] = "b=\"0.091\""
    cmm_dict[line2][6] = "r=\"0.864\""
    cmm_dict[line2][7] = "g=\"0.000\""
    cmm_dict[line2][8] = "b=\"0.091\""
  else:
    line1 = line[0] + "_A:" + line[1] + "-" + line[2]
    line2 = line[0] + "_B:" + line[1] + "-" + line[2]
    cmm_dict[line1][6] = "r=\"0.091\""
    cmm_dict[line1][7] = "g=\"0.045\""
    cmm_dict[line1][8] = "b=\"0.955\""
    cmm_dict[line2][6] = "r=\"0.091\""
    cmm_dict[line2][7] = "g=\"0.045\""
    cmm_dict[line2][8] = "b=\"0.955\""
#  print ' '.join(cmm_dict[line1])
#  print ' '.join(cmm_dict[line2])
#    beads_with_lad.append(line1)
#    beads_with_lad.append(line2)

#print sorted_beads
for line in cmm_file:
  if line.startswith("<marker "):
    line=line.rstrip()
    line=line.split()
    key=line[-1].replace("beadID=","").replace("\"","").replace("/","").replace(">","")
    print ' '.join(cmm_dict[key])
  else:
    print line.rstrip()
