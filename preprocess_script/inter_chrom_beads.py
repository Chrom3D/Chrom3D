import sys
import pybedtools as pb
import string
import random
import itertools

def file_generator(size=8, chars=string.ascii_uppercase + string.digits):
  fname = ''.join(random.choice(chars) for _ in range(size))
  fname = '/tmp/tmp' + fname 
  return fname

#midpoints = []
beads = []
beads_dict = {}
with open(sys.argv[1]) as inbed:
  for line in inbed:
    if line.startswith("#"):
      continue
    else:
      line2=line.rstrip().split('\t')[0:4]
      beads.append(line2)
      beads_dict[line2[3]] = line.rstrip()

end1 = []
end2 = []
with open(sys.argv[2]) as inchiasig:
  for line in inchiasig:
    line = line.split('\t')
    end1.append([line[0],str((int(line[1])+int(line[2]))/2),str(((int(line[1])+int(line[2]))/2)+1),line[0] + ":" + line[1] + "-" + line[2]])
    end2.append([line[3],str((int(line[4])+int(line[5]))/2),str(((int(line[4])+int(line[5]))/2)+1),line[3] + ":" + line[4] + "-" + line[5]])

end1_tf_name = file_generator()
end1_tf = open(end1_tf_name,'w')
for i in end1:
  i='\t'.join(i)
  print >> end1_tf, i
end1_tf.close()

end2_tf_name = file_generator()
end2_tf = open(end2_tf_name,'w')
for i in end2:
  i='\t'.join(i)
  print >> end2_tf, i
end2_tf.close()

bead_tf_name = file_generator()
bead_tf = open(bead_tf_name,'w')
for i in beads:
  i='\t'.join(i)
  print >> bead_tf, i
bead_tf.close()

end1_bed = pb.BedTool(end1_tf_name)
end2_bed = pb.BedTool(end2_tf_name)
bead_bed = pb.BedTool(bead_tf_name)

end1_with_bead = end1_bed.intersect(bead_bed,wao=True)
end2_with_bead = end2_bed.intersect(bead_bed,wao=True)
#print end1_with_bead

pre_interact_beads = {}
for i in xrange(0,len(end1_with_bead)):
  left = str(end1_with_bead[i]).split()[7]
  right = str(end2_with_bead[i]).split()[7]
  pre_interact_beads.setdefault(left,[]).append(right)
  pre_interact_beads.setdefault(right,[]).append(left)

interchr_interact_beads = {k: list(set(v)) for k, v in pre_interact_beads.items()}

#for key, value in beads_dict.iteritems():
  #print key, value

for key in beads_dict.iterkeys():
  if key in interchr_interact_beads.keys():
    if '.' in beads_dict[key]:
      print beads_dict[key].replace(".","")  + ';'.join(interchr_interact_beads[key])
    else:
      print beads_dict[key] + ";" + ';'.join(interchr_interact_beads[key])
  else:
    print beads_dict[key]
