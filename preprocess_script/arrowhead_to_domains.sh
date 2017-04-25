
# Arg1 - Arrowhead domainlist
# Arg2 - chromosome size --- Make sure it is sorted and chrY removed 


filename="${1%.*}"
tail -n +2 $1 | sort -k1,1 -k2,2n - | mergeBed -i - | awk '{print "chr"$1 "\t" $2 "\t" $3 "\tdomain"}' - > ${filename}.merged.bed
complementBed -i ${filename}.merged.bed -g $2 | cat - ${filename}.merged.bed | sort -k1,1 -k2,2n - | awk '{if($4=="domain") print $0; else print $1 "\t" $2 "\t" $3 "\tgap"}'  > ${filename}.domains

for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X}
do
  chromosome=chr$i
  awk '$1==chr' chr="$chromosome" ${filename}.domains >  ${filename}.$chromosome.domains
done
