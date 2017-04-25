gtrack_w_lad=$1
inter_chr_sig=$2
out=$3

python inter_chrom_beads.py $gtrack_w_lad $inter_chr_sig | sort -k1,1 -k2,2n > $out
sed -i  '1s/^/##gtrack version: 1.0\n##track type: linked segments\n###seqid\tstart\tend\tid\tradius\tlamin\tedges\n/' $out
