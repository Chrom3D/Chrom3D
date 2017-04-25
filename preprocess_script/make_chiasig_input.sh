python pre_aggregate_count.py $1 $2 | awk '{mod[$1" "$2]+=$3}END{for (i in mod){print i,mod[i]}}' - | sed 's/:/\t/g' | sed 's/-/\t/g' | sed 's/ /\t/g' | sort -k1,1 -k2,2n | awk '$2!=$5' - 
