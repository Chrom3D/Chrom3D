python make_gtrack.py $1 $2 | sort -k1,1 -k2,2n > $3
sed -i  '1s/^/##gtrack version: 1.0\n##track type: linked segments\n###seqid\tstart\tend\tid\tradius\tedges\n/' $3
