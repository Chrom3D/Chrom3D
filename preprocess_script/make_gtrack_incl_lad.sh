bedtools intersect -c -a $1 -b $2 | awk '{if($7>=1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t1\t" $6; else  print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t0\t" $6}' > $3
sed -i  '1s/^/##gtrack version: 1.0\n##track type: linked segments\n###seqid\tstart\tend\tid\tradius\tlamin\tedges\n/' $3
