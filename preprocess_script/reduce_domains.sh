cat $1 | groupBy -g 1,4 -c 2,3 -o min,max -full | awk -v OFS="\t" ' { print $1,$5,$6,$4}'

