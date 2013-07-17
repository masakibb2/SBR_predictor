#!/bin/bash

DIR=$1
echo "calclate PSSM in $DIR"

F=($(ls -1 ${DIR}/*.fasta))
cnt=0

echo "$F"


while [ $cnt -lt ${#F[@]} ]
do
  fas=${F[cnt]}
  fname=$(echo $fas | cut -d"/" -f2)
  blastpgp -d nr -i ${F[cnt]} -Q $1/pssm/$fname.pssm -O $1/out/$fname.out -h 0.001 -j 2 -m 7 -a 4 -J > ${F[cnt]}.out.xml
  cnt=$((cnt+1))
done 
