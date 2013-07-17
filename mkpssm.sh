#!/bin/sh
# set user name of tosto.

mkdir .mkpssm
mkdir .mkpssm/out .mkpssm/pssm

python split_fasta.py $1 .mkpssm
./calc_pssm.sh .mkpssm 

python joinpsm.py .mkpssm/pssm > $1.pssm

rm -r .mkpssm