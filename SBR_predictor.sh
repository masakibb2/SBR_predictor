#!/bin/sh

./mkpssm.sh $1
python predict_pssm.py -i $1.pssm 

