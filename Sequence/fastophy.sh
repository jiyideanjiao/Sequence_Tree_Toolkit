#!/bin/bash -ve

#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1

for i in `ls *.codon`
do
id=$(basename $i .codon)
python convert_fasta2phylip.py $i $id.phy
done;
