#!/bin/bash
#SBATCH --job-name="genin"
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=michele.wyler@ivi.admin.ch
#SBATCH --array=1-13376

# module loading
TEMPDIR="$HOME/temp_Lysander"
export PATH="$PATH:$HOME"
mkdir -p $TEMPDIR/genin2_out

# file name
IDX=$(cat ${TEMPDIR}/nomiFasta.txt | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $0}')
NAME=`basename $IDX .fa`

# because of bad programming, segments need to be named _PB1 and all 8 with the same name
sed -e "s/>.*/>$NAME/g" $IDX | perl $HOME/geninHeaderPrep.pl - > $TEMPDIR/fasta_geninHeader/$NAME.fa 

genin2 -o $TEMPDIR/genin2_out/genin_$NAME.tsv $TEMPDIR/fasta_geninHeader/$NAME.fa

# remove again fasta
rm $TEMPDIR/fasta_geninHeader/$NAME.fa
