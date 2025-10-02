#!/bin/bash
#SBATCH --job-name="USgeno"
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=michele.wyler@ivi.admin.ch

# module loading
TEMPDIR="$HOME/temp_Lysander"
export PATH="$PATH:$HOME"

ls $TEMPDIR/singleGeno/A_*.fa | while read IDX; do
  NAME=`basename $IDX .fa`
  genoflu.py -f $IDX -n $TEMPDIR/genoflu_out/USclass_$NAME
  #rm $IDX
done
