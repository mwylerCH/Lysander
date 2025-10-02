#!/bin/bash
#SBATCH --job-name="ALN"
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=michele.wyler@ivi.admin.ch
#SBATCH --array=1-8

# module loading
module load Anaconda3
source activate famsa

TEMPDIR="$HOME/temp_Lysander"
export PATH="$PATH:$HOME"

# file name
IDX=$(cat ${TEMPDIR}/fileToAlign.txt | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $0}')
NAME=`basename $IDX .fa`
DIST=`echo $IDX | sed "s/.aln/.distMatrix/"`

# make alignment first
srun famsa -t 20 $IDX $TEMPDIR/${NAME}.aln

# trim alignment (aka remove UTRs)
trimal -in $TEMPDIR/${NAME}.aln -out $TEMPDIR/${NAME}_trimmed.aln -cons 60 -gt 0.9


# run dist matrix
srun famsa -t 20 -dist_export -square_matrix $TEMPDIR/${NAME}_trimmed.aln $TEMPDIR/${NAME}_2.DistMatrix
#srun mafft --thread 30 --retree 1 --adjustdirection --maxiterate 0 $IDX > $TEMPDIR/${NAME}.aln

# run dist matrix
#srun /storage/homefs/mw23o336/gotree matrix -t 60 -i $IDX -o $DIST
