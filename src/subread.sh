#!/usr/bin/env bash
#SBATCH --mem 120G
#SBATCH --time 12:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --gres lscratch:500


TMPDIR="/lscratch/${SLURM_JOB_ID}"
SIF=$(realpath ~/scrpio/src/scrpio.sif)
RFILE=$(realpath ~/scrpio/subread.R)
BCFILE=$(realpath ~/scrpio/737K-arc-v1.txt)
cd $TMPDIR
mkdir ref
mkdir fastq
mkdir db

cp ${RFILE} .
cp ${BCFILE} .

module load singularity/3

ref_dir='/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A'
db_dir=$(realpath ~/scrpio/db)
fastq_dir='/data/CARD_singlecell/Brain_atlas/NABEC_multiome/batch2/AN00013415_10X_RawData_Outs/UMARY-1544_scrn/HTWVGDSX5/'

singularity exec \
    -B ${ref_dir}:${PWD}/ref:ro \
    -B ${db_dir}:${PWD}/db:ro \
    -B ${fastq_dir}:${PWD}/fastq:ro \
    -H ${PWD} \
    ${SIF} \
    Rscript subread.R

cp UMARY-1544* /data/wellerca/scrpio/