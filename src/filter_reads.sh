#!/usr/bin/env bash
#SBATCH --mem 40G
#SBATCH --time 2:59:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH --gres lscratch:50
#SBATCH --partition quick,norm

START=$SECONDS


## Directory definitions #############################################

# pre-computed reference (and associated bwa indices)
REF='/fdb/bwa/indexes/hg38.fa'

# Input fastq.gz files
INDIR='/data/CARD_AA/data/ROSMAP_sample_FASTQs/'

# Output directory to save filtered fastq.gz
OUTDIR='/data/CARD_AA/data/ROSMAP_sample_FASTQs/nonhuman/'
mkdir -p ${OUTDIR}

# File containing list of sample IDs
SAMPLEFILE=${PWD}/samples.txt


## Define job array ##################################################
if [ ! -f "${SAMPLEFILE}" ]; then
    find ${INDIR} -maxdepth 1 -name '*.r1.fastq.gz' | \
    sort | sed 's/\.r1.fastq.gz$//g' > ${SAMPLEFILE}
fi

# Working directory in /lscratch
TMPDIR="/lscratch/${SLURM_JOB_ID}/"




## Sample Definitions ################################################

# Define N-th sample name
if [ "${SLURM_ARRAY_TASK_ID}" == '' ]; then
    N=${1}
else
    N=${SLURM_ARRAY_TASK_ID}
fi

if [ "${N}" == '' ]; then
    echo 'ERROR: No job ID specified as argument 1 or job array number'
    exit 1
fi

# Get N-th sample name and filestem
FS=$(sed -n "${N}"p ${SAMPLEFILE})
ID=$(basename ${FS})
echo "Working with sample #${N}: ${ID}"

# Define full sample name for read pair
READ1=${FS}.r1.fastq.gz
READ2=${FS}.r2.fastq.gz

# Stop if filtered fastq already exists
if [ -f "${OUTDIR}/${ID}.nonhuman.r1.fastq.gz" ]; then
    echo INFO: SKIPPING SAMPLE ${ID}
    echo INFO: "${OUTDIR}/${ID}.nonhuman.r1.fastq.gz" already exists
    exit 0
fi




## RUN ###############################################################

# Load modules
module load bwa/0.7.17
module load samtools/1.17


# Work in /lscratch tmp directory
cd ${TMPDIR}


# bwa index ref if not already done.
# these SHOULD already be pre-computed on biowulf
if [ ! -f "${REF}" ]; then
    bwa index ${REF}
else
    echo 'reference already indexed by bwa'
fi


# run bwa
bwa mem -t 12 ${REF} ${READ1} ${READ2} -o ${ID}.bam
# ~40 minutes with 12 CPUs

# Extract unmapped reads
samtools view -h -b -f 12 -F 256 ${ID}.bam | \
samtools sort > ${ID}.unmapped.bam
# ~1 minute with 12 CPUs


# Convert back to fastq
samtools fastq \
    -1 ${ID}.unmapped.r1.fastq \
    -2 ${ID}.unmapped.r2.fastq \
    ${ID}.unmapped.bam

# Gzip and output to permanent storage
gzip -c ${ID}.unmapped.r1.fastq > ${ID}.nonhuman.r1.fastq.gz ${OUTDIR}
gzip -c ${ID}.unmapped.r2.fastq > ${ID}.nonhuman.r2.fastq.gz ${OUTDIR}

END=$SECONDS
ELAPSED=$(( end - start ))
echo "Completed in ${ELAPSED} seconds"