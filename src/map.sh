#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks 10
#SBATCH --mem 40G
#SBATCH --time 0-12:00:00
#SBATCH --partition norm


module load bwa/0.7.17
module load samtools/1.17
module load GATK/4.3.0.0
module load pear/0.9.11



TMPDIR="/scratch/${SLURM_JOB_ID}/"
OUTDIR=${PWD}


map_paired_end() {
    local fastq_f=${1}
    local fastq_r=${2}
    local reference=${3}
    local tmpdir=${4}
    local outdir=${5}

    pear -f

    mkdir -p ${tmpDir}/${filestem} && cd ${tmpDir}/${filestem}
    
    # run PEAR to assemble paired end reads
    ${projectDir}/bin/pear-0.9.11-linux-x86_64/bin/pear -k \
    -f ${projectDir}/02_simulate_reads/${population}/${filestem}.F.fq.bz2 \
    -r ${projectDir}/02_simulate_reads/${population}/${filestem}.R.fq.bz2 \
    -o ${filestem}
    
    # map assembled reads
    bwa mem ${projectDir}/input_data/dgrp2.reference.fasta \
    ${filestem}.assembled.fastq | \
    samtools view -b \
    -o ${filestem}.assembled.bam -
    
    # map unassembled reads
    bwa mem ${projectDir}/input_data/dgrp2.reference.fasta \
    ${filestem}.unassembled.forward.fastq ${filestem}.unassembled.reverse.fastq | \
    samtools view -b \
    -o ${filestem}.unassembled.bam -
    
    # merge bam files with samtools
    samtools merge ${filestem}.merged.bam ${filestem}.assembled.bam ${filestem}.unassembled.bam
    
    # Sort bam file
    samtools sort ${filestem}.merged.bam > ${filestem}.sorted.bam
    
    # Add read groups and index .bam file
    echo "adding read groups"
    gatk AddOrReplaceReadGroups \
    --INPUT ${filestem}.sorted.bam \
    --OUTPUT ${filestem}.bam \
    --RGLB library1 \
    --RGPL illumina \
    --RGPU platform1 \
    --RGSM ${filestem}
    
    # index final output
    samtools index ${filestem}.bam
    
    # write final output to permenant disk location
    cp ${filestem}.bam* ${projectDir}/02_simulate_reads/${population}/
    
    # clean up
    cd ${projectDir}/02_simulate_reads && \
    rm -rf ${tmpDir}/${filestem}
}

export -f mapIndividual

for individual_n in $(seq ${first_ind} ${last_ind}); do
  filestem=${individual_n}.${chromosome}
  echo "mapping ${filestem}"
  mapIndividual ${filestem} ${tmpDir}
done

rm -rf ${tmpDir}