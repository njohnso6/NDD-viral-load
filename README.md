# NDD-viral-load

# Retrieve and build viral database
NCBI hosts a curated collection of complete viral genomes to serve as our database for aligning reads.
```bash
mkdir -p DB
wget -O DB/viral.fasta.gz https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
wget -O DB/viral.features.gz https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz
module load bwa/0.7.17
bwa index DB/viral.fasta.gz
```

# Retrieve Reads
Configuration of the account used to retrieve data via the `Synapse`
client must be done interactively before submitting the retrieval job:

```bash
module load synapseclient/2.7.0

# Before continuing, ensure you have an authentication token generated
# at synapse.org which will be required before proceeding
# Visit synapse.org and follow the menus. 
# Be sure to include VIEW and DOWNLOAD permissions.
# Your Account
# └─ Account Settings
#    └─ Manage Personal Access Tokens
#       └─ Create New Tooken

# Authenticate synapse account via bash command line:
synapse config

# Leave login name BLANK
# Paste authentication token

# Submit job to retrieve data. Authentication will be inherited
# from the current interactive session's shell.
sbatch get-synapse-data.sh
```

Tada!

```bash
module load bwa/0.7.17

bwa mem -o 01_120405.sam DB/viral.fasta.gz 01_120405.r1.fastq.gz 01_120405.r2.fastq.gz 
```

## bmtagger setup

```bash
mkdir bmtagger && (cd bmtagger
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/README
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/README.bmfilter.txt
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/README.bmtagger.txt
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/README.bmtool.txt
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/README.standalone
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/bmfilter
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/bmtagger.sh
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/bmtool
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/extract_fullseq
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/srprism
chmod +x srprism
chmod +x bmtool
chmod +x extract_fullseq
chmod +x bmfilter
chmod +x bmtagger.sh
dos2unix bmtagger.sh)

export PATH=${PATH}:bmtagger

# make index for bmfilter
## first generate symlink to reference fasta
ln -s /fdb/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa reference.fa
bmtool -d reference.fa -o reference.bitmask -w 18

# make index for srprism
srprism mkindex -i reference.fa -o reference.srprism -M 7168
# generates files with prefix reference.srprism

## make db for blast
module load blast
makeblastdb -in reference.fa -dbtype nucl
```

Separately, retrieve `BBMap` from `sourceforge.net` and unpack the `.tar.gz`
```bash
tar -zxvf BBMap_39.01.tar.gz
```

# running bmtagger to remove human reads
```bash
SID='155_120423'
read1="/data/CARD_AA/data/ROSMAP_sample_FASTQs/${SID}.r1.fastq.gz"
read2="/data/CARD_AA/data/ROSMAP_sample_FASTQs/${SID}.r2.fastq.gz"
bash bmtagger.sh -b reference.bitmask -x reference.srprism \
            -T tmp -q1 -1 ${read1} -2 ${read2} -o ${SID}
```

# pull out reads not present in the bmtagger blacklist
```bash
bbmap/filterbyname.sh \
    in=${SID}.assembled.fastq \
    out=${SID}.merged.filter.fq \
    names=${SID}.filter

bbmap/filterbyname.sh \
    in=${SID}.assembled.fastq \
    out=${SID}.merged.filter.fq \
    names='H099WADXX130403:1:1101:10002:9880'



bbmap/filterbyname.sh \
    in=${SID}.unassembled.forward.fastq \
    in2=${SID}.unassembled.reverse.fastq \
    out=${SID}.unassembled.filter.fq \
    names=${SID}.filter

module load pear
pear -f ${read1} -r ${read2} -o ${SID}

```

2. take reads and run through bmtagger to generate filter list
3. Count number of human reads
4. filter out human reads with BBMap
5. Take remaining fastq and align to viral genome
6. count remaining unaligned reads

|   Sample  | Human reads |   Viral reads  | Other reads |
|-----------|-------------|----------------|-------------|
| 155_120423|  12345      |   123          |   512       |
etc



 # Map to human and convert to bam
bwa mem hg38 > ${ID}.bam


# Filter such that both reads unmapped
samtools view -b -f 12 -F 256 \
  ${ID}.bam > ${ID}.neither_mapped.bam

# split reads in bam
samtools sort -n  ${ID}.neither_mapped.bam ${ID}.neither_mapped.sorted.bam

# convert to fastq
bedtools bamtofastq \
-i ${ID}.neither_mapped.sorted.bam \
-fq ${ID}.nonhuman_r1.fastq \
-fq2 ${ID}.nonhuman_r2.fastq