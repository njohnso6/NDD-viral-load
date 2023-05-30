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


# Filter out reads aligning to human genome
```bash
sbatch --array=1-639%5 src/filter_reads.sh
```
