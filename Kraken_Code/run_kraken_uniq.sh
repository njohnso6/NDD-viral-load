#!/usr/bin/env bash

export KRAKEN2_DB_PATH="/data/CARD_AA/user/johnsonnicl/ndd_virus/kraken_uniq_db"

./krakenuniq_installed/krakenuniq --db ./krakenuniq_db/ --report reports/R24_131017.kreport --preload-size=50G --threads=10 --gzip-compressed --paired /data/CARD_AA/data/ROSMAP_sample_FASTQs/nonhuman/R24_131017.nonhuman.r1.fastq.gz /data/CARD_AA/data/ROSMAP_sample_FASTQs/nonhuman/R24_131017.nonhuman.r2.fastq.gz  > R24_131017.kraken
