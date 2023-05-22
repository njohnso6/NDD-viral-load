#!/usr/bin/env bash
#SBATCH --mem 40G
#SBATCH --time 24:00:00

DIR='/data/CARD_AA/users/wellerca/ROSMAP_sample_FASTQs'

mkdir -p ${DIR} && cd ${DIR}


synapse get -r syn8612097 
