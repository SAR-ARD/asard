#!/usr/bin/env bash
#SBATCH -J ers-asar-nrb-tests
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --cpus-per-task=16
#SBATCH --mem=45gb
#SBATCH --mail-type=all
#SBATCH --mail-user=<your mail address>
#SBATCH --export=NONE
#SBATCH --time=01:00:00
#SBATCH --array=0-7

set -euo pipefail

module load slurm_setup
module load micromamba/1.5.7
module load glib
source $HOME/.bashrc
micromamba activate s1ard

export ASARD_TESTDATA="<test data storage location>"

PARAM_IDS=("ASAR-APP" "ASAR-APS" "ASAR-IMP" "ASAR-IMS" "ERS1-IMP" "ERS1-IMS" "ERS2-IMP" "ERS2-IMS")
ID="${PARAM_IDS[$SLURM_ARRAY_TASK_ID]}"

coverage run --parallel-mode --module pytest -k "$ID" test_processing.py
