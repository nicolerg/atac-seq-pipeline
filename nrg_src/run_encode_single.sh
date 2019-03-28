#!/bin/bash

#SBATCH --account=smontgom
#SBATCH --time=48:00:00
#SBATCH --job-name=_JID_
#SBATCH --mail-type=END,FAIL
#SBATCH --error=_OUTDIR_/slurm_out/slurm-%j-%x.out
#SBATCH --output=_OUTDIR_/slurm_out/slurm-%j-%x.out

# max number of cpus for each pipeline
#  should be <= NUM_CONCURRENT_TASK x "atac.bowtie2_cpu" in input JSON file
#  since bowtie2 is a bottlenecking task in the pipeline
#SBATCH --cpus-per-task=8

# total amount of memory
#  depends on the size of your FASTQs
#  but should be <= NUM_CONCURRENT_TASK x 20GB for big samples
#  or <= NUM_CONCURRENT_TASK x 10GB for small samples
#SBATCH --mem=40G

# do not touch these settings
#  number of tasks and nodes are fixed at 1
#SBATCH --ntasks=1
#SBATCH --nodes=1

set -e 

module load java
module load miniconda/3
source activate encode-atac-seq-pipeline

# check pysam
python -c "import pysam"

NUM_CONCURRENT_TASK=2

srcdir=/labs/smontgom/shared/atac_encode_pl/atac-seq-pipeline
INPUT=_INDIR_/_JID_.json
PIPELINE_METADATA=_OUTDIR_/_JID__metadata.json 

cd _OUTDIR_

java -jar -Dconfig.file=${srcdir}/backends/backend.conf \
/home/nicolerg/cromwell-34.jar run ${srcdir}/atac.wdl -i ${INPUT} -m ${PIPELINE_METADATA}


