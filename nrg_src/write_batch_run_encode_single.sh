#!/bin/bash

set -e

outdir=/labs/smontgom/nicolerg/AFRICAN_OMICS # pipeline outputs will be written to ${outdir}/cromwell-executions/atac
indir=/labs/smontgom/nicolerg/AFRICAN_OMICS/JSON # path to JSON files; where sbatch scripts will be written

mkdir -p ${outdir}
mkdir -p ${outdir}/slurm_out

for json in `ls ${indir} | grep "json"`; do 
	jid=$(echo $json | sed "s/\.json//")
	sed -e "s:_OUTDIR_:$outdir:g" -e "s:_JID_:$jid:g" -e "s:_INDIR_:$indir:g" run_encode_single.sh > ${indir}/run_encode_single_$jid.sh
	
	#sbatch ${indir}/run_encode_single_${jid}.json # uncomment this line if you also want this script to submit all jobs to the queue
done


