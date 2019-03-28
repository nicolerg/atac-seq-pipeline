# Run the ENCODE ATAC-seq pipeline on SCG4

This repository is forked from the [ENCODE ATAC-seq pipeline repository](https://github.com/ENCODE-DCC/atac-seq-pipeline) to provide helper scripts and additional documentation specific for Montgomery Lab implementation of the pipeline on SCG. Please refer to the [ENCODE ATAC-seq pipeline repository](https://github.com/ENCODE-DCC/atac-seq-pipeline) for original documentation and use outside of this scope. 

## 1. Set up the environment  

On SCG4, execute the following code:
```bash
# download cromwell to $HOME
cd 
wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
chmod +rx cromwell-34.jar

screen -S install-encode # this takes a little while, so do it within a screen session 
module load miniconda/3 # load conda module

srcdir=/labs/smontgom/shared/atac_encode_pl/atac-seq-pipeline # points to pipeline code
bash ${srcdir}/conda/uninstall_dependencies.sh # to remove any existing pipeline env
bash ${srcdir}/conda/install_dependencies.sh # install dependencies 
```
Once the enviroment setup is complete, you may exit the screen session. You now have a `conda` environment called `encode-atac-seq-pipeline` with dependencies specific for running the ENCODE ATAC-seq pipeline. 

## 2. Run a test sample 

Run a small test sample with the following code. It takes about 1 hour to complete. This repository uses a modified version of [the original](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/examples/scg/ENCSR356KRQ_subsampled_scg_conda.sh) to allow for specification of an `${outdir}`. Otherwise, outputs are by default written to `${srcdir}` (where this repository resides). 

```bash
srcdir=/labs/smontgom/shared/atac_encode_pl/atac-seq-pipeline
outdir=/path/to/output
cd ${outdir} # pipeline outputs are writted to pwd
sbatch ${srcdir}/examples/scg/ENCSR356KRQ_subsampled_scg_conda.sh
```
Once the test finishes running, outputs for this sample are located in `${outdir}/cromwell-executions/atac/${random-hash}`. Read more about the output directory specs [here](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/output.md) or [here (my own more detailed version)](https://github.com/nicolerg/motrpac_atac_mop#43-output-directory-structure). 

## 3. Run the pipeline for a batch of samples (starting from FASTQ files)

Note that [Step 3.1](https://github.com/nicolerg/atac-seq-pipeline#31-make-a-json-template-for-your-project) and [Step 3.2](https://github.com/nicolerg/atac-seq-pipeline#32-generate-json-files-for-all-samples) are specific for starting the pipeline from FASTQ files (the most upstream step). It is also possible to start the pipeline from the following intermediate files:  
* Raw BAM files  
* Filtered/deduped BAM files  
* tagAlign files (downstream of filtered BAMs)  
Refer to https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/input.md to see how these inputs are defined in the JSON file if you do not want to start from FASTQs.    

### 3.1 Make a JSON template for your project

The input JSON file includes all of the specs needed to run the pipeline for a given sample. See documentation of all of the available parameters [here](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/input.md). One input JSON file is required for each sample (or set of replicates of the same sample). 

Xin and I have written a script to facilitate this process. The first step is to define a set of parameters that will be standard for all samples in this batch. Find a template with all possible options for paired-end samples [here](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/examples/template_pe.json).  

One of the parameters is `"atac.genome_tsv`", which defines the reference genome. ENCODE currently supports 4 default genomes, and I have installed a custom rat genome. All 5 genomes are available on SCG:

|genome|source|built from|TSV path|
|-|-|-|-|
|hg38|ENCODE|[GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)|`/reference/ENCODE/pipeline_genome_data/hg38_scg.tsv`|
|mm10|ENCODE|[mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)|`/reference/ENCODE/pipeline_genome_data/mm10_scg.tsv`|
|hg19|UCSC|[GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)|
|mm9|UCSC|[mm9, NCBI Build 37](<http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit>)|`/reference/ENCODE/pipeline_genome_data/mm9_scg.tsv`|
|rn6|Ensembl|[Rnor_6.0 (GCA_000001895.4), Release 95](<ftp://ftp.ensembl.org/pub/release-95/fasta/rattus_norvegicus/dna/>)|`/labs/smontgom/shared/atac_encode_pl/REFERENCES/rn6/rn6.tsv`|

Follow these steps to generate a standard set of parameters for your batch of samples:

1. Copy `/labs/smontgom/shared/atac_encode_pl/atac-seq-pipeline/examples/template_se.json` (single-end reads) or `/labs/smontgom/shared/atac_encode_pl/atac-seq-pipeline/examples/template_pe.json` (paired-end reads) to new file.  
2. Remove or adjust the parameters as desired.  
* **Do** include the following important parameters, at a minimum: 
    * `"atac.genome_tsv"`: Choose one of the options from the table above.
    * `"atac.auto_detect_adapter"`: `false` is the default. I like to set this to `true` and let the pipeline determine the adapter sequence itself.
* `"atac.paired_end": true,`: `make_json.sh` assumes the inputs are paired-end reads; you will have to make some adjustments to `make_json.sh` if this is not the case. 
    * `"atac.mito_chr_name"`: `"chrM"` is the default. If you are using the `rn6` rat genome, define this as `"MT"`. Otherwise, you can ignore this.
    * `"atac.keep_irregular_chr_in_bfilt_peak"`: `false` is the default. If you are using the `rn6` rat genome, set this to `true`. Otherwise, you can ignore this. 
    * I would recommend defining all of the `"atac.*_cpu"` and `"atac.*_mem_mb"` values yourself since it makes it easier to request the necessary resources when you run the pipeline with `sbatch`. 
    * If you do not care about looking at sample-level ATAC-seq QC metrics, set `"atac.disable_ataqc"` to `true`. This module is computationally intensive.
* Do **not** include the following parameters, as they will be specified by `make_json.sh`:
    * `"atac.fastqs_rep1_R1"`
    * `"atac.fastqs_rep1_R2"`
    * `"atac.title"`
3. Delete the outside brackets (`{` and `}`) but none of the leading whitespace.
4. Save this file somewhere relevant to your project. I will refer to it as `${json_base}`. Find an example here: `/labs/smontgom/shared/atac_encode_pl/atac-seq-pipeline/nrg_src/example_base.json`.  

### 3.2 Generate JSON files for all samples 

Now that you have defined the runtime parameters, use `/labs/smontgom/shared/atac_encode_pl/src/make_json.sh` to generate input JSON files for all of your samples:

1. Copy `/labs/smontgom/shared/atac_encode_pl/atac-seq-pipeline/nrg_src/make_json.sh` to one of your local folders.
2. Define `fastq_dir` (the path to your FASTQ files) and `json_dir` (the desired output directory for the input JSON files) in `make_json.sh`. 
3. Run `make_json.sh` and spot-check one of the JSON files in `${json_dir}`.  

Note that this script will have to be adjusted if your FASTQ files do not follow the default naming convention from `bcl2fastq` (`${SAMPLE_NAME}_L00*_R*_001.fastq.gz`)

### 3.3 Make an `sbatch` script for each JSON file  

`sbatch` submits a batch script to SLURM. `sbatch` parameters cannot contain variables, but I find it useful to specify `sbatch` options like `--job-name` and the path to the error/output log file. Therefore, I have written an `sbatch` template script with palceholders for the `sbatch` parameters I want to change for each job (`/labs/smontgom/shared/atac_encode_pl/atac-seq-pipeline/nrg_src/run_encode_single.sh`). Then I use a short script to substitute the placeholders with the real values and write an `sbatch` script for each sample (`/labs/smontgom/shared/atac_encode_pl/atac-seq-pipeline/nrg_src/write_batch_run_encode_single.sh`). 

In a copy of `/labs/smontgom/shared/atac_encode_pl/atac-seq-pipeline/nrg_src/write_batch_run_encode_single.sh`, do the following:
* Define `${outdir}` as the output directory for the pipeline (`cromwell-executions`). SLURM logs (with both `stderr` and `stdout`) will also be written to `${outdir}/slurm_out`. 
* Define `${indir}` as the path to the JSON files. This script is currently written to write the individual `sbatch` scripts to `${indir}` as well. If that is not desired, adjust the code accordingly. 

Spot-check one of the `sbatch` scripts, named `${indir}/run_encode_single_${SAMPLE_NAME}.json`. If everything looks good, go ahead a submit jobs for all of your samples.

### 3.4 Run the pipeline  

Run `sbatch ${indir}/run_encode_single_${SAMPLE_NAME}.json` for every job you want to submit to the queue. 

## 4. Get familiar with the outputs 

Refer to the [ENCODE documentation of pipeline outputs](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/output.md) or the [documentation I provided](https://github.com/nicolerg/motrpac_atac_mop#43-output-directory-structure) for the MoTrPAC implementation of the pipeline.

## 5. Resume the pipeline in the case of failure

Read [here](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/utils/resumer/README.md) about how you can resume the pipeline in the case of a failed run (instead of re-running the whole thing). 

## 6. Compile QC metrics across a batch of samples into a single TSV file

Follow the instructions [here](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/utils/qc_jsons_to_tsv/README.md). 

