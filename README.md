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
Once the test finishes running, outputs for this sample are located in `${outdir}/cromwell-executions/atac/${random-hash}`. Read more about the output directory specs [here](https://github.com/nicolerg/atac-seq-pipeline#4-output-directory-structure).

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

## 4. Output directory structure 

The outputs of the ENCODE ATAC-seq pipeline are extensive. **Table 4.3.1** briefly describes the logic behind the output directory structure. The subdirectories are listed in roughly the same order as they are executed within the pipeline. **Table 4.3.2** highlights some important files when the pipeline input is a single biological replicate. 

**Table 4.3.1 ENCODE ATAC-seq pipeline output directory structure.**

| Subdirectory                              | Description                             |
|-------------------------------------------|-----------------------------------------|
| `call-read_genome_tsv` | Locate the reference database files |
| `call-trim_adapter` | Trim adapter from raw FASTQ files with cutadapt |
| `call-bowtie2` | Align reads to genome with bowtie2 |
| `call-filter` | Filter raw BAM file (read quality, duplicates, and chrM) |
| `call-bam2ta` | Convert BAM to tagAlign (https://genome.ucsc.edu/FAQ/FAQformat.html#format15) |
| `call-spr` | Make 2 pseudoreplicates from single biological replicate |
| `call-macs2` | Call peaks with MACS2; make p-value and fold-change signal tracks |
| `call-macs2_pr1` | Call peaks for pseudoreplicate 1 (pr1) with MACS2 |
| `call-macs2_pr2` | Call peaks for pr2 with MACS2 | 
| `call-overlap_pr` | *Ignore when pipeline input is a single biological replicate (same as call-reproducibilty_overlap)* |
| `call-idr_pr` | *Ignore when pipeline input is a single biological replicate (same as call-idr_pr)* |
| `call-reproducibility_overlap` | Call high-confidence peaks using the “overlap” method |
| `call-reproducibility_idr` | Call high-confidence peaks using the “IDR” method (more stringent than the “overlap” method) |
| `call-ataqc` | Calculate QC metrics |
| `call-qc_report` | Compile QC metrics |

**Table 4.3.2 Important files in ENCODE ATAC-seq pipeline output directories.**

| Path | File | Description |
|------|------|-------------|
|call-trim_adapter/shard-0/execution | `${SID}_L001_R1_001.trim.merged.fastq.gz` `${SID}_L001_R2_001.trim.merged.fastq.gz` | Adapter-trimmed FASTQ files |
|call-bowtie2/shard-0/execution | ${SID}\_L001_R1_001.trim.merged.bam | Raw BAM file |
|call-bowtie2/shard-0/execution | ${SID}\_L001_R1_001.trim.merged.align.log | Alignment log file |
|call-filter/shard-0/execution | ${SID}\_L001_R1_001.trim.merged.nodup.bam | Filtered BAM file |
|call-filter/shard-0/execution | ${SID}\_L001_R1_001.trim.merged.dup.qc | Picard MarkDuplicates log file |
|call-macs2/shard-0/execution|${SID}\_L001_R1_001.trim.merged.nodup. tn5.pval0.01.300K.narrowPeak.gz | Raw peaks |
|call-macs2/shard-0/execution|${SID}\_L001_R1_001.trim.merged.nodup. tn5.pval0.01.300K.bfilt.narrowPeak.gz | Blacklist-filtered raw peaks (note that this file is identical to the “raw peaks” file for rat because the rat blacklisted region file is empty) |
|call-macs2/shard-0/execution|${SID}\_L001_R1_001.trim.merged.nodup. tn5.fc.signal.bigwig|Fold-change signal bigWig file|
|call-macs2/shard-0/execution|${SID}\_L001_R1_001.trim.merged.nodup. tn5.pval.signal.bigwig|P-value signal bigWig file|
|call-reproducibility_overlap/execution|optimal_peak.narrowPeak.gz|High-confidence peak set called with “overlap” method. **This will be the primary peak set used for downstream analysis**|
|call-reproducibility_overlap/execution|conservative_peak.narrowPeak.gz|Identical to `optimal_peak.narrowPeak.gz` in the same path when pipeline input is a single biological replicate|
|call-reproducibility_idr/execution|optimal_peak.narrowPeak.gz|High-confidence peak set called with “IDR” method. |
|call-reproducibility_idr/execution|conservative_peak.narrowPeak.gz|Identical to `optimal_peak.narrowPeak.gz` in the same path when pipeline input is a single biological replicate|
|call-ataqc/shard-0/execution|${SID}\_L001_R1_001.trim.merged.inserts.hist_graph.pdf|Insert size histogram. **Important for QC.**|
|call-ataqc/shard-0/execution|${SID}\_L001_R1_001.trim.merged_ large_tss-enrich.png|TSS enrichment plot. **Important for QC.**|
|call-ataqc/shard-0/execution|${SID}\_L001_R1_001.trim.merged_ tss-enrich.png|TSS enrichment plot. **Important for QC.**| 
|call-qc_report/execution|qc.html|Pretty HTML report of QC metrics. **Important for QC.**|
|call-qc_report/execution|qc.json|Detailed compilation of QC metrics. **Important for QC.**|

## 5. Resume the pipeline in the case of failure

Read [here](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/utils/resumer/README.md) about how you can resume the pipeline in the case of a failed run (instead of re-running the whole thing). 

## 6. Compile QC metrics across a batch of samples into a single TSV file

The pipeline generates a `qc.json` QC report for each sample. Follow the instructions [here](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/utils/qc_jsons_to_tsv/README.md) to compile QC metrics across samples into a single TSV file. 

**Table 5.1** provides definitions for a limited number of metrics included in the JSON QC reports. The full JSON report includes ~120 lines per sample; some lines are duplicates, and many metrics are irrelevant for running the pipeline with a single biological replicate. Therefore, this table includes a subset of relevant and important metrics when running a single replicate (as will usually be the case for Montgomery Lab ATAC-seq analysis).  

**Table 5.1 Description of relevant QC metrics.**

| Header | Metric | Definition/Notes |
|--------|--------|------------------|
|flagstat_qc|total|Read count from the trimmed FASTQ files|
|flagstat_qc|mapped_pct|Percent of reads that mapped|
|flagstat_qc|paired_properly_pct|Percent of reads that are properly paired|
|dup_qc|dupes_pct|Fraction (not percent) of paired reads that are duplicates|
|pbc_qc|NRF, PBC1|Both of these are measures of library complexity. Ideally, these values should be at least 0.8. If a sample performs poorly in the pipeline, a value of >0.8 indicates that better performance could result from sequencing the library more deeply.|
|nodup_flagstat_qc|total|Total number of reads after filtering the BAM file; number of reads input into peak calling|
|overlap_reproducibility_qc|N1|Number of peaks in the “overlap” peak set|
|idr_reproducibility_qc|N1|Number of peaks in the “IDR” peak set|
|overlap_frip_qc|FRiP_rep1-pr|Fraction of filtered reads in “overlap” peak set|
|ataqc|Read count from sequencer|Read count from the trimmed FASTQ files|
|ataqc|Read count successfully aligned||
|ataqc|Read count after filtering for mapping quality|Primarily, all alignments for a read pair are removed if each pair aligns the maximum number of allowed times, as dictated by `"atac.multimapping"`|
|ataqc|Read count after removing duplicate reads|Reads are removed after marking with picard MarkDuplicates|
|ataqc|Read count after removing mitochondrial reads (final read count)|Identical to `nodup_flagstat_qc:total`. This is misnamed, as the final filtered BAM file still includes a fraction of reads aligning to chrM. chrM reads are removed immediately before peak calling.|
|ataqc|picard_est_library_size|Estimated number of unique molecules in the library based on PE duplication|
|ataqc|Presence of NFR peak|Presence of nucleosome-free region|
|ataqc|Presence of Mono-Nuc peak|Presence of mono-nucleosomal peak (read lengths ~120-250)|
|ataqc|TSS_enrichment|Enrichment of reads in transcription start sites|

**NOTE:** As of Feb. 14, 2019, the following QC metrics are misreported when the multimapping parameter is set to a non-zero number. An issue has been submitted, and this MOP will be updated to reflect when the issue is fixed. Until then, please disregard these metrics. 
- flaqstat_qc: total (inflated) 
- flagstat_qc: mapped (inflated)    
- ataqc: read count from sequencer (inflated)   
- ataqc: read count successfully aligned (inflated) 
- ataqc: read count after filtering for mapping quality (inflated)  
- ataqc: read count after removing duplicate reads (inflated)   
- ataqc: read count after removing mitochondrial reads (final read count) - this does correspond to the number of reads in the filtered BAM file, but it misrepresents the number of reads removed due to chrM. chrM reads are not actually removed until immediately before peak calling (i.e. they are still present in the filtered BAM file `${SID}_L001_R1_001.trim.merged.nodup.bam`) 



