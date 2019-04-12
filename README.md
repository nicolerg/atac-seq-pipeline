# MoTrPAC ATAC-Seq QC and Analysis Pipeline 

This repository provides the MoTrPAC-specific implementation of the [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline), which mirrors the most recent version of the MoTrPAC ATAC-seq QC and Analysis Pipeline MOP, available [here](https://docs.google.com/document/d/1vnB7ITAKnaZYc3v_FCdaDu3z-JXeDncRk5GnqzQVwRw/edit#heading=h.tjbixx8yyd33). Updates to the ENCODE pipeline are merged with this pipeline as necessary.  

### Important references:
- GitHub repository for the ATAC pipeline: https://github.com/ENCODE-DCC/atac-seq-pipeline
- ENCODE ATAC-seq pipeline documentation: https://www.encodeproject.org/atac-seq/
- ENCODE data quality standards: https://www.encodeproject.org/atac-seq/#standards 
- ENCODE terms and definitions: https://www.encodeproject.org/data-standards/terms/

### Table of Contents:
1. [Prepare FASTQ files](https://github.com/nicolerg/motrpac-atac-seq-pipeline#1-prepare-fastq-files)     
    1.1. [BCL to FASTQ conversion with bcl2fastq](https://github.com/nicolerg/motrpac-atac-seq-pipeline#11-bcl-to-fastq-conversion-with-bcl2fastq)     
    1.2. [FastQC on raw FASTQ files](https://github.com/nicolerg/motrpac-atac-seq-pipeline#12-fastqc-on-raw-fastq-files)   
2. [Collect metadata associated with each sample](https://github.com/nicolerg/motrpac-atac-seq-pipeline#2-collect-metadata-associated-with-each-sample)      
3. [Set up ENCODE ATAC-seq pipeline](https://github.com/nicolerg/motrpac-atac-seq-pipeline#3-set-up-encode-atac-seq-pipeline)  
    3.1 [Install and test ENCODE ATAC-seq pipeline and dependencies](https://github.com/nicolerg/motrpac-atac-seq-pipeline#31-install-and-test-encode-atac-seq-pipeline-and-dependencies)  
    3.2 [Build genome databases](https://github.com/nicolerg/motrpac-atac-seq-pipeline#32-build-genome-databases)  
4. [Run the ENCODE ATAC-seq pipeline](https://github.com/nicolerg/motrpac-atac-seq-pipeline#4-run-the-encode-atac-seq-pipeline)     
    4.1 [Generate JSON files](https://github.com/nicolerg/motrpac-atac-seq-pipeline#41-generate-json-files)   
    4.2 [Run the pipeline](https://github.com/nicolerg/motrpac-atac-seq-pipeline#42-run-the-pipeline)     
    4.3 [Output directory structure](https://github.com/nicolerg/motrpac-atac-seq-pipeline#43-output-directory-structure)
5. [Compile important QC metrics](https://github.com/nicolerg/motrpac-atac-seq-pipeline#5-compile-important-qc-metrics)   
6. [Flag problematic samples](https://github.com/nicolerg/motrpac-atac-seq-pipeline#6-flag-problematic-samples)   

## 1. Prepare FASTQ files
### 1.1 BCL to FASTQ conversion with bcl2fastq
Each GET site (Stanford and Mount Sinai) is responsible for sequencing the library and obtaining the demultiplexed FASTQ files for each sample. If sequencing is performed with NovaSeq, raw data is output as BCL files, which must be demultiplexed and converted to FASTQ files with `bcl2fastq` (version 2.20.0). 

`bcl2fastq` must be run with a `SampleSheet.csv` file that does not include the `Adapter` and `AdapterRead2` settings. This will prevent `bcl2fastq` from automatically performing adapter trimming, which provides us with FASTQ files that include the fullest form of the raw data. Adapter trimming is performed downstream. 

`bcl2fastq v2.20.0` can be downloaded directly from Illumina [here](https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html). Use the following command, where `$seqDir` is the path to the sequencing output directory (e.g `181205_NB551514_0071_AHFHLGAFXY`) and `$outDir` is the path to the desired output directory:

```bash
bcl2fastq \   
     --sample-sheet /path/to/SampleSheet.csv
     --runfolder-dir $seqDir \
     --output-dir $outDir 
```
This command will generate two FASTQ files (one for each read direction) per sample per lane, e.g. `${sample_ID}_L${lane_number}_R{1,2}_001.fastq.gz`. 

Collect the file `laneBarcode.html` file (`$seqDir/Reports/html/*/all/all/all/laneBarcode.html`) for primary inspection at each GET site to decide if the whole library needs to be resequenced or some samples need to be re-pooled in a future library due to relatively small number of reads.  

### 1.2 FastQC on raw FASTQ files 
We will utilize FastQC and MultiQC to analyze the quality of the sequencing. The intent is not to hold the analysis pipeline until after this process is finished. It is to capture the metrics so that problematic runs can be flagged and explored. 

#### 1.2.1 Run FastQC on the FASTQ files with the following command:
```bash
fastqc \
    -o $output_directory \
    -t $num_threads \
    -f fastq \
    -a $file_with_adapters \
    ${sample_ID}_L001_R1_001.fastq.gz
    ${sample_ID}_L001_R2_001.fastq.gz
    ${sample_ID}_L002_R1_001.fastq.gz
    ...
    ${sample_ID}_L00N_R2_001.fastq.gz
```
The above command will produce a single zipped folder for each FASTQ file. The folder contains metrics about read quality. These metrics could be made visible for general consumption or they could be held for QA/QC review when required (i.e. in case of run failures, etc.). 

#### 1.2.2 Run MultiQC on the results from FastQC to compile the metrics for all FASTQ files with the following command:
```bash
multiqc \
    -n $report_name \
    -o $output_directory \
    $fastqc_output_directory
```
The above command searches the files in the the target directory (`$fastqc_output_directory`) and outputs a single folder and a single HTML report that summarizes the FastQC results from all FASTQ files into easy-to-read tables and figures.

## 2. Collect metadata associated with each sample 
The GET sites are responsible for submitting metadata associated with each sequencing run to the BIC. These metadata include sample-level metrics, run-level metrics, and metrics required for the BIC to process the ATAC-seq data with the ENCODE pipeline. Hence, the metadata should include (at minimum) the information in the attribute dictionary available [here](https://docs.google.com/document/d/1xlFiax4MTSzNZS3SpG6i3Z3XwuGkMPgONV1QPTObjcA/edit#heading=h.sqhy9p63uf9b). 

## 3. Set up ENCODE ATAC-seq pipeline
All steps in this section must only be performed once. After dependencies are installed and genome databases are built, skip to [here](https://github.com/nicolerg/motrpac_atac_mop#4-run-the-encode-atac-seq-pipeline).

### 3.1 Install and test ENCODE ATAC-seq pipeline and dependencies

This pipeline supports many cloud platforms and cluster engines. It also supports `docker`, `singularity` and `Conda` to resolve complicated software dependencies for the pipeline. A tutorial-based instruction for each platform is provided to understand how to run pipelines. There are special instructions for two major Stanford HPC servers (SCG4 and Sherlock).  

The BIC will be using the Google Cloud Platform implementation of the ENCODE ATAC-seq pipeline, but the GET/CAS sites may preferentially run the pipeline on a different platform before the BIC implementation is production-ready. For example, Stanford has run the pipeline on both the “Local System with Conda” and "Stanford SCG4" platforms. The analysis steps of the pipeline are the same regardless of the platform on which the pipeline runs. The only difference is in the environment-dependent installation of software dependencies and genome databases.  

Follow platform-specific instructions to clone this repository, install dependencies, and run a test sample on the platform of your choice:  

* Cloud platforms
  * Web interface
    * [DNAnexus Platform](docs/tutorial_dx_web.md)
  * CLI (command line interface)
    * [Google Cloud Platform](docs/tutorial_google.md)
    * [DNAnexus Platform](docs/tutorial_dx_cli.md)
* Stanford HPC servers (CLI)
  * [Stanford SCG4](docs/tutorial_scg.md)
  * [Stanford Sherlock 2.0](docs/tutorial_sherlock.md)
* Cluster engines (CLI)
  * [SLURM](docs/tutorial_slurm.md)
  * [Sun GridEngine (SGE/PBS)](docs/tutorial_sge.md)
* Local Linux computers (CLI)
  * [Local system with `singularity`](docs/tutorial_local_singularity.md)
  * [Local system with `docker`](docs/tutorial_local_docker.md)
  * [Local system with `Conda`](docs/tutorial_local_conda.md)
* Local Windows computers (CLI)
  * [Windows 10 Pro with `docker`](docs/tutorial_windows_docker.md)
  * [Windows 10 Pro/Home with `Conda`](docs/tutorial_windows_conda.md)

### 3.2 Build MoTrPAC-specific genome databases
#### 3.2.1 Build a human genome (hg38) database 

1. [Install Conda](https://conda.io/miniconda.html). Skip this if you already have equivalent Conda alternatives (Anaconda Python). Download and run the [installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh). Agree to the license term by typing `yes`. It will ask you about the installation location. On most cluster systems, we recommend installing it outside of your `$HOME` directory since its filesystem is slow and has very limited space. At the end of the installation, choose `yes` to add Miniconda's binary to `$PATH` in your BASH startup script. Note that it is not recommended to install your own version of Conda if you are working on Stanford SCG4. Instead, load Conda with `module load miniconda/3`. 
    ```bash
    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh
    ```
    OR, if working on Stanford SCG4:
    ```bash
    $ module load miniconda/3
    ```
    
2. Install Conda dependencies.
    ```bash
    $ bash conda/uninstall_dependencies.sh  # to remove any existing pipeline env
    $ bash conda/install_dependencies.sh
    ```

3. Specify a destination directory and install the MoTrPAC hg38 reference with the following command. We recommend not to run this installer on a login node of your cluster. It will take >8GB memory and >2h time.
    ```bash
    $ outdir=/path/to/reference/genome
    $ bash conda/build_genome_data.sh motrpac_hg38 ${outdir}
    ```
    
#### 3.2.2 Build a rat (rn6) database 

1. Perform steps 1 and 2 from 3.2.1 if you have not done so once already. (If you have already done this to build the human reference, you do not need to do it again). If you are working on SCG4, run `module load miniconda/3` to load Conda. 

2. Specify a destination directory and install the MoTrPAC rn6 reference with the following command. We recommend not to run this installer on a login node of your cluster. It will take >8GB memory and >2h time.
    ```bash
    $ outdir=/path/to/reference/genome
    $ bash conda/build_genome_data.sh motrpac_rn6 ${outdir}
    ```
    








## 4. Run the ENCODE ATAC-seq pipeline
### 4.1 Generate JSON files
A JSON file specifying input parameters to run the pipeline is required for each sample. Find documentation of definable parameters [here](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/input.md).

Two files are provided to help:

1. `make_json.sh` found [here](https://github.com/nicolerg/motrpac_atac_mop/blob/master/files/make_json.sh) is a script to quickly generate one JSON file per sample. 

- Some parameters must be customized. These include:
    - `fastq_dir`: path to raw FASTQ files
    - `json_dir`: directory where JSON files will be written
    - `genome_ref`: path to reference genome built in step 3.2 (**note that this must correspond to the correct species, i.e. human or rat**)
    - `lanes_per_sample`: number of lanes per sample
- Other parameters may be changed. These include all of those related to memory and CPU usage. 
- A couple of parameters should not be changed for consistency. These include:
    - `"atac.auto_detect_adapter": true`
    - `"atac.multimapping": 4`

2. `example.json` found [here](https://github.com/nicolerg/motrpac_atac_mop/blob/master/files/example.json) provides an example JSON file output by `make_json.sh`, specifically for a sample with a single biological replicate that had paired-end FASTQ files for four lanes.

### 4.2 Run the pipeline 
Actually running the pipeline is straightforward. However, the command is different depending on the environment in which you set up the pipeline. See details for every environment by following the corresponding link under the **Installation and tutorial heading** [here](https://github.com/ENCODE-DCC/atac-seq-pipeline).

For example, for the Google Cloud Platform implementation, use the following code:
```bash
PROJECT=[YOUR_PROJECT_NAME]
BUCKET=gs://[YOUR_BUCKET_NAME]/ENCSR356KRQ_subsampled
INPUT=examples/google/ENCSR356KRQ_subsampled.json

java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=google -Dbackend.providers.google.config.project=${PROJECT} -Dbackend.providers.google.config.root=${BUCKET} cromwell-34.jar run atac.wdl -i ${INPUT} -o workflow_opts/docker.json
```
For the “Local system with Conda” implementation, use the following code, where `${SRCDIR}` is defined as the path in which you installed the pipeline in Step 3.1:
```
source activate encode-atac-seq-pipeline # IMPORTANT!
INPUT=/path/to/json
OUTDIR=/path/to/output/directory
mkdir -p $OUTDIR
cd $OUTDIR
SRCDIR=/path/to/atac-seq-pipeline
java -jar -Dconfig.file=${SRCDIR}/backends/backend.conf ${SRCDIR}/cromwell-34.jar run ${SRCDIR}/atac.wdl -i ${INPUT}
```
A `cromwell-executions` directory containing all of the pipeline outputs is created in whatever directory from which you run the above command (specified in the “Local system with Conda” code as `$OUTDIR`), so choose wisely. One arbitrarily-named subdirectory for each JSON file (assuming the command is run in a loop for several samples) will be written in `cromwell-executions/atac`.  

### 4.3 Output directory structure
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

## 5. Compile important QC metrics
### 5.1 Review QC metrics in qc.json reports and figures generated by the pipeline 
This pipeline generates sample-level QC reports (qc.json). It is helpful to merge individual reports across samples in order to review and plot metrics across a group of samples. The ENCODE ATAC-seq pipeline GitHub provides a script to compile the QC metrics in the JSON QC reports for all samples into a single TSV file. Follow the instructions [here](https://github.com/ENCODE-DCC/atac-seq-pipeline/tree/master/utils/qc_jsons_to_tsv). 

At a minimum, the following files should be reviewed for QC purposes, assuming that the pipeline was run for multiple samples (e.g. all samples in a sequencing run). See **Table 4.3.2** for more details.     
- Insert size histogram
    - `call-ataqc/shard-0/execution/${SID}_L001_R1_001.trim.merged.inserts.hist_graph.pdf`
- TSS enrichment plots
    - `call-ataqc/shard-0/execution/${SID}_L001_R1_001.trim.merged_large_tss-enrich.png`
    - `call-ataqc/shard-0/execution/${SID}_L001_R1_001.trim.merged_tss-enrich.png`
- Sample-level QC reports
    - `call-qc_report/execution/${SID}.qc.html`
    - `call-qc_report/execution/${SID}.qc.json`

The generation of the merged QC TSV file as described above is useful for personal reference but not necessary since a larger number of sample-level JSON reports can be merged downstream.

**Table 5.1** provides definitions for a limited number of metrics included in the JSON QC reports. The full JSON report includes ~120 lines per sample; some lines are duplicates, and many metrics are irrelevant for running the pipeline with a single biological replicate. Therefore, this table includes a subset of relevant and important metrics.  

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

Until this issue is resolved, use the following scripts to quickly calculate the number of raw reads for each sample (alternatively, use the outputs of FastQC):
- [batch_fastq_read_count.sh](https://github.com/nicolerg/motrpac_atac_mop/blob/master/files/batch_fastq_read_count.sh)
    - `${fastq_read_count}`: path to fastq_read_count.sh 
    - `${indir}`: path to FASTQ files 
- [fastq_read_count.sh](https://github.com/nicolerg/motrpac_atac_mop/blob/master/files/fastq_read_count.sh)
    - Called by `batch_fastq_read_count.sh` to calculate the number of raw reads per sample in parallel

The number of raw reads and nodup_flagstat_qc:total should be used to calculate total read loss for each sample:
```
pct_reads_lost = (raw_reads - nodup_flagstat_qc:total) / raw_reads * 100 
```
### 5.2 Get counts of alignments across the genome
Use the script provided [here](https://github.com/nicolerg/motrpac_atac_mop/blob/master/files/chr_info.sh) to get a summary of the number and percent of reads in the final filtered BAM file that map to autosomal, sex, or mitochondrial chromosomes or contigs for each sample in `cromwell-executions/atac`.  Within the script, change `${base}` to point to `cromwell-executions/atac`; change `${outdir}` to reflect the desired output directory. After the script finishes running (it only takes a few seconds), sample-level samtools idxstats files as well as a summary file called `idxstats.summary.txt` can be found in the specified `${outdir}`. 

At a minimum, the values in the `pct_mt` column (percent of filtered reads that align to chrM) in `idxstats.summary.txt` should be reported for every sample. None of the metrics in `idxstats.summary.txt` are currently reported by the ENCODE ATAC-seq pipeline. 

## 6. Flag problematic samples
The following metrics are not strictly exclusion criteria for MoTrPAC samples, but samples should be flagged if any of these conditions are met. Some of these metrics reflect the ENCODE ATAC-seq data standards (available [here](https://www.encodeproject.org/atac-seq/#standards)). 

**Table 6.1 Criteria to flag problematic samples.**

| Description | In terms of Table 5.1 metrics | Comments |
|-------------|-------------------------------|----------|
|< 50 million filtered, non-duplicated, non-mitochondrial paired-end reads in the filtered BAM file (i.e. 25M pairs)|nodup_flagstat_qc:total < 50M|This is the most stringent criterion and may be relaxed|
|Alignment rate < 80%|flagstat_qc:mapped_pct < 80%||
|Fraction of reads in “overlap” peaks < 0.1|overlap_frip_qc:FRiP_rep1-pr < 0.1|This is more relaxed than the ENCODE recommendation|
|Number of peaks in “overlap” peak set < 80,000|overlap_reproducibility_qc:N1 < 80000|This is more relaxed than the ENCODE recommendation|
|A nucleosome-free region is not present|ataqc:Presence of NFR peak != ‘OK’|This should be enforced more strictly|
|A mononucleosome peak is not present|ataqc:Presence_of_Mono-Nuc_peak != ‘OK’|This should be enforced more strictly|
|TSS enrichment < 7|ataqc:TSS_enrichment < 7|It is more difficult to interpret rat TSS enrichment since there are no ENCODE standards for rat |












# ENCODE ATAC-seq pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.156534.svg)](https://doi.org/10.5281/zenodo.156534)[![CircleCI](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master)

## Introduction

This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq or DNase-seq data. The pipeline can be run on compute clusters with job submission engines or stand alone machines. It inherently makes uses of parallelized/distributed computing. Pipeline installation is also easy as most dependencies are automatically installed. The pipeline can be run end-to-end i.e. starting from raw FASTQ files all the way to peak calling and signal track generation; or can be started from intermediate stages as well (e.g. alignment files). The pipeline supports single-end or paired-end ATAC-seq or DNase-seq data (with or without replicates). The pipeline produces formatted HTML reports that include quality control measures specifically designed for ATAC-seq and DNase-seq data, analysis of reproducibility, stringent and relaxed thresholding of peaks, fold-enrichment and pvalue signal tracks.  The pipeline also supports detailed error reporting and easy [resuming](utils/qc_jsons_to_tsv/README.md) of runs. The pipeline has been tested on human, mouse and yeast ATAC-seq data and human and mouse DNase-seq data.

The ATAC-seq pipeline specification is also the official pipeline specification of the Encyclopedia of DNA Elements (ENCODE) consortium. The ATAC-seq pipeline protocol definition is [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing). Some parts of the ATAC-seq pipeline were developed in collaboration with Jason Buenrostro, Alicia Schep and Will Greenleaf at Stanford.

### Features

* **Flexibility**: Support for `docker`, `singularity` and `Conda`.
* **Portability**: Support for many cloud platforms (Google/DNAnexus) and cluster engines (SLURM/SGE/PBS).
* **Resumability**: [Resume](utils/qc_jsons_to_tsv/README.md) a failed workflow from where it left off.
* **User-friendly HTML report**: tabulated quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/cross-correlation measures).
  - Examples: [HTML](https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR889WQX/example_output/qc.html), [JSON](docs/example_output/v1.1.5/qc.json)
* **ATAqC**: Annotation-based analysis including TSS enrichment and comparison to Roadmap DNase.
* **Genomes**: Pre-built database for GRCh38, hg19, mm10, mm9 and additional support for custom genomes.

## Installation and tutorial

This pipeline supports many cloud platforms and cluster engines. It also supports `docker`, `singularity` and `Conda` to resolve complicated software dependencies for the pipeline. A tutorial-based instruction for each platform will be helpful to understand how to run pipelines. There are special instructions for two major Stanford HPC servers (SCG4 and Sherlock).

* Cloud platforms
  * Web interface
    * [DNAnexus Platform](docs/tutorial_dx_web.md)
  * CLI (command line interface)
    * [Google Cloud Platform](docs/tutorial_google.md)
    * [DNAnexus Platform](docs/tutorial_dx_cli.md)
* Stanford HPC servers (CLI)
  * [Stanford SCG4](docs/tutorial_scg.md)
  * [Stanford Sherlock 2.0](docs/tutorial_sherlock.md)
* Cluster engines (CLI)
  * [SLURM](docs/tutorial_slurm.md)
  * [Sun GridEngine (SGE/PBS)](docs/tutorial_sge.md)
* Local Linux computers (CLI)
  * [Local system with `singularity`](docs/tutorial_local_singularity.md)
  * [Local system with `docker`](docs/tutorial_local_docker.md)
  * [Local system with `Conda`](docs/tutorial_local_conda.md)
* Local Windows computers (CLI)
  * [Windows 10 Pro with `docker`](docs/tutorial_windows_docker.md)
  * [Windows 10 Pro/Home with `Conda`](docs/tutorial_windows_conda.md)

## Input JSON file

[Input JSON file specification](docs/input.md)

## Output directories

[Output directory specification](docs/output.md)

## Useful tools

There are some useful tools to post-process outputs of the pipeline.

### qc_jsons_to_tsv

[This tool](utils/qc_jsons_to_tsv/README.md) recursively finds and parses all `qc.json` (pipeline's [final output](docs/example_output/v1.1.5/qc.json)) found from a specified root directory. It generates a TSV file that has all quality metrics tabulated in rows for each experiment and replicate. This tool also estimates overall quality of a sample by [a criteria definition JSON file](utils/qc_jsons_to_tsv/criteria.default.json) which can be a good guideline for QC'ing experiments.

### resumer

[This tool](utils/resumer/README.md) parses a metadata JSON file from a previous failed workflow and generates a new input JSON file to start a pipeline from where it left off.

### ENCODE downloader

[This tool](https://github.com/kundajelab/ENCODE_downloader) downloads any type (FASTQ, BAM, PEAK, ...) of data from the ENCODE portal. It also generates a metadata JSON file per experiment which will be very useful to make an input JSON file for the pipeline.
