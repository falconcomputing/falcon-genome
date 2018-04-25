# Falcon Accelerated Genomics Pipeline User Guide

<!-- TOC depthFrom:2 depthTo:3 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Introduction](#introduction)
- [Quick Start on Public Clouds](#quick-start-on-public-clouds)
- [FAQ](#faq)
- [System Requirements and Installation](#system-requirements-and-installation)
	- [Software Prerequisites](#software-prerequisites)
	- [System Setup](#system-setup)
	- [Preparation](#preparation)
- [Synopsis](#synopsis)
	- [Common Options Among Methods](#common-options-among-methods)
	- [fcs-genome align](#fcs-genome-align)
	- [fcs-genome markdup](#fcs-genome-markdup)
	- [fcs-genome indel](#fcs-genome-indel)
	- [fcs-genome bqsr](#fcs-genome-bqsr)
	- [fcs-genome baserecal](#fcs-genome-baserecal)
	- [fcs-genome printreads](#fcs-genome-printreads)
	- [fcs-genome htc](#fcs-genome-htc)
	- [fcs-genome ug](#fcs-genome-ug)
	- [fcs-genome joint](#fcs-genome-joint)
	- [fcs-genome gatk](#fcs-genome-gatk)
- [Quick Start](#quick-start)
	- [Generating a Marked Duplicates BAM file from Paired-End FASTQ files](#generating-a-marked-duplicates-bam-file-from-paired-end-fastq-files)
	- [Performing Indel Re-alignment from a Marked Duplicates BAM file](#performing-indel-re-alignment-from-a-marked-duplicates-bam-file)
	- [Generating Base Quality Recalibration Report (BQSR) from a BAM file with known sites](#generating-base-quality-recalibration-report-bqsr-from-a-bam-file-with-known-sites)
	- [Generating Genomic VCF (gVCF) file from a BAM file with Haplotype Caller](#generating-genomic-vcf-gvcf-file-from-a-bam-file-with-haplotype-caller)
- [Tuning Configurations](#tuning-configurations)
	- [Reference Table for Configurations](#reference-table-for-configurations)

<!-- /TOC -->

## Introduction
The Falcon Accelerated Genomics Pipelines (FAGP) comprising the `fcs-genome` software allows for variant calling for both germline and somatic mutations based on the GATK Best Practices pipelines. The performance of the pipelines is significantly improved with Falcon's acceleration technologies.

Symmetric to the GATK Best Practices pipelines, the typical workflow starts with raw FASTQ sequence paired-end reads and proceeds to obtain a filtered set of variants that can be annotated for further analysis. The figure below depicts the flow of the germline variant calling pipeline. Beginning with paired-end FASTQ sequence files, the first step is to map the sequences to the reference. The resulting mapped BAM file is sorted, and duplicates are marked. This step performed using the command fcs-genome align, is equivalent to BWA-MEM, samtools sort and picard MarkDuplicates of the GATK Best Practices pipelines.

The second step is to recalibrate base quality score to account for biases caused by the sequencing machine. The Falcon pipeline command for this is `fcs-genome bqsr`. Its GATK equivalent first runs the GATK BaseRecalibrator, which produces a table of recalibrated reads, followed by GATK PrintReads which implements the table of recalibrated reads to produce a new, analysis-ready BAM file. The final step is germline variant calling, implementing the command fcs-genome htc which corresponds to GATK HaplotypeCaller.

![Falcon Workflow](resources/fcs-genome-workflow.jpeg)
Figure 1. Side-by-side analysis of the Falcon Accelerated Pipeline and the GATK Best Practices Pipeline: The middle panel indicates the general workflow starting with 1. Mapping the FASTQ sequences to the reference 2. Recalibrating base quality score and finally 3. Calling germline variants. The upper and lower panels illustrate the command-line implementation of the workflow using the Falcon Accelerated Pipeline and GATK Best Practices Pipeline respectively.

The table below shows which of the components of the GATK best practices have a Falcon accelerated counterpart and which ones are left in their original forms:

| Original Tool | Original Version | Command | Falcon Accelerated Command |
| --- | --- | --- | --- |  
| BWA | 0.7.13 | mem | align |
| samtools | 1.3 | view, sort | |
| picard | 1.141 | MarkDuplicates | |
| GATK | 3.8 | BaseRecalibrator | baserecal |
| | | PrintReads | printreads |
| | | HaplotypeCaller | htc |
| | | IndelRealigner | indel |
| | | UnifiedGenotyper | ug |
| | | CombinedGVCFs | joint |
| | | GenotypeGVCFs | |

This User Guide provides details on the setup of the Falcon Genome pipeline, command-line usage and a step-by-step example to run the variant calling pipeline.

## Quick Start on Public Clouds
- [AWS](aws/README.md)
- [Huawei Cloud](hwcloud/README.md)
- [Alibaba Cloud](aliyun/README.md)

## FAQ
[FAQ for Falcon Accelerated Genomics Pipelines.](FAQ.md)

## System Requirements and Installation
### Software Prerequisites
The software package of the Falcon Accelerated Genomics Pipelines is self-contained with required software. Please refer to the release notes inside each software distribution for each component and its version. The recommended operating system and required packages are listed as follows:
+ CentOS Linux 7.x
+ epel-release, boost, glog, gflags, java

### System Setup
+ The software for fcs-genome is installed in /usr/local/falcon
+ System information can be modified and is stored in /usr/local/fcs-genome.conf. Details on tuning configuration parameters is explained in a later section.
+ Export fcs-genome and other required tools to the PATH: source /usr/bin/falcon/setup.sh

### Preparation
+ Working folder: Paths to the reference genome and the input data are required parameters for the pipeline to run.  Setting up a working folder containing this data and allowing it to be readable is a mandatory step before the start of the pipeline.
+ Temporary folder: Most steps in the pipeline produce intermediate files that need to be stored at a temporary location. It must be ensured that this location has free disk space between 3-5X times the size of the input files. The location of the temporary folder can be modified in /usr/local/fcs-genome.conf.
+ Falcon License: A valid license needs to be setup in the environment variable $LM_LICENSE_PATH. If the license file is improperly configured, an error message is reported:
[fcs-genome] ERROR: Cannot connect to the license server: -15
[fcs-genome] ERROR: Please contact support@falcon-computing.com for details.
+ Obtaining the Reference and its index: The reference and its index can be downloaded from the Broad Institute website using the following FTP link:
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/  
To take full advantage of the FPGA acceleration provided by the Falcon Genome image, the reference genome needs to be preprocessed by running the script /usr/local/falcon/prepare-ref.sh <path-to-fasta>. This step is optional, and the regular reference genome files (FASTA) will still work without processing. The processed reference genome, on the other hand, will also work for other software such as BWA, Picard, GATK, etc.
+ Optional arguments: GATK relies on files with known variants in its processing. For example, known variant files including the 1000 Genome indel sites, the Mills indel sites and the dbSNP sites can be given as additional parameters for the pipeline steps. These can also be downloaded from the Broad Institute website.

## Synopsis
This section provides all the methods available in the fcs-genome command with their respective options settings.

```
fcs-genome align -r ref.fasta -1 input_1.fastq -2 input_2.fastq \
  -o aln.sorted.bam  --rg RG_ID --sp sample_id \
  --pl platform --lb library
fcs-genome markdup -i aln.sorted.bam -o aln.marked.bam
fcs-genome indel -r ref.fasta -i aln.sorted.bam -o indel.bam
fcs-genome bqsr -r ref.fasta -i indel.bam -o recal.bam
fcs-genome baserecal -r ref.fasta -i indel.bam -o recalibration_report.grp
fcs-genome printreads -r ref.fasta -b recalibration_report.grp -i indel.bam \
  -o recal.bam
fcs-genome htc -r ref.fasta -i recal.bam -o final.gvcf
fcs-genome joint -r ref.fasta -i final.gvcf.gz -o final.vcf
fcs-genome ug -r ref.fasta -i recal.bam -o final.vcf
```
For additional options, type in the command-line `fcs-genome [method]`.

In addition, the option `--extra-options` can be used to apply options in the original GATK that is not included in `fcs-genome`. For example:  
```
fcs-genome printreads -r ref.fasta -b recalibration_report.grp -i indel.bam \
  -o recal.bam --extra-options "-n 100000"
```

Please check the GATK documentation for all extra options available. The tables below show all options available in each method.  

(\*): Required

### Common Options Among Methods

| Option | Alternative | Argument | Description |
| --- | --- | --- | --- |
| -h | --help | | print help messages |
| -f | --force | | overwrite output files if they exist |
| -O | --extra-options | String(\*) | extra options in GATK for the command. Use " " to enclose the GATK command. Example "--TheOption arg" |

### fcs-genome align
Perform alignment using the Burrows-Wheeler Algorithm. It is the equivalent of bwa-mem. By default, mark duplicates are performed. If `--align-only` is set, no mark duplicate will be performed.

| Option | Alternative | Argument | Description |
| --- | --- | --- | --- |
| -r | --ref | String(\*) | reference genome path |
| -1 | --fastq1 | String(\*) | input pair-end Read 1 FASTQ file |
| -2 | --fastq2 | String(\*) | input pair-end Read 2 FASTQ file |
| -o | --output | String(\*) | output BAM file (if --align-only is set, the output will be a directory of BAM files) |
| -R | --rg | String(\*) | read group ID ('ID' in BAM header) |
| -S | --sp | String(\*) | sample ID ('SM' in BAM header) |
| -P | --pl | String(\*) | platform ID ('PL' in BAM header) |
| -L | --lb | String(\*) | library ID ('LB' in BAM header) |
| -l | --align-only | | skip mark duplicates |

### fcs-genome markdup
Takes a BAM file and mark duplicates the reads.

| Option | Alternative | Argument | Description |
| --- | --- | --- | --- |
| -i | --input | String(\*) | input BAM file |
| -o | --output | String(\*) | output BAM file |

### fcs-genome indel
Take a BAM file and perform indel re-alignment.

| Option | Alternative | Argument | Description |
| --- | --- | --- | --- |
| -r | --ref | String(\*) | reference genome path |
| -i | --input | String(\*) | input BAM file |
| -o | --output | String(\*) | output BAM file |
| -K | --known | String(\*) | known indels for realignment(VCF format). If more VCF are considered, add -K for each file |

### fcs-genome bqsr
Take a BAM file and perform Base Quality Score Recalibration. It can be performed within a region defined in the `--knownSites` option. If --bqsr is set, a report is generated.

| Option | Alternative | Argument | Description |
| --- | --- | --- | --- |
| -r | --ref | String(\*) | reference genome path |
| -b | --bqsr | String(\*) | output BQSR file (if left blank, no file will be produced) |
| -i | --input | String(\*) | input BAM file |
| -o | --output | String(\*) | output BAM file |
| -K | --knownSites | String(\*) | known indels for realignment (VCF format). If more VCF are considered, add -K for each file |

### fcs-genome baserecal
Take a BAM file and generate a Base Quality Score Recalibration.  

| Option | Alternative | Argument | Description |
| --- | --- | --- | --- |
| -r | --ref | String(\*) | reference genome path |
| -i | --input | String(\*) | input BAM file |
| -o | --output | String(\*) | output BQSR file |
| -K | --knownSites | String(\*) | known indels for realignment (VCF format). If more VCF are considered, add -K for each file |

### fcs-genome printreads
Take a BAM file and filter reads according to some settings defined in `--extra-options`.

| Option | Alternative | Argument | Description |
| --- | --- | --- | --- |
| -r | --ref | String(\*) | reference genome path |
| -b | --bqsr | String(\*) | input BQSR file |
| -i | --input | String(\*) | input BAM file or directory |
| -o | --output | String(\*) | output BAM files |

### fcs-genome htc
Take a BAM file and generate a gVCF file by default.  If --produce-vcf is set, a VCF file is generated instead of gVCF.

| Option | Alternative | Argument | Description |
| --- | --- | --- | --- |
| -r | --ref | String(\*) | reference genome path |
| -i | --input | String(\*) | input BAM file or directory |
| -o | --output | String(\*) | output gVCF/VCF file (if --skip-concat is set the output will be a directory of gVCF files) |
| -v | --produce-vcf | | produce VCF files from HaplotypeCaller instead of gVCF |
| -s | --skip-concat | | (deprecated) produce a set of gVCF/VCF files instead of one |

### fcs-genome ug
This method is the equivalent of UnifiedGenotype in GATK. It takes a BAM file as an input and generates a VCF file.  It accepts options from GATK through `--extra-option`

| Option | Alternative | Argument | Description |
| --- | --- | --- | --- |
| -r | --ref | String(\*) | reference genome path |
| -i | --input | String(\*) | input BAM file or directory |
| -o | --output | String(\*) | output a compressed VCF file |
| -s | --skip-concat | String(\*) | produce a set of vcf files instead of one |

### fcs-genome joint
This method performs a joint variant calling from a set of VCF files.

| Option | Alternative | Argument | Description |
| --- | --- | --- | --- |
| -r | --ref | String(\*) | reference genome path |
| -i | --input-dir | String(\*) | input dir containing compressed gVCF files |
| -o | --output | String(\*) | output compressed gVCF files |
| -c | --combine-only | | combine GVCFs only and skip genotyping |
| -g | --skip-combine | | (deprecated) perform genotype gVCFs only and skip combine gVCF |

### fcs-genome gatk
This method emulates the original GATK command.  Please refer the GATK documentation for additional details.

## Quick Start
The examples below were written in BASH script and quickly tested using an instance of 16-cores (Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz, 2 threads per core).   Each example below can be saved in a file and be submitted to the server as follows:
```
chmod a+x myscript.sh ; nohup ./myscript.sh &
```
For illustration purposes, the FASTQ files (small_1.fastq.gz and small_2.fastq.gz) used in the examples below contain 10K paired-end reads. They can be generated easily from any paired-end reads FASTQ files using the following Linux commands:  
```
zcat originalFASTQ_R1.fastq.gz | head -n 40000 > small_1.fastq ; gzip small_1.fastq
zcat originalFASTQ_R2.fastq.gz | head -n 40000 > small_2.fastq ; gzip small_2.fastq
```
In FASTQ format, each DNA read consists of 4 lines. Therefore, to get 10,000 DNA reads, 40,000 lines need to be extracted from the original FASTQ file.
For more exhaustive test, the platinum pedigree samples (NA12878, NA12891 and NA12892) can be used as examples. They can be downloaded from http://www.internationalgenome.org/data-portal/sample/.  Alternatively, Illumina BaseSpace (account required) provides Public Data sequenced with the most recent technology.

### Generating a Marked Duplicates BAM file from Paired-End FASTQ files
fcs-genome align performs alignment to the reference, sorts, marks duplicates, and save the mapped reads in a BAM file. If --align-only is set, no marking duplicates is performed. The BASH script below illustrates the usage of the align method:
```
SAMPLE_ID="small"
R1=${SAMPLE_ID}_1.fastq.gz
R2=${SAMPLE_ID}_2.fastq.gz"
REF="/local/ref/human_g1k_v37.fasta
BAMFILE=${SAMPLE_ID}_marked_sorted.bam
RG_ID="H0BA0ADXX"
PLATFORM="Illumina"
LIB="RD001"

fcs-genome align \
  -r $REF -1 $R1 -2 $R2 \
  -o ${BAMFILE} \
  --rg $RG_ID --sp ${SAMPLE_ID} \
  --pl ${PLATFORM} --lb ${LIB}
```
For 10K paired-reads contained in the FASTQ files, it took 13 seconds for alignment and 1 second for marking duplicates.  The BAM file is generated with its respective index.

### Performing Indel Re-alignment from a Marked Duplicates BAM file
Once the alignment is completed,  indel-realignment is perfomed.  The BASH script below demonstrates the usage of the indel method:
```
REF="/local/ref/human_g1k_v37.fasta
SAMPLE_ID="small"
BAM_INPUT=${SAMPLE_ID}_marked_sorted.bam
BAM_OUTPUT=${SAMPLE_ID}_marked_sorted_indel_realign

fcs-genome indel \
  -r $REF \
  -i ${BAM_INPUT} \
  -o ${BAM_OUTPUT}
```
A folder called ${SAMPLE_ID}_marked_sorted_indel_realign/ is created with a set of BAM and bai files with indels re-aligned. It takes 100 seconds to perform.
Performing Base Quality Score Recalibration (BQSR) from BAM file with pre-defined known sites
fcs-genome bqsr performs GATK's Base Quality Score Recalibration and Print Reads in a single command. Per-base quality scores produced by the sequencing machine are checked for errors and corrected. The recalibrated reads are written into a folder that contains a BAM files set. During the process, a recalibration report is generated. The script below illustrates the usage of bqsr method:
```
REF="/local/ref/human_g1k_v37.fasta
ThousandGen="/local/ref/1000G_phase1.indels.b37.vcf"
Mills="/local/ref/Mills_and_1000G_gold_standard.indels.b37.vcf"
SNP="/local/ref/dbsnp_138.b37.vcf"
SAMPLE_ID="small"
BAM_INPUT=${SAMPLE_ID}_marked_sorted_indel_realign
BAM_OUTPUT=${SAMPLE_ID}_recalibrated

fcs-genome bqsr \
  -r $REF \
  -i ${BAM_INPUT} -o ${BAM_OUTPUT} \
  -b recalibration_report.grp \
  -K $ThousandGen -K $Mills -K $SNP"
```
For this example, it took 1203 seconds to complete.

### Generating Base Quality Recalibration Report (BQSR) from a BAM file with known sites
In this example, the BQSR analysis was performed using as an input a folder that contained BAM files and their respective bai files.  A base recalibration report is generated.
```
REF="/local/ref/human_g1k_v37.fasta
SAMPLE_ID="small"
BAM_INPUT=${SAMPLE_ID}_marked_sorted_indel_realign
ThousandGen="/local/ref/1000G_phase1.indels.b37.vcf"
Mills="/local/ref/Mills_and_1000G_gold_standard.indels.b37.vcf"
SNP="/local/ref/dbsnp_138.b37.vcf"

fcs-genome baserecal \
  -r $REF \
  -i ${BAM_INPUT} -o recalibration_report.grp \
  -K $ThousandGen -K $Mills -K $SNP"
```
The command also works with a single BAM file.  It takes around 1177 seconds to complete.

### Generating Genomic VCF (gVCF) file from a BAM file with Haplotype Caller
fcs-genome htc performs germline variant calling using the input BAM file with default output format as gVCF. if --produce-vcf is set, a VCF file is produced.
```
SAMPLE_ID="small"
REF="/local/ref/human_g1k_v37.fasta
BAM_INPUT=${SAMPLE_ID}_recalibrated.bam
OutputVCF=${SAMPLE_ID}_final.gvcf

fcs-genome htc \
  -r ${REF} \
  -i ${BAM_INPUT} \
  -o ${OutputVCF}
```
For this example, it takes 415 seconds to complete. The htc option accepts multiple BAM files as input.

## Tuning Configurations
Configurations can be tuned to define the settings for each command-line option during the run. The default configuration settings are stored in /usr/local/fcs-genome.conf. If a file with the same name `fcs-genome.conf` is presented in the present directory, its values will be used to overwrite the default values. In addition, environmental variables can be used to overwrite both default configurations and the configurations in `fcs-genome.conf` in the present directory.
An example of the configuration settings for the germline variant calling pipeline is as below:
```
temp_dir = /local/temp
gatk.ncontigs = 32
gatk.nprocs = 16
gatk.nct = 1
gatk.memory = 8
```
The key `temp_dir` specifies the system folder to store temporary files. Some steps in `fcs-genome`, including `align`, will write large files to a temporary folder. Please ensure this configuration is set to a location with enough space. The recommended free space is 3~5x the input data size.
The GATK steps, such as BaseRecalibratior, PrintReads and HaplotypeCaller, are run in parallel. By default, 32 total processes will be used for each GATK step. To change the default number, the key `gatk.ncontigs` be set. The configuration key `gatk.nprocs` is used to specify the number of concurrent processes in each step. `gatk.memory` the memory consumed by each process. Ideally, `gatk.nprocs` be less than or equal to the total number of CPU cores, and the product of gatk.nprocs and gatk.memory would be less than or equal to the total memory. The number of concurrent process number and memory per process can be changed to individual steps with the following format: [step-name].nprocs, [step-name].memory

### Reference Table for Configurations
| Configuration key | Argument Type | Default Value | Description |
| --- | --- | --- | --- |
| bwa.verbose | int | 0 | verbose level of bwa output |
| bwa.nt | int | -1 | number of threads for bwa, default is set to use all available threads in the system |
| bwa.num_batches_per_part | int | 20 | max num records in each BAM file |
| bwa.use_fpga | bool | true | option to enable FPGA for bwa-mem |
| bwa.use_sort | bool | true | enable sorting in bwa-mem |
| bwa.enforce_order | bool | true | enforce strict sorting ordering |
| bwa.fpga.bit_path | string | "" | path to FPGA bitstream for bwa |
| bwa.scaleout_mode | bool | | enable scale-out mode for bwa |
| markdup.max_files | int | 4096 | max opened files in markdup |
| markdup.nt | int | 16 | thread num in markdup |
| markdup.overflow-list-size | int | 2000000 | Overflow list size in markdup |
| gatk.scalout_mode | bool | | enable scale-out mode for gatk |
| gatk.intv.path | string | "" | default path to existing contig intervals |
| gatk.ncontigs | int | 32 | default contig partition num in GATK steps |
| gatk.nprocs | int | | default process num in all GATK steps, set to cpu num or gatk.ncontics whichever is the lesser value |
| gatk.nct | int | 1 | default thread number in GATK steps |
| gatk.memory | int | 8 | default heap memory in GATK steps |
| gatk.skip_pseudo_chr | bool | | skip pseudo chromosome intervals |
| gatk.bqsr.nprocs | int | | default process num in GATK BaseRecalibrator |
| gatk.bqsr.nct | int | | default thread num in GATK BaseRecalibrator |
| gatk.bqsr.memory | int | | default heap memory in GATK BaseRecalibrator |
| gatk.pr.nprocs | int | | default process num in GATK PrintReads |
| gatk.pr.nct | int | | default thread num in GATK PrintReads |
| gatk.pr.memory | int | | default heap memory in GATK PrintReads |
| gatk.htc.nprocs | int | | default process num in GATK HaplotypeCaller |
| gatk.htc.nct | int | | default thread num in GATK HaplotypeCaller |
| gatk.htc.memory | int | | default heap memory in GATK HaplotypeCaller |
| gatk.indel.nprocs | int | | default process num in GATK IndelRealigner |
| gatk.indel.memory | int | | default heap memory in GATK IndelRealigner |
| gatk.ug.nprocs | int | | default process num in GATK UnifiedGenotyper |
| gatk.ug.nt | int | | default thread num in GATK UnifiedGenotyper |
| gatk.ug.memory | int | | default heap memory in GATK UnifiedGenotyper |
| gatk.rtc.nt | int | 16 | default thread num in GATK UnifiedGenotyper |
| gatk.rtc.memory | int | 48 | default heap memory in GATK UnifiedGenotyper |
| gatk.joint.ncontigs | int | 32 | default contig partition num in joint genotyping |
| gatk.combine.nprocs | int | 16 | default process num in GATK CombineGVCFs |
| gatk.genotype.nprocs | int | 32 | default process num in GATK GenotypeGVCFs |
| gatk.genotype.memory | int | 4 | default heap memory in GATK GenotypeGVCFs |
