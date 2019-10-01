# FAQ for Falcon Accelerated Genomics Pipelines
## Navigation
<!-- TOC depthFrom:3 depthTo:3 withLinks:1 updateOnSave:1 orderedList:0 -->

- [1. What are the accelerated steps available in the Falcon Accelerated Genomics Pipelines?](#1-what-are-the-accelerated-steps-available-in-the-falcon-accelerated-genomics-pipelines)
- [2. `fcs-genome align` fails, and `bwa-flow` log reports `[E::bwa_idx_load_from_disk]` error.](#2-fcs-genome-align-fails-and-bwa-flow-log-reports-ebwaidxloadfromdisk-error)
- [3 `fcs-genome bqsr` fails with the missing fasta index file.](#3-fcs-genome-bqsr-fails-with-the-missing-fasta-index-file)
- [4. Falcon accelerated GATK steps fails will `Killed` message.](#4-falcon-accelerated-gatk-steps-fails-will-killed-message)
- [5. `fcs-genome bqsr`, `fcs-genome ir`, `fcs-genome baserecals` slow or stuck issue.](#5-fcs-genome-bqsr-fcs-genome-ir-fcs-genome-baserecals-slow-or-stuck-issue)
- [6. `fcs-genome mutect2`, contiguous chromosome issue.](#5-fcs-genome-mutect2-contiguous-chromosome-issue)


<!-- /TOC -->

### 1. What are the accelerated steps available in the Falcon Accelerated Genomics Pipelines?
#### Answer
Please see the table below for a list of accelerated steps:

| Original Tool | Original Version | Command    | Falcon Accelerated Command |
| ------------- | ---------------- | ---------- | -------------------------- |  
| BWA           | 0.7.13           | mem        | align                      |
| SAMTOOLS      | 1.3              | view, sort | align                      |
| PICARD        | 1.141            | MarkDuplicates  | align                 |
| GATK | 3.8 | BaseRecalibrator | baserecal |
|      |     | PrintReads       | printreads |
|      |     | HaplotypeCaller  | htc |
|      |     | IndelRealigner   | indel |
|      |     | UnifiedGenotyper | ug |
|      |     | CombinedGVCFs    | joint |
|      |     | GenotypeGVCFs    | joint |

### 2. `fcs-genome align` fails, and `bwa-flow` log reports `[E::bwa_idx_load_from_disk]` error.
Error message:
```
[bwa-flow]: [E::bwa_idx_load_from_disk] fail to locate the index files
```
#### Answer
This error is caused by missing index of the reference genome. Assuming the reference genome filename is *ref.fa*, then `bwa-flow` will look for index file at *ref.fa.fai*. In addition, the following files related to *ref.fa* will be required:
- *ref.fa.amb*
- *ref.fa.ann*
- *ref.fa.bwt*
- *ref.fa.fai*
- *ref.fa.pac*
- *ref.fa.sa*

These files can be generated using */usr/local/falcon/prepare-ref.sh* helper script provided as part of the solution.

### 3 `fcs-genome bqsr` fails with the missing fasta index file.
Error Message:
```
##### ERROR ------------------------------------------------------------------------------------------
##### ERROR A USER ERROR has occurred (version 3.8-falcon-v0.4.1-1-g91287ae):
##### ERROR
##### ERROR This means that one or more arguments or inputs in your command are incorrect.
##### ERROR The error message below tells you what is the problem.
##### ERROR
##### ERROR If the problem is an invalid argument, please check the online documentation guide
##### ERROR (or rerun your command with --help) to view allowable command-line arguments for this tool.
##### ERROR
##### ERROR Visit our website and forum for extensive documentation and answers to
##### ERROR commonly asked questions https://software.broadinstitute.org/gatk
##### ERROR
##### ERROR Please do NOT post this error to the GATK forum unless you have really tried to fix it yourself.
##### ERROR
##### ERROR MESSAGE: Fasta index file /local/ref/human_g1k_v37.fasta.fai for reference /local/ref/human_g1k_v37.fasta does not exist. Please see https://software.broadinstitute.org/gatk/documentation/article?id=1601 for help creating it.
##### ERROR ------------------------------------------------------------------------------------------
```

#### Answer
The error is the same as the previous case, which is caused by missing reference index files. According to GATK documentation shown in the [link](https://software.broadinstitute.org/gatk/documentation/article?id=1601) in the error messsage, GATK requires `.dict` and `.fai` files of the reference genome.

### 4. Falcon accelerated GATK steps fails will `Killed` message.
Error message:
```
[2018-04-24 12:01:20 fcs-genome] INFO: Start doing Haplotype Caller
sh: line 1: 11906 Killed java -d64 -Xmx4g -jar /usr/local/falcon/bin/../tools/package/GenomeAnalysisTK.jar -T HaplotypeCaller -R ...
```
#### Answer
The error message indicates that the JAVA process running GATK was killed by the OS. More information can be obtained by running `dmesg`.
For the majority of cases, the issue is due to memory issue. `fcs-genome` will automatic decide the memory usage for each GATK step, but due to the unpredictability of JVM the prediction can be wrong sometime. To fix the issue, the first thing to try is specifying in the *fcs-genome.conf* configuration file, lowering the resource consumption of `fcs-genome`.

For example, if the machine has 32 core and 128Gb memory, we can modify the configuration file with the following two lines:
```
gatk.nprocs = 16
gatk.memory = 4
```

### 5. `fcs-genome bqsr`, `fcs-genome ir`, `fcs-genome baserecals` slow or stuck issue.
Affected commands all have `--known-sites` argument enabled. GATK may sometimes reports error message similar to this:
```
WARN 20:46:58,630 RMDTrackBuilder - Index file dbsnp_138.hg19.vcf.idx is out of date (index older than input file), falling back to an in-memory index
```
#### Answer
The reason is that GATK does not recognize valid index file for the VCF database, either because the index is missing or because it's older than the VCF file. This will result in GATK rebuilding the index in-memory, which causes the process to be very slow. The fix is to build the index before a run, or update the index file by `touch` if it's already built.
To build an index if it's not present, simply use the original GATK to run a command:
```
fcs-genome gatk -T BaseRecalibrator \
    -R hg19.fa \
    -I any-input.bam \
    -knownSites dbsnp_138.hg19.vcf \
    -o output.rpt
```
The `*.idx` file will be automatically build in the same directory of `dbsnp_138.hg19.vcf`.

### 5. `fcs-genome mutect2` contiguous chromosome issue.
Affected command displays ERROR Message: 
```
The chromosome block chr19 is not contiguous, consider running with -a
```
#### Answer
This ERROR message indicates that the data contains unsorted inputs or entries that were not generated sequentially. Current pipeline assumes that parts BAM files were generated following the coordinates from the reference fasta file. If the sample was sequenced using a capture kit, then the BED file that covers the regions defined by the capture should be set in every fcs-genome tool used in the pipeline where the interval option (-L) is available such as BQSR, HTC, mutect2, etc. It is recommended that BED file should be sorted by chromosome and position before execution:

sort -k1,1V -k2,2n IntervalFile.bed > IntervalFile_sorted.bed




