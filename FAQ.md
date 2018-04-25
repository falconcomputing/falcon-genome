# FAQ for Falcon Accelerated Genomics Pipelines
## Navigation
<!-- TOC depthFrom:3 depthTo:3 withLinks:1 updateOnSave:1 orderedList:0 -->

- [1. What are the accelerated steps available in the Falcon Accelerated Genomics Pipelines?](#1-what-are-the-accelerated-steps-available-in-the-falcon-accelerated-genomics-pipelines)
- [2. `fcs-genome align` fails, and `bwa-flow` log reports `[E::bwa_idx_load_from_disk]` error.](#2-fcs-genome-align-fails-and-bwa-flow-log-reports-ebwaidxloadfromdisk-error)
- [3 `fcs-genome bqsr` fails with the missing fasta index file.](#3-fcs-genome-bqsr-fails-with-the-missing-fasta-index-file)
- [4. Falcon accelerated GATK steps fails will `Killed` message.](#4-falcon-accelerated-gatk-steps-fails-will-killed-message)

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
