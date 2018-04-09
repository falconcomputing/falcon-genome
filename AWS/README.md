# Quickstart on Falcon Genome Image 

The Falcon Genome image includes the following:
- Software distributions for Falcon Genomics Solutions, installed at
`FALCON_DIR=/usr/local/falcon`
- Automated setup scripts, installed at `/root/setup-scripts`
- Dependent libraries and relevant environment setup at `/etc/profile.d`

There are three simple steps for using the Falcon Genome image on AWS
1. Request an EC2 instance with sufficient EBS volume;
2. Login to the EC2 instance and run the automated setup scripts following the
prompts;
3. Run the accelerated genomics data processing software using the `fcs-genome`
command from the ECS instance.

# System Setup

To use the automated setup script, you need to obtain the following information:
- path to the storage (can be either a storage device or network location)
- (optional) path to the reference genome

The automated system setup scripts will setup local EBS storage, and 
store system information into `/usr/local/falcon/fcs-genome.conf`, which will
be used by the `fcs-genome` software.

# Preparation

1. Temporary folder: some steps in `fcs-genome`, including `align`,
will write large files to a temporary folder. Please ensure there is enough
space inside the temporary folder. The recommended free space is 3~5x the
input data size
1. Falcon License: valid license need to be setup in the environment variable
`$LM_LICENSE_PATH`. It can either point to a file or a network position in the
form of `<port#>@<license-host>`. If the license file is not configured
correctly, the following error will be reported:
    ```
    [fcs-genome] ERROR: Cannot connect to the license server: -15
    [fcs-genome] ERROR: Please contact support@falcon-computing.com for details.
    ```
1. (Optional) Reference genome: to take full advantage of the FPGA acceleration
provided by the Falcon Genome image, the reference genome needs to be
preprocessed. This step is optional, and regular reference genome files (FASTA)
will still work without processing. The processed reference genome, on the
other hand, will also work for other softwares such as BWA, Picard, GATK,
etc.    
The preparation can be done with `$FALCON_DIR/prepare-ref.sh <path-to-fasta>`
1. (Optional) Known database indexing: GATK relies on VCF files for the known
variants such as dbSNP in its core processing. If index file for these VCF
databases do not exist, performance will be degraded. To generate these
indexes, either GATK or the `tabix` tool in the
[samtools distribution](http://www.htslib.org/download/) can be used.

# Start the Application
Please refer to the User Guide for `fcs-genome` for details on pipeline
execution.
