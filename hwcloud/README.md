# 峰科加速基因分析镜像设置说明

峰科基因镜像包括如下组件：
- 峰科基因分析软件*fcs-genome*，安装在`FALCON_DIR=/usr/local/falcon`
- 自动配置脚本，安装在`/root/setup.sh`
- 峰科软件所需系统组件，以及FPGA调用所需系统配置。

# 准备工作

1. 申请一台FPGA云服务器
1. （可选）挂载云盘，峰科加速基因软件对本地存储有较高存取吞吐率要求，建议使用超高IO本地云盘存储输入输出和临时文件以最大化分析性能。
1. 确认临时文件存储目录，并确保临时目录有足够可用空间，一般情况下，峰科加速软件运行过程中会产生临时文件，大小为输入*FASTQ.GZ*文件的3到5倍。
1. 配置峰科软件License的环境变量`$LM_LICENSE_PATH`。此环境变量既可以指向文件，也可以是远程地址：`<port>@<hostname>`。默认情况下本镜像已包含已有License配置。如果License配置错误或者License已过期，将会出现如下提示：
    ```
    [fcs-genome] ERROR: Cannot connect to the license server: -15
    [fcs-genome] ERROR: Please contact support@falcon-computing.com for details.
    ```
    如需获取最新License请联系峰科销售，或者发邮件至support@falcon-computing.com

# 镜像配置
镜像配置可通过运行自动配置脚本完成，此脚本包含峰科基因软件所需临时目录配置以及参考基因预处理配置。登陆云服务器后，按如下命令运行配置脚本：
```
> cd /root/
> ./setup.sh
```

配置脚本会执行如下几步：
1. 配置临时文件夹。如果用户已配置好本地存储，此步可直接提供存储路径；如果还未配置本地存储，用户需提供云盘位置（如/dev/vdb），此脚本会初始化并加载本地云盘。脚本提示如下：
    ```
    Setting up working dir...
    If you already have the working directory ready please enter
    the dir path, otherwise, please enter 'continue' or 'c':
    ```
    如果用户已提前配置好本地存储，此处直接输入本地存储地址，按脚本提示如下：
    ```
    If you already have the working directory ready please enter
    the dir path, otherwise, please enter 'continue' or 'c': /local
    ```
    如果用户需要配置本地存储，此处输入`continue`或者`c`，之后按照提示输入本地存储挂载位置和工作目录地址，以下假设挂载路径工作文件夹为`/local`：
    ```
    If you already have the working directory ready please enter
    the dir path, otherwise, please enter 'continue' or 'c': c

    Please enter the storage device: /dev/vdb
    Please enter the path to the work dir: /local
    ```
1. 预处理参考基因组。为了性能最大化，峰科软件提供参考基因组预处理。预处理为可选步骤，不影响正常使用，普通参考基因组未经预处理也可以使用，处理过的参考基因组也与标准文件匹配，不影响其他基因分析软件如*BWA*, *GATK*, *Picard*等的使用。脚本提示如下：
    ```
    Preparing reference genome...
    Please enter the path to the reference genome, or leave it blank to skip this step:
    ```
    如需跳过此步，直接输入回车。如需执行预处理步骤，此处应输入参考基因组位置，例如：`/local/ref/human_g1k_v37.fasta`。根据参考基因数据特性，预处理可能需要十几分钟到一两个小时。
1. 配置完毕会出现如下信息：
    ```
    Configuration is successful.
    ```
    配置成功后，查看`/usr/local/falcon/fcs-genome.conf`配置文件中将会有如下条目：
    ```
    temp_dir = /local/temp
    ```

# 运行软件
本镜像提供人重基因组分析示例脚本，位置在`$FALCON_DIR/example-wgs-germline.sh`。此示例脚本完成人重基因组分析中基因比对，去冗余，碱基质量校对和突变检测几个步骤。使用示例脚本前，需要配置脚本内如下几个参数，以下按之前提供的示例参数配置：
```
local_dir=/local
fastq_dir=/local/fastq
ref_dir=/local/ref
ref_genome=$ref_dir/human_g1k_v37.fasta
db138_SNPs=$ref_dir/dbsnp_138.b37.vcf
```

脚本使用方式为：
```
./example-wgs-germline.sh sample_id
```
其中sample_id为基因样本名称，给定之后，示例脚本将会在`$fastq_dir`中查找如下输入文件：
```
$fastq_dir/${sample_id}_1.fastq.gz
$fastq_dir/${sample_id}_2.fastq.gz
```
如果输入文件命名不符合此结构，示例脚本无法正常执行。示例脚本需要`$local_dir`内有输入文件大小5倍以上的空间，输出文件为：
- 去冗余BAM文件：`${sample_id}.bam`
- 重校对BAM文件：`${sample_id}.recal.bam`
- 突变检测结果VCF文件：`${sample_id}.vcf`

峰科基因分析软件具体使用方法请参考使用文档（英文）。
