# fq2count
- a pipeline from bulk RNA-seq data to expresssion count matrix 
- 参考学习自：https://github.com/dongxuemin666/RNA-combine 加以自己的理解进一步整理
- 支持基于Bulk RNAseq fastq数据的质控、比对、定量流程的集成化处理，具体使用方式参考如下。
![pipeline.png](https://s6.jpg.cm/2022/02/20/LfLExt.png)
# 1、准备分析环境

## 1.1 linux文件夹环境

```shell
git clone https://github.com/lishensuo/fq2count.git
# gitee备用(快)： git clone https://gitee.com/li-shensuo/fq2count.git
tar -xvf fq2count.tar
work_path=/home/ssli/fq2count
cd $work_path

#为script脚本增加执行权限
##前一个为批量比对的脚本文件；后四个为去除rRNA相关的脚本文件
chmod u+x ${work_path}/scripts/BulkRNAseq.sh
chmod u+x ${work_path}/scripts/indexdb_rna
chmod u+x ${work_path}/scripts/merge-paired-reads.sh
chmod u+x ${work_path}/scripts/sortmerna
chmod u+x ${work_path}/scripts/unmerge-paired-reads.sh

#新建将会用到的文件夹
mkdir ${work_path}/0.rawfq
mkdir ${work_path}/1.rm_rrna
mkdir ${work_path}/2.trim
mkdir ${work_path}/3.align
mkdir ${work_path}/4.count]
mkdir ${work_path}/rrna_index
```



## 1.2 conda环境

```shell
cat requirement.txt
# bioconda::samtools==1.14
# bioconda::subread==2.0.1
# bioconda::refgenie==0.12.1
# bioconda::sra-tools==2.11.0
# bioconda::hisat2==2.2.1

conda create -n fq2count -y
conda activate fq2count
conda install -c conda-forge mamba -y
conda install --file=requirement.txt -y
```



# 2、准备数据

## 2.1 参考基因组gtf文件

- 通过refgenie下载，这里以hg38版本为例

```shell
#第一次使用refgenie需要运行下面两行代码
mkdir ~/refgenie
refgenie init -c ~/refgenie/genome_config.yaml

#下载hg38 gtf
refgenie pull hg38/gencode_gtf -c ~/refgenie/genome_config.yaml
echo $(refgenie seek hg38/gencode_gtf -c ~/refgenie/genome_config.yaml)
```

## 2.2 比对软件(hisat2)索引文件

- 通过refgenie下载，这里以hg38版本为例

```shell
refgenie pull hg38/hisat2_index -c ~/refgenie/genome_config.yaml
echo $(refgenie seek hg38/hisat2_index -c ~/refgenie/genome_config.yaml)
```

## 2.3 构建rRNA索引文件（optinonnal）

- 如果考虑去除fastq文件里的rRNA reads，需要运行这一步

```shell
#使用脚本文件scripts/indexdb_rna创建核糖体索引
cd ${work_path}/rrna_index
${work_path}/scripts/indexdb_rna --ref ${work_path}/scripts/rRNA_databases/silva-bac-16s-id90.fasta,./silva-bac-16s-db:\
${work_path}/scripts/rRNA_databases/silva-bac-23s-id98.fasta,./silva-bac-23s-db:\
${work_path}/scripts/rRNA_databases/silva-arc-16s-id95.fasta,./silva-arc-16s-db:\
${work_path}/scripts/rRNA_databases/silva-arc-23s-id98.fasta,./silva-arc-23s-db:\
${work_path}/scripts/rRNA_databases/silva-euk-18s-id95.fasta,./silva-euk-18s-db:\
${work_path}/scripts/rRNA_databases/silva-euk-28s-id98.fasta,./silva-euk-28s:\
${work_path}/scripts/rRNA_databases/rfam-5s-database-id98.fasta,./rfam-5s-db:\
${work_path}/scripts/rRNA_databases/rfam-5.8s-database-id98.fasta,./rfam-5.8s-db
```

## 2.4 测序数据

```shell
#如下为示例分析文件
cat SRR_list.txt
# SRR12720999
# SRR12721000
# SRR12721001
SRR_file=SRR_list.txt
cd ${work_path}/0.rawfq
for srr in $(cat $SRR_file)
do
#下载sra文件
prefetch -p -X 35G ${srr} -O .
#拆分为fastq文件
fasterq-dump --split-files ${srr} -p -O  ./
rm -rf ${srr}
done

#注意是未压缩的fastq文件，主要是考虑需要去除rRNA的情况
ls
-rw-r--r-- 1 ssli  7.5G Feb 20 10:34 SRR12720999_1.fastq
-rw-r--r-- 1 ssli  7.5G Feb 20 10:34 SRR12720999_2.fastq
-rw-r--r-- 1 ssli  7.4G Feb 20 10:45 SRR12721000_1.fastq
-rw-r--r-- 1 ssli  7.4G Feb 20 10:45 SRR12721000_2.fastq
-rw-r--r-- 1 ssli  8.0G Feb 20 10:53 SRR12721001_1.fastq
-rw-r--r-- 1 ssli  8.0G Feb 20 10:53 SRR12721001_2.fastq
```



# 3、批量比对

- 指定相关参数

```shell
work_path=/home/ssli/rnaseq/fq2count
#交代是否去除核糖体rRNA：0表示不去除，1表示去除(比较耗时)
Trim_rRNA=0
#交代测序读长（根据实际fastq数据）
read_length=100
#交代线程数
threads=20
#交代参考文件的路径（根据需要选择合适的版本）
hisat_index=$(refgenie seek hg38/hisat2_index -c ~/refgenie/genome_config.yaml)
gtf_path=$(refgenie seek hg38/gencode_gtf -c ~/refgenie/genome_config.yaml)
```

- 批量比对

```shell

${work_path}/scripts/BulkRNAseq.sh $work_path $Trim_rRNA $read_length \
$threads $hisat_index $gtf_path

head ${work_path}/4.count/expression_matrix.txt
# Geneid ../3.align/SRR12720999_nsorted.bam ../3.align/SRR12721000_nsorted.bam ../3.align/SRR12721001_nsorted.bam
# ENSG00000223972.5 0 0 0
# ENSG00000227232.5 24 15 14
# ENSG00000278267.1 12 5 8
# ENSG00000243485.5 0 0 1
# ENSG00000284332.1 0 0 0
# ENSG00000237613.2 0 0 0
# ENSG00000268020.3 0 0 0
# ENSG00000240361.2 0 0 0
# ENSG00000186092.6 0 0 0
```
