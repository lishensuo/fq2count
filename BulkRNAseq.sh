#!/bin/bash

#交代主工作路径
work_path=$1
#交代是否去除核糖体rRNA(1: Yes; 0: No)
Trim_rRNA=$2
#交代测序读长
read_length=$3
#交代线程数
threads=$4
#交代参考文件的路径
hisat_index=$5
gtf_path=$6


# work_path=/home/ssli/fq2count
# Trim_rRNA=0
# read_length=100
# threads=20
# hisat_index=$(refgenie seek hg38/hisat2_index -c ~/refgenie/genome_config.yaml)
# gtf_path=$(refgenie seek hg38/gencode_gtf -c ~/refgenie/genome_config.yaml)

# ${work_path}/scripts/BulkRNAseq.sh $work_path $Trim_rRNA $read_length $threads $hisat_index $gtf_path

# ${work_path}/scripts/test.sh $work_path $Trim_rRNA $read_length \
# $threads $hisat_index $gtf_path



#先检查一遍fastq数据
cd ${work_path}/0.rawfq
flag=1
files=$(ls *_1.fastq 2> /dev/null | wc -l)
if [ "$files" = "0" ] ; then echo "file should be like ***_1.fastq and ***_2.fastq";  fi

#开始批量比对
for file in $(ls | grep _1.fastq$)
do
pre_name=${file%_1.fastq}
fq1="${pre_name}_1.fastq"
fq2="${pre_name}_2.fastq"

# step 1
# remove ribosome genes
cd ${work_path}/1.rm_rrna 
if [ $Trim_rRNA == 1  ]; then
${work_path}/scripts/merge-paired-reads.sh \
${work_path}/0.rawfq/$fq1 ${work_path}/0.rawfq/$fq2 \
${pre_name}_merge.fastq
${work_path}/scripts/sortmerna \
--ref ${work_path}/scripts/rRNA_databases/silva-bac-16s-id90.fasta,${work_path}/rrna_index/silva-bac-16s-db:\
${work_path}/scripts/rRNA_databases/silva-bac-23s-id98.fasta,${work_path}/rrna_index/silva-bac-23s-db:\
${work_path}/scripts/rRNA_databases/silva-arc-16s-id95.fasta,${work_path}/rrna_index/silva-arc-16s-db:\
${work_path}/scripts/rRNA_databases/silva-arc-23s-id98.fasta,${work_path}/rrna_index/silva-arc-23s-db:\
${work_path}/scripts/rRNA_databases/silva-euk-18s-id95.fasta,${work_path}/rrna_index/silva-euk-18s-db:\
${work_path}/scripts/rRNA_databases/silva-euk-28s-id98.fasta,${work_path}/rrna_index/silva-euk-28s:\
${work_path}/scripts/rRNA_databases/rfam-5s-database-id98.fasta,${work_path}/rrna_index/rfam-5s-db:\
${work_path}/scripts/rRNA_databases/rfam-5.8s-database-id98.fasta,${work_path}/rrna_index/rfam-5.8s-db --reads ${pre_name}_merge.fastq \
--num_alignments 1 --fastx --aligned ${pre_name}_rRNA --other ${pre_name}_non_rRNA --paired_in --log -v -a $threads
# paired-in: If one of the paired-end reads is Aligned, put both reads into Aligned FASTA/Q file
${work_path}/scripts/unmerge-paired-reads.sh ${pre_name}_non_rRNA.fastq \
${pre_name}_non_rRNA_1.fastq ${pre_name}_non_rRNA_2.fastq
rm *merge.fastq
rm *non_rRNA.fastq
else
ln -s ../0.rawfq/${pre_name}_1.fastq ./${pre_name}_non_rRNA_1.fastq
ln -s ../0.rawfq/${pre_name}_2.fastq ./${pre_name}_non_rRNA_2.fastq
fi


# Step 2
# Trim reads using Trimmomatic
echo ${pre_name} Start to TRIM by trimmomatic..
cd ${work_path}/2.trim
java -jar ${work_path}/scripts/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads $threads \
../1.rm_rrna/${pre_name}_non_rRNA_1.fastq \
../1.rm_rrna/${pre_name}_non_rRNA_2.fastq \
-baseout ${pre_name}_nonrRNA_trimmed.fastq \
ILLUMINACLIP:${work_path}/scripts/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$read_length \
2> ${pre_name}_Trimmomatic.log
rm ../1.rm_rrna/${pre_name}_non_rRNA_1.fastq
rm ../1.rm_rrna/${pre_name}_non_rRNA_2.fastq

# Step 3
# Map reads to hg19 reference
echo ${pre_name} Start to MAP by hisat2..
cd ${work_path}/3.align
#$hisat2_path/
hisat2 -t -p $threads -q -x $hisat_index \
--rg-id ${pre_name} \
--rg SM:${pre_name} --rg LB:${pre_name} --rg PL:ILLUMINA \
-1 ../2.trim/${pre_name}_nonrRNA_trimmed_1P.fastq \
-2 ../2.trim/${pre_name}_nonrRNA_trimmed_2P.fastq \
-S ${pre_name}.sam \
2> ${pre_name}_alignment.log

# sam to bam
samtools view -S ${pre_name}.sam -b > ${pre_name}.bam
# 10min
time samtools sort  ${pre_name}.bam -o ${pre_name}_nsorted.bam
samtools index ${pre_name}_nsorted.bam

echo ${pre_name} start to COUNT by featureCounts..
cd ${work_path}/4.count
# Step 4
# Count number of reads on genes
featureCounts \
-T $threads \
-p \
-t exon \
-g gene_id \
-a $gtf_path \
-o ${pre_name}_featureCounts.txt \
../3.align/${pre_name}_nsorted.bam \
2> ${pre_name}_featureCounts.log

sed -i '1d' ${pre_name}_featureCounts.txt
awk '{print $1,$7}' ${pre_name}_featureCounts.txt >${pre_name}.count
        
if [ $flag == 1 ]; then
awk '{print $1}' ${pre_name}_featureCounts.txt >gene_expression_matrix.txt
join gene_expression_matrix.txt ${pre_name}.count --nocheck-order >gene_expression_matrix_${flag}.txt
else
flag_1=$((flag-1))
join gene_expression_matrix_${flag_1}.txt ${pre_name}.count --nocheck-order >gene_expression_matrix_${flag}.txt
fi

rm ../3.align/${pre_name}_nsorted.bam
rm ../3.align/${pre_name}_nsorted.bam.bai
rm ../3.align/${pre_name}.bam
rm ../3.align/${pre_name}.sam

flag=$((flag+1))
done

#最后整理表达矩阵
cd ${work_path}/4.count
flag=$((flag-1))
mv gene_expression_matrix_${flag}.txt expression_matrix.txt
rm gene_expression_matrix* *count
