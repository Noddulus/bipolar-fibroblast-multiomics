# Data Preprocessing 

**Note:** The original batch submission scripts executed on the computing cluster are no longer available. However, the exact tools, versions, and parameters used to process the raw `.fastq` files for this project are documented below for full reproducibility.

### Removing low quality reads
* **Tools:** TrimGalore, Cutadapt
* **Execution:**
```bash
/apps/TrimGalore-0.6.4/trim_galore --paired -o /trimmed_seqs --path_to_cutadapt /bin/cutadapt /untrimmed_seqs/*gz
```
### Alignment
* **Tools:** samtools, bowtie2
* **Execution:**
```bash
#Defining the trimmed fastq files for this alignment
fastq_files=(*R1_001_val_1.fq.gz)
this_fastq=${fastq_files[(SGE_TASK_ID-1)]}
#Replace the ending of the validated trimmed file for the second read
fastq1=`echo ${this_fastq}`
fastq2=`echo ${this_fastq} | sed 's/R1_001_val_1.fq.gz/R2_001_val_2.fq.gz/g'`
prefix=`echo ${this_fastq} | sed 's/_R1_001_val_1.fq.gz//g'`

#Export the BOWTIE2_INDEXES variable leading to directory with index
BOWTIE2_INDEXES=/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index
export BOWTIE2_INDEXES

#allow up to 2kb insert size
#up to 4 alignments
#the best paired alignment is retained in future filtering
bowtie2 \
-x genome \
-1 $fastq1 \
-2 $fastq2 \
-k 4 \
-X 2000 \
--local \
-p 4 \
-S $prefix.sam \
2>$prefix.align.log

samtools view -S -h -b $prefix.sam | samtools sort -T $prefix.raw -o $prefix.raw.bam

samtools index ${prefix}.raw.bam

#get stats for mtDNA and other alignment info
samtools idxstats ${prefix}.raw.bam > ${prefix}.raw.bam.stats
```
### Removing PCR duplicates
* **Tools:** samtools, picard tools, Anaconda
* **Execution:**
```bash
bam_files=(./*raw.bam)

this_bam=${bam_files[($SGE_TASK_ID-1)]}

prefix=`echo ${this_bam} | sed 's/.raw.bam$//'`

#keep proper pair
#filter out mate unmapped
#filter out poor quality
samtools view -bh -o ${prefix}.filt.bam -f 2 -F 524 -q 30 ${this_bam}

#sort by read name, fixmate for marking duplicates
samtools sort -o ${prefix}.filt.qnameSort.bam -n -T ${prefix}.raw ${prefix}.filt.bam

samtools fixmate ${prefix}.filt.qnameSort.bam ${prefix}.fixmate.bam

samtools view -bh -o ${prefix}.fixmate.filt.bam -f 2 -F 1804 ${prefix}.fixmate.bam

samtools sort -o ${prefix}.toMark.bam -T ${prefix}.filt.raw ${prefix}.fixmate.filt.bam

rm ${prefix}.fixmate.filt.bam

rm ${prefix}.fixmate.bam

rm ${prefix}.filt.qnameSort.bam

rm ${prefix}.filt.bam

#Mark duplicates

#Previous working picard INPUt=${this_bam}
java -Xmx4G -jar /u/local/apps/picard-tools/2.9.0/picard.jar MarkDuplicates \
INPUT=${prefix}.toMark.bam \
OUTPUT=${prefix}.markedDup.bam \
METRICS_FILE=${prefix}.markedDup.metrics.txt \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=LENIENT

rm ${prefix}.toMark.bam
```
### Assigning read counts to Peaks (ATAC-seq data)
