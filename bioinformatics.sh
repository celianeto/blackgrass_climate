### trimming of fastq with trimmomaticPE 
trimmomaticPE - "LEADING:5" - "TRAILING:5" - "SLIDINGWINDOW:4:20" - "MINLEN:50"


### mapping trimmed fastq to genome using bwa 
### this needs to be done in loop for all samples to be proccessed
### replace SAMPLENAME with sample name and path to trimmed files (output from above)
bwa mem -R '@RG\tID:SAMPLENAME\tSM:SAMPLENAME\tPL:ILLUMINA' PATH_TO_REF_FASTA PATH_TO_SAMPLE/SAMPLENAME1.1.fastq.gz PATH_TO_SAMPLE/SAMPLENAME1.2.fastq.gz | samtools sort -m 180G -o PATH_TO_OUTPUT/SAMPLENAME1.sorted.sam


### sam to bam conversion using samtools 
samtools view -@ 40 -q 20 -S -b PATH_TO_OUTPUT/SAMPLENAME1.sorted.sam > PATH_TO_OUTPUT/SAMPLENAME1.bam
samtools sort -@ 40 PATH_TO_OUTPUT/SAMPLENAME1.bam -o PATH_TO_OUTPUT/SAMPLENAME1.sorted.bam


### mark and remove dupliactes
picard MarkDuplicates INPUT=PATH_TO_OUTPUT/SAMPLENAME1.sorted.bam OUTPUT=PATH_TO_OUTPUT/SAMPLENAME1.dedup.bam METRICS_FILE=PATH_TO_OUTPUT/SAMPLENAME1.metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
### re-index
samtools index PATH_TO_OUTPUT/SAMPLENAME1.dedup.bam


### converting bam to mpileup 
### create different files as needed -- samples to include should follow one after another
### these files will be used as input in npstat 
### they will also be converted to sync files for other downstream analyses (see below)
### do for all chrs
samtools mpileup -r Chr1 -f PATH_TO_REF_FASTA PATH_TO_OUTPUT/SAMPLENAME1.dedup.bam PATH_TO_OUTPUT/SAMPLENAME2.dedup.bam PATH_TO_OUTPUT/SAMPLENAME3.dedup.bam PATH_TO_OUTPUT/SAMPLENAMEn.dedup.bam > allPops.chr1.mpileup


### converting mpileup to sync 
### sync files are needed as input for poolfstat, see below
### do for all chrs
java -ea -Xmx180g -jar PATH_TO_/popoolation2-code/mpileup2sync.jar --input allPops.chr1.mpileup --output allPops.chr1.sync --fastq-type sanger --min-qual 20 --threads 40


#### converting sync to genobaypass format
#### in R
### see poolfstat manual for details
library(poolfstat)

# sync to pool object
pool = popsync2pooldata("allPops.chr1.sync",poolsizes=rep(50,64), min.rc = 5, min.maf = 0.05)

# pool object to genobaypass format 
pooldata2genobaypass(pool, writing.dir = getwd())



