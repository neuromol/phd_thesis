# RNA-seq Protocol

*method adapted from https://www.nature.com/articles/nprot.2012.016*

**After checking Fastqs with FastQC and removing over represented reads and adaptors with Cutadapt**


## Mapping and sorting

**The following commands were forwarded to an automated script to produce one mapped and sorted BAM file for each sample (which were in moved so they were in their own folder) :**

**note that the IDs have been changed to sampleID**

```
#!/bin/bash
 #PBS -l select=1:ncpus=4:mem=16G,walltime=5:00:00
cd  /home/genomescratch1/RNAseq
hisat2 -p 4 --dta-cufflinks -x /home/c3114203/ref/genome -1  sampleID_L1_1.fq -2 sampleID_L1_2.fq -S sampleID.sam 
samtools view -@ 4 -S -b sampleID.sam >  sampleID.bam 
samtools sort -@ 4 sampleID.bam sampleID.sort
samtools index sampleID.sort.bam
```

## CuffLinks

**Transcript assembling to generate gtf for each case with cufflinks, this step and the cuffmerge step can be skipped if you are only interested in known transcripts (just need to supply .gtf file for the reference that you aligned to)**

```
#!/bin/bash
 #PBS -l select=1:ncpus=4:mem=15G,walltime=24:00:00
cd /home/genomescratch1/RNAseq/sampleID/
cufflinks -p 4 -g ~/ref/genes.gtf sampleID.sort.bam
```

## CuffMerge

**The location for each sample's .gtf file was detailed into a text file**

```
#!/bin/bash
 #PBS -l select=1:ncpus=8:mem=30G,walltime=24:00:00
cd /home/genomescratch1/RNAseq/
cuffmerge -p 8 -g ~/ref/genes.gtf list.txt
```

## CuffDiff

**Differential Expression**

```
#!/bin/bash
 #PBS -l select=1:ncpus=8:mem=30G,walltime=500:00:00
cd /home/genomescratch1/RNAseq/
cuffdiff -p 8 -o results/ -L SZ,control /home/genomescratch1/RNAseq/merged_asm/merged.gtf *each_case.bam (separated by , ) *each_control
```

## Individual analysis

**Each case's gtf was merged with ever controls to create individual gtf files, with cuffdiff performed on each case compared to every control**











