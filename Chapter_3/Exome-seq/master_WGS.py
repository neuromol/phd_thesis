

#########################################################
#   WGS full script (and Exome)
#
#   Set varibles in folders.py
#   Set sample numbers in all.txt
#
##########################################################




import os, sys

from folders import *


# Stage 1 - Mapping with Bowtie2
f = open('ids.txt','r')
string = ""
while 1:

    line = f.readline()
    if not line:break
    string = line


    fileA = line[0:-1] + BATCHNUM  +  ".sh"



    D = open(fileA,"w")
    D.write("#!/bin/bash \n")
    D.write("#PBS -l select=1:ncpus=8:mem=32G,walltime=500:00:00 \n")
    D.write("cd ")
    D.write(MYPATH)
    D.write(" \n")
    D.write("#PBS -N ")
    D.write(line[0:-1])
    D.write("\n")
    D.write("export PATH=/usr/local/jdk1.8.0_45/:$PATH \n")


#### starting time

    D.write("echo \" ")
    D.write(line[0:-1])
    D.write(" $(date) \"  >> ")
    D.write(MYPATH)
    D.write("startime.txt \n")

########################################################
##### Stage 1 and 2

#bowtie mapping
    D.write("bowtie2 -p 8 --very-sensitive -x /home/c3114203/ref/genome  -1 ")
    D.write(RAWPATH)
    D.write(line[0:-1])
    D.write("_R1.fastq.gz -2 ")
    D.write(RAWPATH)
    D.write(line[0:-1])
    D.write("_R2.fastq.gz -S ")
    D.write(STAGE1)
    D.write(line[0:-1])
    D.write(".sam \n")

#writing to log confirm Mapping
    D.write("echo \" ")
    D.write(line[0:-1])
    D.write(" $(date) \"  >> ")
    D.write(MYPATH)


    D.write("stage1log.txt \n")

#Samtools SAM2BAM
    D.write("samtools view -@ 8 -S -b ")
    D.write(STAGE1)
    D.write(line[0:-1])
    D.write(".sam > ")
    D.write(PREBAM)
    D.write(line[0:-1])
    D.write(".bam \n")


#Samtools sort
    D.write("samtools sort -@ 8 ")
    D.write(PREBAM)
    D.write(line[0:-1])
    D.write(".bam  ")
    D.write(STAGE2)
    D.write(line[0:-1])
    D.write(".sort \n")


#writing to log confirm Mapping
    D.write("echo \" ")
    D.write(line[0:-1])
    D.write(" $(date) \"  >> ")
    D.write(MYPATH)
    D.write("stage2log.txt \n")



#remove samfile
   # D.write("rm ")
   # D.write(STAGE1)
   # D.write(line[0:-1])
   # D.write(".sam  \n")

#Mark Dups
    D.write("java -Xmx32g -jar  /home/c3114203/apps/picard-tools-1.99/MarkDuplicates.jar INPUT= ")
    D.write(STAGE2)
    D.write(line[0:-1])
    D.write(".sort.bam  OUTPUT= ")
    D.write(STAGE2)
    D.write(line[0:-1])
    D.write(".dups.bam METRICS_FILE=")
    D.write(STAGE2)
    D.write(line[0:-1])
    D.write(".metrics.txt \n")


############################################


############################################
### Stage 3

#Picard tools
    D.write("java -Xmx32g -jar  /home/c3114203/apps/picard-tools-1.99/AddOrReplaceReadGroups.jar I= ")
    D.write(STAGE2)
    D.write(line[0:-1])
    D.write(".dups.bam  O= ")
    D.write(STAGE3)
    D.write(line[0:-1])
    D.write(".header.bam RGPL=ILLUMINA RGCN=Garvin RGDT=2015-12 RGID=")
    D.write(line[0:-1])
    D.write(" RGLB=")
    D.write(line[0:-1])
    D.write(" RGPU=")
    D.write(line[0:-1])
    D.write(" RGSM=")
    D.write(line[0:-1])
    D.write(" \n")


#Samtools Index
    D.write("samtools index ")
    D.write(STAGE3)
    D.write(line[0:-1])
    D.write(".header.bam \n")

#writing to log confirm Mapping
    D.write("echo \" ")
    D.write(line[0:-1])
    D.write(" $(date) \"  >> ")
    D.write(MYPATH)
    D.write("stage3log.txt \n")

#remove sort.bam
    D.write("rm ")
    D.write(STAGE2)
    D.write(line[0:-1])
    D.write(".sort.bam \n")

#Remove pre_bam
    D.write("rm ")
    D.write(PREBAM)
    D.write(line[0:-1])
    D.write(".bam \n ")



#remove dup.bam
    D.write("rm ")
    D.write(STAGE2)
    D.write(line[0:-1])
    D.write(".dups.bam \n")

###########################################

###########################################
#### Stage 4




#indel target Creat0r - stage 4
    D.write("java -Xmx32g -jar /home/c3114203/apps/GATK/GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R /home/c3114203/ref/genome.fa   -I ")
    D.write(STAGE3)
    D.write(line[0:-1])
    D.write(".header.bam -known /home/c3114203/ref/indels.vcf -o ")
    D.write(STAGE4)
    D.write(line[0:-1])
    D.write(".intervals \n")


#indel target Creat0r
    D.write("java -Xmx32g -jar /home/c3114203/apps/GATK/GenomeAnalysisTK.jar  -T IndelRealigner -R /home/c3114203/ref/genome.fa  -I ")
    D.write(STAGE3)
    D.write(line[0:-1])
    D.write(".header.bam -known /home/c3114203/ref/indels.vcf -targetIntervals ")
    D.write(STAGE4)
    D.write(line[0:-1])
    D.write(".intervals -o ")
    D.write(STAGE4)
    D.write(line[0:-1])
    D.write(".indel.bam \n")

#writing to log confirm Mapping
    D.write("echo \" ")
    D.write(line[0:-1])
    D.write(" $(date) \"  >> ")
    D.write(MYPATH)
    D.write("stage4log.txt \n")

#remove header.bam
    D.write("rm ")
    D.write(STAGE3)
    D.write(line[0:-1])
    D.write(".header.bam \n")


# Stage 5 Part one
    D.write("java -Xmx32g -jar /home/c3114203/apps/GATK/GenomeAnalysisTK.jar -et NO_ET -K /home/c3114203/apps/GATK/joshua.atkins_uon.edu.au.key -nct 8 -T BaseRecalibrator -R /home/c3114203/ref/genome.fa   -I ")
    D.write(STAGE4)
    D.write(line[0:-1])
    D.write(".indel.bam  -knownSites /home/c3114203/ref/indels.vcf -knownSites /home/c3114203/ref/dbsnp_137.hg19.vcf -o ")
    D.write(STAGE5)
    D.write(line[0:-1])
    D.write(".grp \n")
        
#Stage 5 Part two

    D.write("java -Xmx32g -jar /home/c3114203/apps/GATK/GenomeAnalysisTK.jar -et NO_ET -K /home/c3114203/apps/GATK/joshua.atkins_uon.edu.au.key -nct 8 -T PrintReads -R /home/c3114203/ref/genome.fa -I ")
    D.write(STAGE4)
    D.write(line[0:-1])
    D.write(".indel.bam  -BQSR ")
    D.write(STAGE5)
    D.write(line[0:-1])
    D.write(".grp -o ")
    D.write(STAGE5)
    D.write(line[0:-1])
    D.write(".recal.bam \n")
                
                
#writing to log confirm Mapping
    D.write("echo \" ")
    D.write(line[0:-1])
    D.write(" $(date) \"  >> ")
    D.write(MYPATH)
    D.write("stage5log.txt \n")
                
#remove indel.bam
    D.write("rm ")
    D.write(STAGE4)
    D.write(line[0:-1])
    D.write(".indel.bam \n")
                
## Stage 6
#Varient Calling - stage 6
    D.write("java -Xmx32g -jar /home/c3114203/apps/GATK/GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R /home/c3114203/ref/genome.fa  -I ")
    D.write(STAGE5)
    D.write(line[0:-1])
    D.write(".recal.bam  -minPruning 5 -o ")
    D.write(STAGE6)
    D.write(line[0:-1])
    D.write(".vcf \n")
                
                
#writing to log confirm Mapping
    D.write("echo \" ")
    D.write(line[0:-1])
    D.write(" $(date) \"  >> ")
    D.write(MYPATH)
    D.write("stage6log.txt \n")
        
####### Stage 7 and 8
                
#Varient Annotation dbsnp stage 7
    D.write("java -Xmx32g -jar /home/c3114203/apps/GATK/GenomeAnalysisTK.jar  -T VariantAnnotator -R /home/c3114203/ref/genome.fa  -I ")
    D.write(STAGE5)
    D.write(line[0:-1])
    D.write(".recal.bam -V ")
    D.write(STAGE6)
    D.write(line[0:-1])
    D.write(".vcf -o ")
    D.write(STAGE7)
    D.write(line[0:-1])
    D.write(".dbsnp.vcf --dbsnp /home/c3114203/ref/dbsnp_137.hg19.vcf \n")
                
                
#writing to log confirm Mapping
    D.write("echo \" ")
    D.write(line[0:-1])
    D.write(" $(date) \"  >> ")
    D.write(MYPATH)
    D.write("stage7log.txt \n")
                
                
#Functional Annotation snpEff - stage 8
    D.write("java -Xmx32g -jar /home/c3114203/apps/snpEff/snpEff.jar  -csvStats ")
    D.write(STAGE7)
    D.write(line[0:-1])
    D.write(".csv hg19 ")
    D.write(STAGE7)
    D.write(line[0:-1])
    D.write(".dbsnp.vcf > ")
    D.write(STAGE8)
    D.write(line[0:-1])
    D.write(".ann.vcf \n")
                
                
                
##Pull out High-risk genes and generate high-gene list
    D.write("cat ")
    D.write(STAGE8)
    D.write(line[0:-1])
    D.write(".ann.vcf | java -jar /home/c3114203/apps/snpEff/SnpSift.jar extractFields - CHROM POS ID REF ALT \"ANN[*].GENE\" \"ANN[*].EFFECT\" \"ANN[*].IMPACT\" | grep \"HIGH\"  > " )
    D.write(STAGE8)
    D.write(line[0:-1])
    D.write(".high.txt \n")
                
                
#writing to log confirm Mapping
    D.write("echo \" ")
    D.write(line[0:-1])
    D.write(" $(date) \"  >> ")
    D.write(MYPATH)
    D.write("stage8log.txt \n")
                
            
                
####email
    D.write("mail -s \"Job ")
    D.write(line[0:-1])
    D.write(" completed \" ")
    D.write(EMAILAD)
    D.write(" < /dev/null \n")
                
                
#### Finish time
                
    D.write("echo \" ")
    D.write(line[0:-1])
    D.write(" $(date) \"  >> ")
    D.write(MYPATH)
                
                
    D.write("endtime.txt \n")
                
                
###########################################
                
                
                

