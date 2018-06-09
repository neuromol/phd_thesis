##### The creation of LoF Database


## Each vcf file had the same size header, extract the header information as so; 
head -n 80 100105.ann.vcf 


# use awk length to check if the reference or alt allele has a value greater than 1 = INDEL
tail -n +80 100105.ann.vcf | awk ' length($4) > 1 ||  length($5) > 1 '


# to extract SNPs 
tail -n +80 100105.ann.vcf | awk ' length($4) == 1 ||  length($5)== 1 '


###example of bash 
## all this was done via the bash, no python script was created to automate this (Time restrictions)

#PBS -l select=1:ncpus=1:mem=16G,walltime=01:00:00
cd /home/genomescratch1/old_stage8/new_ann
tail -n +80 100105.ann.vcf | awk ' length($4) == 1 &&  length($5)== 1 ' > 100105.ann.vcf.snps 
cat 100105.ann.vcf.header 100105.ann.vcf.snps > SNPs/100105.ann.vcf.snps


##INDEL filtration

#PBS -l select=1:ncpus=1:mem=16G,walltime=10:00:00
cd /home/genomescratch1/old_stage8/new_ann/INDELs/
cat 100105.ann.vcf.indel | java -jar /home/c3114203/apps/snpEff/SnpSift.jar filter " ( QD > 2.0 )  && ( FS < 200 )  && (ReadP
osRankSum > -20.0) && (SOR < 10.0 )" > filtered/100105.ann.vcf.indel.filtered 


## extraction of "high" impact indels (considered LoF in this thesis) 


#PBS -l select=1:ncpus=1:mem=16G,walltime=10:00:00
cd /home/genomescratch1/old_stage8/new_ann/INDELs/filtered/
cat 100105.ann.vcf.indel.filtered | java -jar /home/c3114203/apps/snpEff/SnpSift.jar filter "ANN[*].IMPACT = 'HIGH' " > LoF/1
00105.ann.vcf.indel.filtered.LoF 



### Annotation to dnNSFP

#PBS -l select=1:ncpus=1:mem=16G,walltime=10:00:00
cd /home/genomescratch1/old_stage8/new_ann/INDELs/filtered/LoF/
java -jar /home/c3114203/apps/snpEff/SnpSift.jar dbnsfp -db /home/c3114203/apps/snpEff/dbNSFP2.9.txt.gz -f SIFT_pred,Polyphen
2_HDIV_pred,Polyphen2_HVAR_pred,MetaLR_pred,ExAC_AC,ExAC_AF,clinvar_clnsig,clinvar_trait -v 100105.ann.vcf.indel.filtered.LoF
 > 100105.ann.vcf.indel.filtered.LoF.dbNSFP


## Put it into a clean format to convert to a database 
#PBS -l select=1:ncpus=1:mem=16G,walltime=10:00:00
cd /home/genomescratch1/old_stage8/new_ann/INDELs/filtered/LoF/
java -jar /home/c3114203/apps/snpEff/SnpSift.jar extractFields -e "." 100105.ann.vcf.indel.filtered.LoF.dbNSFP CHROM POS ID R
EF ALT EFF[0].GENE EFF[0].IMPACT EFF[0].EFFECT dbNSFP_ExAC_AC > clean/100105.ann.vcf.indel.filtered.LoF.dbNSFP.clean

#### All the clean data was merged into one file, using AWK to add the ID as the first field in the merged document  
ASRB_ID	#CHROM	POS	ID	REF	ALT	gene	IMPACT	EFFECT	ExAC_AC	clinvar   




#### LoF SNPS  
##
## Pipeline is similar to INDELS 

#PBS -l select=1:ncpus=1:mem=16G,walltime=01:00:00
cd /home/genomescratch1/old_stage8/new_ann/SNPs/
cat 100105.ann.vcf.snps | java -jar /home/c3114203/apps/snpEff/SnpSift.jar filter " ( QD > 2.0 )  && ( FS < 60 ) && (MQ > 40) && (MQRankSum > -12.5) && (ReadPos
RankSum > -8.0) && (SOR < 3.0 )" > filtered/100105.ann.vcf.snps.filtered


## this command was applied (its the same as the INDEL  ) 
cat 100105.ann.vcf.snps.filtered | java -jar /home/c3114203/apps/snpEff/SnpSift.jar filter "ANN[*].IMPACT = 'HIGH' " > LoF/1
00105.ann.vcf.snps.filtered.LoF 





### LoF annotation 
cd /home/genomescratch1/old_stage8/new_ann/SNPs/filtered/LoF/ 
java -jar /home/c3114203/apps/snpEff/SnpSift.jar dbnsfp -db /home/c3114203/apps/snpEff/dbNSFP2.9.txt.gz -f SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,Met
aLR_pred,ExAC_AC,ExAC_AF,clinvar_clnsig,clinvar_trait -v 100105.ann.vcf.snps.filtered.LoF > 100105.ann.vcf.snps.filtered.LoF.dbNSFP

**LoF vcf to text format** 
```java -jar /home/c3114203/apps/snpEff/SnpSift.jar extractFields -e "." 100105.ann.vcf.snps.filtered.LoF.dbNSFP CHROM POS ID REF ALT EFF[0].GENE EFF[0].IMPACT EFF
[0].EFFECT dbNSFP_ExAC_AC  dbNSFP_clinvar_trait dbNSFP_clinvar_clnsig dbNSFP_MetaLR_pred dbNSFP_Polyphen2_HDIV_pred > clean/100105.ann.vcf.snps.filtered.LoF.dbN
SFP.clean```

**once again this was merged together and went on to the Database_scripts**


##Missense (not part of thesis - Just added for future reference) 

**Missense filtering**
```#PBS -l select=1:ncpus=1:mem=16G,walltime=10:00:00
cd /home/genomescratch1/old_stage8/new_ann/SNPs/filtered/
 cat 100105.ann.vcf.snps.filtered | java -jar /home/c3114203/apps/snpEff/SnpSift.jar filter "ANN[*].EFFECT = 'missense_variant' " > Miss/100105.ann.vcf.snps.fil
tered.clean' ```


```#PBS -l select=1:ncpus=1:mem=16G,walltime=10:00:00
cd /home/genomescratch1/old_stage8/new_ann/SNPs/filtered/Miss/
java -jar /home/c3114203/apps/snpEff/SnpSift.jar dbnsfp -db /home/c3114203/apps/snpEff/dbNSFP2.9.txt.gz -f SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,Met
aLR_pred,ExAC_AC,ExAC_AF,clinvar_clnsig,clinvar_trait -v 100105.ann.vcf.snps.filtered.clean > 100105.ann.vcf.snps.filtered.clean.dbNSFP```


**Then this was merged, went on to Database_scripts**



