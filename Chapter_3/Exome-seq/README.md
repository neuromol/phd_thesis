# Exome-seq/Whole genome pipeline

**This is a python script where the ids for each samples it put in ids.txt and in the folders.py folders, email address, SCP server (in case of this thesis it was neuromols local domain name) and batch number(very important to capture batch effects)** 

## To get the script running 

**The purpose of this script is to show my working in my thesis (not for this script to be used - its outdated)** 

**first create all folders manually that are defined in the folders.py eg stage1 stage2**

**in master_WGS.py change all GATK/JAVA/Reference files (This was lazy coding on my behalf)**

**remove my GATK key from some of the GATK packages - this caused many problems when our HPC was taken off internet access as parts of the GATK pipeline send out stats to their servers, we found that this would cause the walker to hang if it didn't recieve and echo back**

**Update the folders.py**

**run the following**
```
python master_WGS.py
```

to submit to HPC

```
for i in $(cat ids.txt) 
do 
qsub $i.sh
done
```

**This script was developed in 2014, while its been chopped and changed over the years - it is long overdue to be redeveloped**


