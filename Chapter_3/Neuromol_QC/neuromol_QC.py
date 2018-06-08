# -*- coding: utf-8 -*-
"""
Neuromol_QC.py is used to do basic QC on raw plink files
"""
import os
import sys
from subprocess import Popen, PIPE
import shlex
import re
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import cluster
import numpy as np
pd.options.mode.chained_assignment = None
import argparse
from itertools import groupby, cycle
from operator import itemgetter

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Neuromol_QC is a basic QC script that perfroms QC variants and sample QC and plots the results in Results directory - Josh Atkins c3114203 at uon.edu.au ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Input file name eg plink not plink.bim',required='True', default='plink')
    parser.add_argument('-l', '--label', help='The name of the dataset you are working on',required='True', default='default')
    parser.add_argument('-r', '--ref', help='The folder of the 1000g reference files and high LD regions file DEFAULT: ../Reference/ ', default='Reference/')
    parser.add_argument('-o', '--output', help='output file DEFAULT: Results_label/', default='Results')
    parser.add_argument('-1000g', '--diff_ref', help='if using a different name structure than 1000g_ref_panel state it here', default='1000g_ref_panel')
    parser.add_argument('-s', '--sex', help='Check sex from fam file DEFAULT: yes', default='yes')
    parser.add_argument('-b', '--inbreed', help='remove samples with > 0.2 inbreeding coefficent DEFAULT: yes', default='yes')
    parser.add_argument('-p', '--PCA', help='Do a PCA based off 1000g and remove outliers DEFAULT: yes', default='yes')
    parser.add_argument('-g', '--gwas', help='Do a GWAS of raw and clean dataset DEFAULT: yes', default='yes')
    results = parser.parse_args(args)
    return (results.input , results.label , results.ref , results.output , results.diff_ref , results.sex , results.inbreed , results.PCA , results.gwas)

def get_varz_removed(array):
    """ Code to work out how many variants have been removed at each plink stage """
    plink_VP = int(re.findall(r'\d+' ,str([x for x in array if re.search("variants and" , x )]))[0])
    plink_VL = int(re.findall(r'\d+' ,str([x for x in array if re.search("variants loaded" , x )]))[0])
    variants_removed = (plink_VL - plink_VP)
    return variants_removed

def get_ppl_removed(array):
    """ Works out how many people have been removed after each plink stage """
    plink_PL = int(re.findall(r'\d+' ,str([x for x in array if re.search("loaded from .fam" , x )]))[0])
    plink_PP = int(re.findall(r'\d+' ,str([x for x in array if re.search("people pass" , x )]))[1])
    people_removed = plink_PL - plink_PP
    return people_removed

def execute(cmd):
    """ to run commands on the shell and pipe out the display """
    plink_cmd_run = shlex.split(cmd)
    p = Popen(plink_cmd_run, shell=False, stdout=PIPE , stderr=PIPE)
    out, err = p.communicate()
    return (out)

def remove_ambig(raw):
    """ Remove Ambigous SNPs """
    print ("Stage One --- Removing Ambigous SNPs")
    Popen("""awk ' { if  ($5=="-" || $6=="-" ) print $2 }  ' """ + raw  + """.bim  > remove_missing.txt""", shell=True).wait()
    Popen("""awk ' { if (( $5=="T" && $6=="A") || ($5=="A" && $6=="T") || ( $5=="C" && $6=="G") || ( $5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ; } ' """ + raw + """.bim | grep ambig | awk ' {print $1 } ' >> remove_missing.txt""", shell=True).wait()
    plink_cmd = ("plink --bfile " +  raw + " --exclude remove_missing.txt --make-bed --out ambig")
    p = execute(plink_cmd)
    array_1 = p.split("\n")
    ppl_removed = get_ppl_removed(array_1)
    variants_removed = get_varz_removed(array_1)
    return ppl_removed , variants_removed

def basic_qc():
    """ basic QC - MISSING GT, HWE , MAF """
    print ("Stage Two --- Removing SNPs/Samples Geno 0.02 , HWE 0.000001 , MAF 0.01")
    plink_cmd =("plink --bfile ambig --geno 0.02 --maf 0.01 --hwe 0.0000001 --make-bed --out raw_qc")
    p = execute(plink_cmd)
    array_1 = p.split("\n")
    ppl_removed = get_ppl_removed(array_1)
    variants_removed = get_varz_removed(array_1)
    return ppl_removed , variants_removed

def sex_check(plink_file_stage):
    """ Sex Check - Impute Sex """
    print ("Stage Three --- Check the sex from .fam file")
    plink_cmd = ("plink --bfile " + plink_file_stage + " --check-sex --out sextest")
    p = execute(plink_cmd)
    Popen("""grep PROBLEM sextest.sexcheck | awk ' { print $1 , $2 } ' > sex.remove""", shell=True).wait()
    plink_cmd = ("plink --bfile " + plink_file_stage + " --remove sex.remove --make-bed --out sex_check ")
    p = execute(plink_cmd)
    array_1 = p.split("\n")
    ppl_removed = get_ppl_removed(array_1)
    variants_removed = get_varz_removed(array_1)
    return ppl_removed , variants_removed

def inbreeding_co(plink_file_stage):
    """Inbreeding Coefficient"""
    print ("Stage Four --- Removing individual with inbreeding coefficient > 0.2")
    plink_cmd = ("plink --bfile "+ plink_file_stage +" --het --out inbreeding")
    p = execute(plink_cmd)
    Popen(""" awk ' $6 > 0.2 ' inbreeding.het | awk ' { print $1 , $2 } ' | tail -n +2 > inbreeding.remove """, shell=True).wait()
    plink_cmd = ("plink --bfile " + plink_file_stage + " --remove inbreeding.remove --make-bed --out inbreeding_co ")
    p = execute(plink_cmd)
    array_1 = p.split("\n")
    ppl_removed = get_ppl_removed(array_1)
    variants_removed = get_varz_removed(array_1)
    return ppl_removed , variants_removed

def pop_strat(plink_file_stage , ref , diff_ref): #### need to finish this off
    """Population Stratification using PCA on 1000g reference"""
    ### remove high-levels of LD
    print ("Stage Five --- Population Stratification using PCA on 1000g data and K-means clustering to remove outliers")
    plink_cmd = ("plink --bfile " + plink_file_stage + " --exclude " + ref + "high-LD-regions-hg19.txt --range --indep-pairwise 50 5 0.2 --out prune")
    p = execute(plink_cmd)
    plink_cmd = ("plink --bfile " + plink_file_stage + " --extract prune.prune.in --genome --make-bed --out IBS")
    p = execute(plink_cmd)
    Popen("""awk ' $10 > 0.25 ' IBS.genome | awk ' { print $1 , $2 } ' | tail -n +2 > IBS.remove""",shell=True).wait()
    plink_cmd = ("""plink --bfile """+ plink_file_stage + """ --remove IBS.remove --make-bed --out temp5""")
    p = execute(plink_cmd)
    array_1 = p.split("\n")
    ppl_removed = get_ppl_removed(array_1)
    variants_removed = get_varz_removed(array_1)
    plink_cmd = ("plink --bfile IBS --write-snplist --out dataset")
    p = execute(plink_cmd)
    plink_cmd = ("plink --bfile " + ref + diff_ref + " --write-snplist --out 1000g")
    p = execute(plink_cmd)
    plink_cmd = ("plink --bfile IBS --extract 1000g.snplist --make-bed --out dataset1")
    p = execute(plink_cmd)
    plink_cmd = ("plink --bfile "+ ref + diff_ref + " --extract dataset.snplist --make-bed --out 1000g_temp")
    p = execute(plink_cmd)
    plink_cmd = ("plink --bfile 1000g_temp --bmerge dataset1 --make-bed --out combined")
    p = execute(plink_cmd)
    if os.path.isfile("combined.bed") is True:
        plink_cmd = ("plink --bfile combined --pca --out pca")
        p = execute(plink_cmd)
    else:
        plink_cmd = ("plink --bfile 1000g_temp --exclude combined-merge.missnp --make-bed --out 1000g_temp2")
        p = execute(plink_cmd)
        plink_cmd = ("plink --bfile dataset1 --exclude combined-merge.missnp --make-bed --out dataset2")
        p = execute(plink_cmd)
        plink_cmd = ("plink --bfile 1000g_temp2 --bmerge dataset2 --make-bed --out combined")
        p = execute(plink_cmd)
        plink_cmd = ("plink --bfile combined --pca --out pca")
        p = execute(plink_cmd)

    PCA_graph("pca.eigenvec" , label)

    Popen("""awk ' { print $1, $1  }  ' """ + label + """.outliers.txt > PCA.remove""",shell=True).wait()
    plink_cmd = ("""plink --bfile """ + plink_file_stage + """ --remove PCA.remove --make-bed --out """ + label + """_clean""")
    p = execute(plink_cmd)
    array_1 = p.split("\n")
    ppl_removed2 = get_ppl_removed(array_1)
    variants_removed2 = get_varz_removed(array_1)
    return ppl_removed , variants_removed , ppl_removed2 , variants_removed2

def association(plink_file_stage) :
    
    print("Stage Six --- Basic Association to check genomic inflation")
    plink_cmd = ("plink --bfile " + plink_file_stage + " --assoc --adjust --ci 0.95 --out clean")
    p = execute(plink_cmd)
    array_1 = p.split("\n")
    lambda_est = float(re.findall(r'\d+\.\d+' ,str([x for x in array_1 if re.search("Genomic inflation est" , x )]))[0])
    title = ( label + "_clean")
    QQ_plot("clean.assoc" , lambda_est , title)
    plink_cmd = ("plink --bfile "+ raw +" --assoc --adjust --ci 0.95 --out dirty")
    p = execute(plink_cmd)
    array_1 = p.split("\n")
    lambda_est = float(re.findall(r'\d+\.\d+' ,str([x for x in array_1 if re.search("Genomic inflation est" , x )]))[0])
    title = ( label + "_raw")
    QQ_plot("dirty.assoc", lambda_est , title )

def PCA_graph(INPUT_FILE, DATASET_LABEL):
    """uses the PCA plink does and generates plot and uses k-mean cluster to determine outliers to remove"""
    def SuperPop(x):
        if x in ["GBR" , "CEU" , "TSI" , "FIN" , "IBS" ]:
            return "EUR"
        elif x in ["CHB" , "JPT" , "CHS" , "CDX" , "KHV"]:
            return "EAS"
        elif x in ["YRI" , "LWK" , "GWD" , "MSL" , "ESN" , "ASW" , "ACB" ]:
                return "AFR"
        elif x in ["MXL" , "PUR" , "CLM" , "PEL" ]:
                return "AMR"
        elif x in ["GIH" , "PJL" , "BEB" , "STU" , "ITU" ]:
            return "SAS"
        else:
            return "Samples"
    ## Starting to handle big data so bringing in Pandas
    raw = pd.read_csv(INPUT_FILE, sep=" " , header=None )
    ## put 1000g data into superpopulation groups and define dataset
    clean = (raw[list(raw.columns[:4])])
    clean.columns = ['FAM_ID' , 'ID' , 'C1' , 'C2' ]
    clean.set_index(['FAM_ID'], inplace = True)
    ## setting up super population codes to map colours for graph
    clean["POP"] = clean.ID.apply(SuperPop)
    groups = clean.groupby('POP')
    ## Plotting
    fig , ax = plt.subplots()
    ax.margins(0.1)
    for name, group in groups:
        ax.plot(group.C1, group.C2, marker='o', linestyle='', ms=4, label=name)
    ax.legend(numpoints=1, loc='best')
    plt.xlabel('Component 1' )
    plt.ylabel('Component 2' )
    plt.suptitle("PCA on " + DATASET_LABEL , weight= 'bold')
    fig.savefig( DATASET_LABEL +".PCA_results.pdf")
    plt.close()
    ##kmean clustering to find outliers
    find_out = clean[['C1' , 'C2']].copy()
    k_means = cluster.KMeans(n_clusters=5,)
    k_means.fit(find_out)
    centroids = k_means.cluster_centers_
    labels = k_means.labels_
    results = pd.DataFrame([clean.index,labels]).T
    results.columns = ["FAM_ID" , "k_group"]
    results["ID"] = clean[["ID"]].copy()
    results.set_index(['FAM_ID'], inplace = True)
    output_label = (DATASET_LABEL + ".PCA_kmeans.txt")
    ## Display samples that are not Europeans in dataset
    merge_df = pd.merge(clean , results, right_index=True , left_index=True )
    merge_df['k_group'] = merge_df['k_group'].astype(int)
    test = merge_df.loc[merge_df['POP'] == "EUR" , ['k_group']].apply(np.median)
    Euro_group = int(test)
    #print ("European cluster is :" + str(Euro_group))
    your_samples = merge_df.loc[merge_df['POP'] == "Samples" , ['k_group']]
    your_samples['check'] = np.where(your_samples['k_group'] == Euro_group , 'good' , 'bad')
    bad_ids = your_samples[your_samples['check'] =='bad']
    after = (clean[~clean.index.isin(bad_ids.index)])
    count = len(bad_ids.index.get_level_values(0))
    #print (str(count) + " Samples fall outside the European cluster ")
    after_groups = after.groupby('POP')
    ### Plotting with outliers removed
    fig , ax = plt.subplots()
    ax.margins(0.1)
    for name, after_groups in after_groups:
        ax.plot(after_groups.C1, after_groups.C2, marker='o', linestyle='', ms=4, label=name)
    ax.legend(numpoints=1, loc='best')
    plt.xlabel('Component 1' )
    plt.ylabel('Component 2' )
    plt.suptitle("Outliers removed PCA on " + DATASET_LABEL + " - " + str(count) + " Samples were removed" ,weight= 'bold' )
    #print ("Graph saved as " + DATASET_LABEL + ".PCA_results.pdf")
    #print ("Outliers removed Graph saved as " + DATASET_LABEL + ".outliers_removed_PCA_results.pdf")
    fig.savefig( DATASET_LABEL +".outliers_removed_PCA_results.pdf")
    output_id = (DATASET_LABEL + ".outliers.txt")
    #print ("bad IDs exported to text file : " + output_id)
    bad_ids.to_csv(output_id , sep="\t" , header=None  )
    plt.close()

def QQ_plot(ass_results, lambda_est , title ):
    lambda_est = str(round(lambda_est,4))
    #Input Data
    data1 = pd.read_csv(ass_results , delim_whitespace=True , header=0)
    data1 = data1.dropna()
    data = data1.P.copy().values.flatten()
    data.sort()
    #work out n
    n = len(data)
    ##defining function
    def log_trans(x):
        t = -np.log10(x)
        return t
    def log_exp(x):
        t = np.log(x)
        return t
    ##to calculate the expected
    increment = 1.0/np.float64(n)
    current_expected = np.float64(1.0)
    ##empty array
    observed= []
    expected= []
    chi = []
    ##feed array with calculations
    for i in data:
        moo = (log_trans(i))
        observed.append(moo)
        woof = (log_trans(current_expected))
        expected.append(woof)
        current_expected -= increment
    #sort from lowest to highest
    observed.sort()
    expected.sort()
    #Work out axis for graph
    pmin = (min(data))
    pmax = (max(data))
    axisMax = ((log_trans(pmin)+ 1 ))
    # plot the arrays
    fig = plt.figure()
    plt.xlim([0,axisMax])
    plt.xlabel("Expected -log10(p)")
    plt.ylim([0,axisMax])
    plt.ylabel("Observed -log10(p)" )
    plt.title( title + " QQ Plot - est Lambda = " +lambda_est )
    ##    # the observed vs. expected data
    dataAx = fig.add_subplot(1,1,1)
    dataAx.plot(expected,observed,'r.') # red dots
    # a diagonal line for comparison
    lineAx = fig.add_subplot(1,1,1)
    lineAx.plot([0,axisMax],[0,axisMax],'b-') # blue line
    fig.savefig(title + ".QQ-plot.png", dpi=200)
    plt.close()

def main(raw, label , ref , output , diff_ref , sex , inbreeding , PCA , gwas):

    if "no" in (sex , inbreeding, PCA, gwas):
    #lets work out which test to skip
        s1_p, s1_v = (remove_ambig(raw))
        s2_p, s2_v = (basic_qc())
       
        dataset = "raw_qc"
        if sex == "yes":
            s3_p, s3_v = (sex_check(dataset))
            dataset = "sex_check"
        elif sex != "yes":
            s3_v = 0
            s3_p = 0
            
            dataset ="raw_qc"
        if inbreeding =="yes":
          
            s4_p, s4_v = inbreeding_co(dataset)
            dataset == "inbreeding_co"
        elif inbreeding !="yes":
            dataset == "raw_qc"
            s4_v = 0
            s4_p = 0
        if PCA =="yes":
            s5_p1, s5_v1 , s5_p2 , s5_v2  = (pop_strat(dataset , ref , diff_ref))
            dataset == str(label + "_clean")
        elif PCA!="yes":
            s5_v = 0
            s5_p = 0
        if gwas =="yes":
            association(dataset)
    else:

        s1_p, s1_v = (remove_ambig(raw))
        s2_p, s2_v = (basic_qc())
        s3_p, s3_v = (sex_check("raw_qc"))
        s4_p, s4_v = (inbreeding_co("sex_check"))
        s5_p1, s5_v1 , s5_p2 , s5_v2 = (pop_strat("inbreeding_co",ref , diff_ref ))
        dataset = (label + "_clean")
        association(dataset)
    print ("Cleaning Up...........................")
    
    total_var = (s1_v + s2_v + s3_v + s4_v + s5_v1 + s5_v2 )
    total_ppl = ( s1_p + s2_p + s3_p + s4_p + s5_p1 + s5_p2 )
    
    
    
    
    #clean up time
    Popen("""cat *.remove >  IDs_removed.txt  """,shell=True).wait()
    Popen("""rm sex* combined* inbreeding* dataset* IBS* ambig* pca* temp5* prune* 1000g* raw* PCA* clean.hh clean.log dirty.hh dirty.log remove_missing.txt  """,shell=True).wait()


    ## print stats
    print (" ")
    line = ("............................................................................................................")
    ambigstat = (" Number of variants from Ambigious SNPs stage removed: " + str(s1_v)  + " People removed: " + str(s1_p) )
    QCstats = (" Number of variants from basic QC stage SNPs removed: " + str(s2_v)  + " People removed: " + str(s2_p) )
    sexstats = (" Number of Variants from the sex check SNPs removed: " + str(s3_v)  + " People removed: " + str(s3_p) )
    inbreedstats = (" Number of variants from inbreeding coefficient test removed: " + str(s4_v)  + " People removed: " + str(s4_p) )
    IBSstats = (" Number of variants from Pairwise IBS test removed : " + str(s5_v1)  + " People removed: " + str(s5_p1) )
    popstrat = (" Number of variants from Population Stratification test removed : " + str(s5_v2)  + " People removed: " + str(s5_p2) )
   
    totstats = (" Total variants removed: " + str(total_var) + " Total People removed: " + str(total_ppl) )


    print (line)
    print (ambigstat)
    print (QCstats)
    print (sexstats)
    print (inbreedstats)
    print (IBSstats)
    print (popstrat)
    print (line)
    print (totstats)

    f = open ('Output_stats.txt', 'w')
    f.write(line + "\n" )
    f.write(ambigstat  + "\n")
    f.write(QCstats + "\n")
    f.write(sexstats + "\n")
    f.write(inbreedstats + "\n")
    f.write(IBSstats + "\n")
    f.write(popstrat + "\n")
    f.write(line + "\n")
    f.write(totstats + "\n")

    #output directory
    Popen("""mkdir -p """ + output + "_" + label ,shell=True).wait()
    Popen("""mv """ +label + ".* " + output +"_" + label  ,shell=True).wait()
    Popen("""mv """ +label + "_* " + output +"_" + label  ,shell=True).wait()
    Popen("""mv clean* """ + output +"_" + label  ,shell=True).wait()
    Popen("""mv dirty* """ + output +"_" + label ,shell=True).wait()
    Popen("""mv IDs_removed.txt """ + output +"_" + label  ,shell=True).wait()
    Popen("""mv Output_stats.txt """ + output +"_" + label  ,shell=True).wait()



if __name__ == '__main__':
    raw, label , ref , output , diff_ref , sex , inbreeding , PCA , gwas = check_arg(sys.argv[1:])

    main(raw, label , ref , output , diff_ref , sex , inbreeding , PCA , gwas )



