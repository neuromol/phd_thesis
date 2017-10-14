### merge Plink file back together


import os
import sys
from subprocess import Popen, PIPE
import shlex
import re

import argparse
from itertools import groupby, cycle
from operator import itemgetter

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Put imputation file back together - Josh Atkins c3114203 uon.edu.au ' ,formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-c', '--chrom', help='Current chromsome',required='True')

    results = parser.parse_args(args)
    return (results.chrom)


def execute(cmd):
    """ to run commands on the shell and pipe out the display """
    plink_cmd_run = shlex.split(cmd)
    p = Popen(plink_cmd_run, shell=False, stdout=PIPE , stderr=PIPE)
    out, err = p.communicate()
    return (out)


def plink_1():
    """ basic QC - MISSING GT, HWE , MAF """
    print ("First run of plink ")
    plink_cmd =("plink --bfile chunk_2  --merge-list merge.list --make-bed --out temp")
    p = execute(plink_cmd)



def list_files():
    Popen(""" ls *.bed > bed.txt """, shell=True).wait()
    Popen("""ls *.bim > bim.txt  """, shell=True).wait()
    Popen("""ls *.fam > fam.txt """, shell=True).wait()
    Popen("""ls *.fam | cut -f 1 -d "." > chunks.txt """, shell=True).wait()
    Popen("""paste bed.txt bim.txt fam.txt | column -s $'\t' -t > merge.list """, shell=True).wait()

def plink_2():
    """ Remove multiallelic SNPs """
    Popen("""for i in $(cat chunks.txt) ; do plink --bfile $i --exclude temp-merge.missnp --make-bed --out temp/$i ; done """, shell=True).wait()


def main (chrom):
    os.chdir("chr"+ chrom)
    #Popen("""cd chr"""+ chrom, shell=True).wait()
    list_files()
    plink_1()
    Popen("""mkdir -p temp""" , shell=True).wait()
    plink_2()
    #Popen("""cd temp""" , shell=True).wait()
    #Popen("""cp ../merge.list temp/""", shell=True).wait()
    os.chdir("temp/")
    plink_cmd =("plink --bfile chunk_2  --merge-list ../merge.list --make-bed --out chr" + chrom)
    p = execute(plink_cmd)



if __name__ == '__main__':

    chrom = check_arg(sys.argv[1:])
    main(chrom )
