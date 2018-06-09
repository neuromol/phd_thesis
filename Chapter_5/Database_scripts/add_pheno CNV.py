#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import warnings
import sys
warnings.filterwarnings('ignore')
import argparse


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='add Cog_group Counts - Josh Atkins c3114203@uon.edu.au ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='file to add pli to',required='True')
    parser.add_argument('-p', '--pheno', help='phenotype_file', required='True')
    parser.add_argument('-o', '--output', help='output', required='True')
    results = parser.parse_args(args)
    return (results.input , results.pheno , results.output)


def get_phenotype(row , pheno_file):
    IDZ = row['ASRB_ID']
    if ((pheno_file['ID'] == IDZ).any()) :
        search = pheno_file[(pheno_file['ID'] == IDZ)]
        group = search.Group.values[0]
        print (search)
    else:
        group = 0
    return group


def main(file , phenotype_file , output):
    data1 = pd.read_csv(file , delimiter="\t")
    data1.columns=['ASRB_ID' ,'#CHROM',	'POS'	,'ID',	'REF',	'ALT','gene','IMPACT','EFFECT', 'ExAC_AC' ,'clinvar_trait','clinvar_clnsig' , 'MetaLR_pred','Polyphen2_HDIV_pred' , 'pLI']
    #data1.columns=["ASRB_ID" ,'#CHROM',	'POS'	,'ID',	'REF',	'ALT', 'gene' ,'IMPACT','EFFECT', 'ExAC_AC' , 'pLI' ]



    print (data1)
    Pheno_file = pd.read_csv(phenotype_file , delimiter="\t")
   # Pheno_file[['ID']] = Pheno_file[['ID']].astype(int)
    #data1[['ASRB_ID']] = data1[['ASRB_ID']].astype(int)
    data1['Group'] = data1.apply(lambda row: get_phenotype(row, Pheno_file),axis=1)
    #data1['case/con'] = data1['Group'].apply(convert_group)


    data1.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':

    file , phenotype_file , output = check_arg(sys.argv[1:])
    main(file , phenotype_file, output)
