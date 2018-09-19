#!/usr/bin/env python


import pandas as pd
import numpy as np
import warnings
import sys
warnings.filterwarnings('ignore')
import argparse


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='To plot and do stats on PRS - Josh Atkins c3114203@uon.edu.au ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='file to add pli to',required='True')
    parser.add_argument('-p', '--pli', help='pli_file', required='True')
    parser.add_argument('-o', '--output', help='output', required='True')
    results = parser.parse_args(args)
    return (results.input , results.pli , results.output)


def get_pli(row , pheno_file):
    IDZ = row['gene']
    if ((pheno_file['gene'] == IDZ).any()) :
        search = pheno_file[(pheno_file['gene'] == IDZ)]
        group = search.pLI.values[0]
    else:
        group = "."
    return group



def main(file , pli , output ):
    data1 = pd.read_csv(file , sep = "\s+|\t+|\s+\t+|\t+\s+",  header = 0 )
    data1.columns=['ASRB_ID ' ,'#CHROM',	'POS'	,'ID',	'REF',	'ALT',	'gene',	'IMPACT',	'EFFECT',	'ExAC_AC' ,"clinvar_trait","clinvar_clnsig" , "dbNSFP_SIFT_pred","Polyphen2_HDIV_pred" ]
    #data1.columns=['ASRB_ID ' ,'#CHROM',	'POS'	,'ID',	'REF',	'ALT',	'gene',	'IMPACT',	'EFFECT',	'ExAC_AC'  ]

#    data1.columns = df.columns.str.replace('EFF[0].GENE', 'gene')
#    data1.columns = df.columns.str.replace('EFF[0].IMPACT', 'IMPACT')
#    data1.columns = df.columns.str.replace('EFF[0].EFFECT', 'EFFECT')
#


    data2 = pd.read_csv(pli , sep='\t', header = 0 )

    data1['pLI'] = data1.apply(lambda row: get_pli(row,data2),axis=1)

    data1.to_csv(output , sep='\t' , index=False)














if __name__ == '__main__':

    file , pli , output = check_arg(sys.argv[1:])
    main(file , pli, output)
