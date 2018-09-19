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
    parser.add_argument('-p', '--pheno', help='phenotype_file')
    parser.add_argument('-o', '--output', help='output', required='True')
    results = parser.parse_args(args)
    return (results.input , results.pheno , results.output)



def get_var_count_Cont(row , count):
    IDZ = row['ID']
   
    search = count[(count['ID'] == IDZ)]
    cont = search.Cont.values[0]
    CD = search.CD.values[0]
    CS = search.CS.values[0]
    total = search.All.values[0]
#    print (cont + CD + CS)
  #  print (search)
   # print (search)
    return pd.Series([ cont , CD , CS , total ]) 


#def get_var_count_CD(row , count):
#    IDZ = row['ID']
#   
#    search = count[(count['ID'] == IDZ)]
#  
#    CD = search.CD.values[0]
##    CS = search.CS.values[0]
##    print (cont + CD + CS)
#    return CD
#
#
#
#def get_var_count_CS(row , count):
#    IDZ = row['ID']
#   
#    search = count[(count['ID'] == IDZ)]
#
##    CD = search.CD.values[0]
#    CS = search.CS.values[0]
##    print (cont + CD + CS)
#    return CS













def main(file , phenotype_file , output ):
    data1 = pd.read_csv(file , sep='\t' , header = 0 )
    #data1.columns=['ASRB_ID'	,'#CHROM',	'POS'	,'ID',	'REF',	'ALT',	'gene',	'IMPACT','EFFECT','ExAC_AC','clinvar', 'pLI' , 'Group' , 'case/con' ]
    
    print ( len(data1.index))
    
    data1 = data1[data1.IMPACT == "HIGH"]
    
    print ( len(data1.index))
    
    
    #variants = data1.filter(['ASRB_ID' , 'ID' , 'Group'] , axis=1 ) 

#    

    #data1['Count'] = data1.groupby(['Group'])['ID'].transform('count')

   # data1.to_csv(output, sep='\t', index=False)
    
    count = pd.crosstab(data1['ID'], data1['Group'] ,margins=True)
    
    count.to_csv("variant_count.txt" ,sep='\t' )
    
    count = pd.read_csv("variant_count.txt" , sep='\t' , header = 0 )
    count.columns = ['ID' , 'Cont' , 'CD' , 'CS' , 'All' ]
   
    new = pd.merge(data1 , count, left_on="ID" , right_on="ID" )
    #print (count.columns)
   # data1['Cont_var_count'] , data1['CD_var_count'] , data1['CS_var_count'] , data1['total_var_count'] = data1.apply(lambda row: get_var_count_Cont(row, count),axis=1)
    #data1['CD_var_count'] = data1.apply(lambda row: get_var_count_CD(row, count),axis=1)
    #data1['CS_var_count'] = data1.apply(lambda row: get_var_count_CS(row, count),axis=1)
    print (new)
    
   # print( data1.apply(lambda row: get_var_count_Cont(row, count),axis=1))

    new.to_csv(output, sep='\t', index=False)

 

   

if __name__ == '__main__':

    file , phenotype_file , output = check_arg(sys.argv[1:])
    main(file , phenotype_file, output)