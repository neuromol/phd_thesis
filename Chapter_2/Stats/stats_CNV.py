#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 15:57:18 2017

@author: dickie_ho
"""
import matplotlib.pyplot as plt 
import statsmodels.api as sm 
import statsmodels.formula.api as smf
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import scipy.stats as stats
import warnings
import sys
warnings.filterwarnings('ignore')

plt.style.use('ggplot')
import argparse


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='To plot and do stats on Score - Josh Atkins c3114203@uon.edu.au ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Clean_CNV File',required='True')
    parser.add_argument('-l', '--label', help='label for dataset', default='Score')
    parser.add_argument('-o', '--output', help='location for output, default is current path', default='./')
    results = parser.parse_args(args)
    return (results.input , results.label , results.output)



def convert_group(x): 
    #convert labels to cont vars 
    if x == str("Control"):
        return 0
    elif x == str("Spared Cognition"):
        return 1 
    elif x == str("Impaired cognition"):
        return 2

def check_variance(x, y):
    #check if variance is equal between dists 
    z_val , p_val = stats.ranksums(x , y)

    return p_val

def t_test(x , y , p):
    ## work out which T-test to run 
    if p > 0.05: 
        Z , t_val = stats.ttest_ind(x , y, equal_var=True)
        return t_val
    elif p < 0.05:
        z , t_val = stats.ttest_ind(x , y, equal_var=False) 
        return t_val

def get_max_min(df):
    ##first part of making a binary variable is to find out the range
    min_num = df.Case.min()
    max_num = df.Case.max()
    return min_num , max_num 

def make_binary(x , min_num , max_num):
    #converts top cont varible to 1, lowest to 0  
    if x == min_num :
        return 0 
    elif x == max_num :
        return 1 

def main(file , name , output ):
    data1 = pd.read_csv(file , sep='\t', index_col = 0 , header = 0 )
    data1['Case'] = data1['Group'].apply(convert_group)
    scaler = StandardScaler()
    data1[['Score']] = scaler.fit_transform(data1[['Score']])
    control_Score = data1[data1['Case']==0]
    CD_Score = data1[data1['Case']==2]
    CS_Score = data1[data1['Case']==1]
    
    ### check variants between groups
    p_val_var_Cont_CD = check_variance(control_Score['Score'], CD_Score['Score'])
    p_val_var_Cont_CS = check_variance(control_Score['Score'], CS_Score['Score'])
    p_val_var_CS_CD = check_variance(CS_Score['Score'], CD_Score['Score'])

    ## T-test
    p_val_cont_CD = t_test(control_Score['Score'], CD_Score['Score'], p_val_var_Cont_CD )
    p_val_cont_CS = t_test(control_Score['Score'], CS_Score['Score'], p_val_var_Cont_CS )
    p_val_CS_CD = t_test(CS_Score['Score'], CD_Score['Score'], p_val_var_CS_CD )

    CSvsCD = data1[data1['Case'] > 0]
    ConvsCS = data1[data1['Case'] < 2]
    ConvsCD = data1[data1['Case'] != 1]
    
    #ConvsCD
    min_num , max_num = get_max_min(ConvsCD)
    ConvsCD['log'] = ConvsCD.Case.map({ min_num:0 , max_num:1})
    log_reg = (smf.logit(formula = 'log ~ Score' , data=ConvsCD).fit(maxiter=1000))
    p_val = (log_reg.pvalues.Score)
    label = str("Score Cont vs CD t-test p=" + str(format(p_val_cont_CD,'.2e'))+ " Log Reg p=" + str(format(p_val,'.2e')) )
    ConvsCD.groupby('Group').Score.plot.hist(alpha=0.5,  grid=True, legend=True, title=label  )
    plt.xlabel("Scaled Score Score")
    plt.legend(prop={'size':10})
    plt.tight_layout()
    output_name = (output+ name + "_ConvsCD.pdf")
    plt.savefig(output_name)
    plt.close()
    
    #ConvsCS
    log_reg = (smf.logit(formula = 'Case ~ Score' , data=ConvsCS).fit(maxiter=1000))
    p_val = (log_reg.pvalues.Score)
    label = str("Score Cont vs CS t-test p=" + str(format(p_val_cont_CS,'.2e')) + " Log Reg p=" + str(format(p_val,'.2e')) )
    ConvsCS.groupby('Group').Score.plot.hist(alpha=0.5, grid=True, legend=True, title=label)
    output_name = (output+ name + "_ConvsCS.pdf")
    plt.xlabel("Scaled Score Score")
    plt.legend(prop={'size':10})
    plt.tight_layout()
    plt.savefig(output_name)
    plt.close()
    
    ##CDvsCS
    min_num , max_num = get_max_min(CSvsCD) 
    CSvsCD['log'] = CSvsCD.Case.map({ min_num:0 , max_num:1})
    log_reg = (smf.logit(formula = 'log ~ Score' , data=CSvsCD).fit(maxiter=1000))
    p_val = (log_reg.pvalues.Score)
    label = str("Score CD vs CS t-test p=" + str(format(p_val_CS_CD,'.2e'))+ " Log Reg p=" + str(format(p_val,'.2e')))
    CSvsCD.groupby('Group').Score.plot.hist( alpha=0.5, grid=True, legend=True, title=label,)
    output_name = (output+ name + "_CSvsCD.pdf")
    plt.xlabel("Scaled Score Score")
    plt.legend(prop={'size':10})
    plt.tight_layout()
    plt.savefig(output_name)
    plt.close()
    
    ### convert CD and CS to 1 for log reg
    data1['log'] = data1.Case.map({ 0:0 , 1:1 , 2:1})
    log_reg = (smf.logit(formula = 'log ~ Score' , data=data1).fit(maxiter=1000))
    p_val = (log_reg.pvalues.Score)
    label = str("Score Case vs Control Log Reg p=" + str(format(p_val,'.2e') ))
    data1.groupby('Group').Score.plot.hist(alpha=0.5, grid=True, legend=True, title=label)
    output_name = (output+ name + "_group_wise.pdf")
    plt.xlabel("Scaled Score Score")
    plt.legend(prop={'size':10})
    plt.tight_layout()
    plt.savefig(output_name)
    plt.close()

    
if __name__ == '__main__':

    file , label , output = check_arg(sys.argv[1:])
    main(file , label, output)