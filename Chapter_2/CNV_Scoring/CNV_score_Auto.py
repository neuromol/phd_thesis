#!/usr/bin/env python3
# -*- coding: utf-8 -*-


### Get best CNV_SCORE


## importing libraries
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
from decimal import * 
import math
from sklearn.cross_validation import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve

#### Input Parameteres 
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Basic script to add up the Exact ODDs Ratio from the PGC CNV working group - Josh Atkins c3114203@uon.edu.au ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Clean_CNV File',required='True', default='test.cnv')
    parser.add_argument('-m', '--master', help='Master score file with CHR, BP, OR, P_vale' ,required='True')
    parser.add_argument('-s', '--start', help='Starting P_value threshold, Default 0.01 ', default='0.01')
    parser.add_argument('-f', '--finish', help='Finishing P_value threshold, Default 0.7 ', default='0.7')
    parser.add_argument('-n', '--increment', help='Increments of P_value threshold  DEFAULT: 0.01 ', default='0.01')
    parser.add_argument('-p', '--phenotype', help='Phenotype file ID, Group ', required='True')
    parser.add_argument('-o', '--output', help='Output Directory  DEFAULT: ./results ', default='./results')
    parser.add_argument('-l', '--label', help='Label for graphing DEFAULT: Test ', default='Test')
    
    results = parser.parse_args(args)
    return (results.input , results.master , results.start, results.finish , results.increment , results.phenotype , results.output ,results.label)

### Functions 

def make_score(master_file , p_value):
    score = master_file[master_file.P_value.apply(lambda x: x <= p_value)]
    return score 
    
def get_score(row, score_file):
    chr = row['CHR']
    bp1= row['BP1']
    bp2 = row['BP2']
    search = score_file[(score_file['CHR'] == chr) & (score_file['BP'] > bp1) & (score_file['BP'] < bp2)]
    total_score = search.OR.sum()
    return (total_score )

def increase_p(current_P, increment , finish):
    new_p = ( Decimal(current_P) + Decimal(increment) )
    
    if new_p >= finish : 
        return 0 
    else:
        return new_p 

def total_score(row , CNV_file , test ):
    ID = row['ID']
    total = CNV_file[(CNV_file['ID']==ID)]
    score = total[test].sum()
    return score

def get_phenotype(row , pheno_file):
    IDZ = row['ID']
    if ((pheno_file['ID'] == IDZ).any()) :  
        search = pheno_file[(pheno_file['ID'] == IDZ)]
        group = search.Group.values[0]
    else: 
        group = 0 
    return group    


def do_stats(df , score):
    test = score
    formula1 = ("bin1 ~ " + str(score) )
    log_reg = (smf.logit(formula = formula1 , data=df ).fit(maxiter=1000, disp=0))
    p_val = (log_reg.pvalues)
    return p_val 



def do_roc(df, score):
    
    y = df['bin1']
    X = pd.DataFrame()
    X['score'] = df[[score]]
    
    X_train , X_test , y_train , y_test = train_test_split(X , y , test_size = 0.2 , random_state=42)
    
    model = LogisticRegression(penalty='l2' , C=1)
    model.fit(X_train, y_train) 
    
    
    logit_roc_auc = Decimal(roc_auc_score(y_test, model.predict(X_test)) )
    fpr , tpr , thresholds = roc_curve(y_test, model.predict_proba(X_test)[:,1])

#    plt.figure()
#    plt.plot(fpr , tpr , label='ROC curve (area = %0.2f)' % logit_roc_auc)
#    plt.plot([0,1], [0,1], 'k--')
#    plt.xlim([0.0, 1.0])
#    plt.ylim([0.0, 1.05])
#    plt.xlabel('False Positive Rate')
#    plt.ylabel('True Positive Rate')
#    plt.title('Receiver Operating Characteristic PRS')
#    plt.legend(loc="lower right")
#    plt.show
#    
    
    return logit_roc_auc
    


def main(input_file , master_file , start , finish , increment , phenotype_file , output , label): 
    ### dealing with PennCNV output
    CNV_file=pd.read_csv(input_file,delimiter=r"\s+", header=None)
    CNV_file.columns=[ 'pos' , 'SNPs' , 'Size' , 'State' , 'ID' , 'startSNP' , 'endSNP' , 'Confidence_score']
	###uncomment line below if using columns.pl as it adds "sample.ID" this removes the sample.
    #CNV_file['ID'] = CNV_file['ID'].map(lambda x: x.split('.')[1]).astype(str)
    CNV_file['ID'] = CNV_file['ID'].astype(str)
    CNV_file['CHR'] = CNV_file['pos'].map(lambda x: x.split(':')[0])
    CNV_file['BP'] = CNV_file['pos'].map(lambda x: x.split(':')[1])
    CNV_file['BP1'] = CNV_file['BP'].map(lambda x: x.split('-')[0]).astype(int)
    CNV_file['BP2'] = CNV_file['BP'].map(lambda x: x.split('-')[1]).astype(int)
    CNV_file.drop(['pos' , 'SNPs', 'Size' , 'State' , 'startSNP' , 'endSNP' , 'Confidence_score', 'BP'], inplace =True, axis=1)
    
    
    ### dealing with Master_file 
    master_score_file = pd.read_csv(master_file, delimiter="\t") 
    master_score_file[['BP']] = master_score_file[['BP']].astype(int)
  
    ## Phenotype file 
    Pheno_file = pd.read_csv(phenotype_file , delimiter="\t")
    Pheno_file[['ID']] = Pheno_file[['ID']].astype(str)
    print ("cleaning up the input file and removing IDs not in the phenotype file " )
    ## getting rid of IDs we dont want  
    CNV_file['Group'] = CNV_file.apply(lambda row: get_phenotype(row, Pheno_file),axis=1) 
    CNV_file = CNV_file[CNV_file.Group != 0]
    #print (CNV_file)
    
    print (" " )
    print ("starting scoring -- this can take a while" )
    print (" ")
 ##working   
    current_P = Decimal(start) 
    best_p = pd.DataFrame()
    #best_p.columns = ['testID' , 'P_val' ]
    finish = Decimal(finish)
    increment = Decimal(increment)
    count = 1
    while current_P !=0 : 
        print ("Scoring P-value : " + str(current_P))
        test = make_score(master_score_file,current_P)
        name = ("test_" + str(count))
        CNV_file[name] = CNV_file.apply(lambda row: get_score(row, test),axis=1)        
        best_p = best_p.append({'testID' : name , 'P_val_cut_off' : current_P }, ignore_index=True )
        current_P = increase_p(current_P, increment , finish)
        
        count = count + 1

    print (" " )
    print ("Done scoring ... now generating each individual's score for each test" )
    ID_array = CNV_file.ID.unique()
    columns = ['ID']
    scores = pd.DataFrame(ID_array, columns=columns)

    first_test = 1 
    while first_test < count: 
        test_name = ('test_' + str(first_test))        
        scores[test_name] = scores.apply(lambda row: total_score(row, CNV_file, test_name),axis=1)
        first_test = first_test + 1 
    
    scores['Group'] = scores.apply(lambda row: get_phenotype(row, Pheno_file),axis=1) 
    
    counts = scores['Group'].value_counts()
    counts = counts.to_string(header=None )
    
    total_counts = scores['Group'].count()
    
    print (' ')

    print ('Total number of samples : ' + str(total_counts)) 
    print ('Which is made up of : ')
    print (counts)
    print (' ' ) 
    print ("performing logistic regression between groups") 
    
    scaler = StandardScaler()
    
    second_test = 1 
    while second_test < count: 
        test_name = ('test_' + str(second_test)) 
        scores[[test_name]] = scaler.fit_transform(scores[[test_name]])
        second_test = second_test + 1 
    
    
    
    dummies = pd.get_dummies(scores['Group'])
    dummies.columns=["bin1" , "bin2" ]

    scores = scores.join(dummies)
   
 
    P_results = pd.DataFrame()
    third_test = 1 
    while third_test < count: 
        test_name = ('test_' + str(third_test)) 
        
        results = do_stats(scores, test_name) 
        
        auc = do_roc(scores, test_name)
        
        
        p_val = (results[1])
        log = -math.log(p_val,10)
        P_results = P_results.append({'testID' : test_name , 'P_val' : p_val , '-log_P' : log , 'AUC' : auc}, ignore_index=True )  
        third_test = third_test + 1 
      
    
    results =pd.merge(P_results, best_p, on=['testID'])
    results['rank']= results['-log_P'].rank(ascending=False)
    
    best_P_cut_off = results.loc[results['rank'] == 1]['P_val_cut_off'].values
    
    best_P_cut_off = best_P_cut_off[0]
    
    best_P_result = results.loc[results['rank'] == 1]['P_val'].values
    
    best_P_result = best_P_result[0]
    
    best_auc = results.loc[results['rank'] == 1]['AUC'].values
    
    best_auc = best_auc[0]

    
    best_test = results.loc[results['rank'] == 1]['testID'].values
    best_test = best_test[0]
    
    col_list = ['ID' , best_test , "Group" ]

    output_df = scores[col_list]    
    
    print (" ")
    print ("The best P-value threshold for prediction is : " + str(best_P_cut_off))  
    print ("With a logistric Regression P-value of : " + str(best_P_result) )
    print (" " )

    output_name = (output + "_" + label + "_" + str(best_P_cut_off) + ".txt")
    results_output_name = ( output + "_" + label + "_" + "summary.txt")
    
    output_df.to_csv(output_name, sep='\t', index=False)
    results.to_csv(results_output_name, sep='\t', index=False)

if __name__ == '__main__':
    print ("___________________________________________________________________________________________________" )
    print (" " )
    print ("CNV_scoring python script to find the best P-value threshold for prediction by adding up odds ratios of CNV summary stats to individual datasets " )
    print ("   Josh Atkins  -- c3114203@uon.edu.au " )
    print (" ")
    print ("___________________________________________________________________________________________________" )
    print (" ")
    input_file , master_file , start , finish , increment , phenotype_file , output , label = check_arg(sys.argv[1:])
    main(input_file , master_file , start , finish , increment , phenotype_file , output , label)






