#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:19:01 2017

@author: dickie_ho
"""
from decimal import Decimal
import numpy as np
import pandas as pd
from pandas import DataFrame, Series
import statsmodels.formula.api as sm
from sklearn.linear_model import LinearRegression
import scipy, scipy.stats
import matplotlib.pyplot as plt 
from sklearn.preprocessing import StandardScaler

def convert_group(x): 
    #convert labels to cont vars 
    if x == str("Control"):
        return 0
    elif x == str("Spared Cognition"):
        return 1 
    elif x == str("Impaired cognition"):
        return 1


data1 = pd.read_csv("ASRB.txt", sep = "\t",  header = 0 )

data1['Case'] = data1['Group'].apply(convert_group)
scaler = StandardScaler()
data1[['PRS']] = scaler.fit_transform(data1[['PRS']])

CD = data1.loc[(data1['Case'] == 1) ]

#CD = data1


Y = CD['Score']
X = CD['PRS']
result = sm.OLS( Y, X ).fit()
print (result.summary() ) 
#ax = plt.subplot(111)

#CD.plot('PRS', 'Score', ax=ax , kind="scatter" )

#ax.set_xlim(-0.00043,0.000012)
r_val = Decimal(result.rsquared)
r_val = round(r_val, 3)
plt.scatter(X,Y , 1)
z = np.polyfit(X, Y, 1)
p = np.poly1d(z)
plt.plot(X,p(X),"r--")
plt.axis([-2.0003,4 , -2, 7])
plt.xlabel('PRS score scaled')
plt.ylabel("CNV score scaled")
plt.title('OLS Regression of PRS and CNV Score')
plt.text(2, 6, "R-square = " + str(r_val) )

plt.savefig("linReg.pdf")
plt.show()