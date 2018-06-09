#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import statsmodels.api as sm 
import statsmodels.formula.api as smf
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cross_validation import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt


data1 = pd.read_csv("ASRB_train_PRS_0.48_SEX.txt", sep='\t', index_col = 0 , header = 0 )
data1['Case'] = np.where(data1['Group']=='Control' , 0 , 1)
y = data1['Case']
X = pd.DataFrame()
scaler = StandardScaler()
data1[['PRS']] = scaler.fit_transform(data1[['PRS']])
X['PRS'] = data1['PRS']
X['Sex'] = data1['Sex']
#['CNV'] = data1['Score']
X_train , X_test , y_train , y_test = train_test_split(X , y , test_size = 0.8 , random_state=42)
model = LogisticRegression(penalty='l2' , C=1)
model.fit(X_train, y_train) 

print ("log is %2.2f" % accuracy_score(y_test,model.predict(X_test)))

logit_roc_auc = roc_auc_score(y_test, model.predict(X_test))

print ("AUC = %2.2f" % logit_roc_auc)


print (classification_report(y_test, model.predict(X_test)))


fpr , tpr , thresholds = roc_curve(y_test, model.predict_proba(X_test)[:,1])

plt.figure()
plt.plot(fpr , tpr , label='ROC curve (area = %0.2f)' % logit_roc_auc)
plt.plot([0,1], [0,1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic PGRS')
plt.legend(loc="lower right")
plt.show
#plt.savefig("ASRB_0.48_CNV_0.9.pdf")