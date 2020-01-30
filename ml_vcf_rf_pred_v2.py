#python ml_vcf_rf_pred_v1.py ann_matrix_cor.csv ML_model 
#python ml_vcf_rf_pred.py_v1 ann_matrix_cor.csv rf_model.scale.sav 
#python ../ml_vcf_rf_pred.py TPF-18-270_tu_2_KAM56A86_S7_001_twicefiltered.ann_matrix_cor.csv Mutect2_rf_model.2019-08-20-16-59.sav 


import pandas as pd
import numpy as np
import sys
import os
import re
import math
import time
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
import numpy as np
import matplotlib.pyplot as plt

from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import CCA
from sklearn.externals import joblib

from sklearn.metrics import confusion_matrix

from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

print ("This is the name of the script: ", sys.argv[0])
print ("Number of arguments: ", len(sys.argv))
print ("The arguments are: " , str(sys.argv))

f1name=sys.argv[1]

ml_model=sys.argv[2]

#decision tree model do not need scale
#scaler_model=sys.argv[3]

f1name_re=re.sub('csv','',f1name)
ml_model_re=re.sub('sav','',ml_model)
print(f1name_re)

try:
    f1hand=open(f1name)
    print(f1name)
except:
    print('Input file matrix_cor.csv no exist')

try:
    f2hand=open(ml_model)
    print(ml_model)
except:
    print('Input ML model no exist')

#try:
#    f3hand=open(scaler_model)
#    print(scaler_model)
#except:
#    print('Input scaler model no exist')


df = pd.read_csv(f1name)
df.shape
df.head()
#pd.value_counts(df['TRUTH'])
df.describe().transpose()

if '%INFO/TRUTH' in df.columns:
    X = df.iloc[:,7:-1]
    y = df['%INFO/TRUTH']
    y_test_c = np.asarray(y).reshape(-1)
else:
    X = df.iloc[:,7:]

##transform POP_AF and P_CONTAM result to log
#X['POP_AF'] = [math.log(i,10) for i in X['POP_AF']]
#X['P_CONTAM'] = [math.log(i,10) for i in (X['P_CONTAM']+1e-100)]

#scaler = joblib.load(scaler_model)
#scaler.set_params(with_mean=False)
#X_test_a = scaler.transform(X)
X_test_c = np.asarray(X)

loaded_model = joblib.load(ml_model)
y_pred_c = loaded_model.predict(X_test_c)

if '%INFO/TRUTH' in df:
    conf_matrix = confusion_matrix(y_test_c,y_pred_c, labels=[0,1])
    print(conf_matrix)
    tn, fp, fn, tp = confusion_matrix(y_test_c,y_pred_c).ravel()
    rf_precision = tp/(tp+fp)
    rf_sensitivity = tp/(tp+fn)
    rf_auc = rf_precision*rf_sensitivity
    print('Precision: {0}; Sensitivity: {1}; AUC: {2}'.format(rf_precision, rf_sensitivity, rf_auc))
    rf_conf = pd.DataFrame(conf_matrix)
    rf_conf['PPV/TPR']=[rf_precision,rf_sensitivity]
    rf_conf.to_csv(f1name_re+'_Mutect2_rf_'+ml_model_re+'conf.tsv',sep='\t')
    #pd.DataFrame(conf_matrix).to_csv(f1name_re+ml_model_re+'conf.tsv',sep='\t')
    #X_test['Truth']=y_test

Y_pred = y_pred_c.tolist()
df['rf']=Y_pred  
#df.to_csv(f1name_re+ml_model_re+'pred.tsv',sep='\t')
#temp = df[df['rf']==1]
#out_bed = temp[['CHROM','POS']]
out_bed = df[['%CHROM','%POS','rf']]
out_bed['POS1'] = out_bed['%POS'] - 1
#out_bed = out_bed[['CHROM','POS1','POS']]
out_bed = out_bed[['%CHROM','POS1','%POS','rf']]
out_bed.to_csv(f1name_re+ml_model_re+'pred.bed',sep='\t',header=False,index=False)