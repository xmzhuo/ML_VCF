#python ml_vcf_rf_train_v2.py ann_matrix_cor.csv test_split(within 0 to 1;float)? balanced(TN:TP, 0 without adjustment; float)? warm_model (optional)
#python ml_vcf_rf_v1.py ann_matrix_cor.csv normalized([0,1])? test_split(within 0 to 1;float)? balanced(TN:TP, 0 without adjustment; float)? warm_model (optional)
#python ml_vcf_rf_v1.py ann_matrix_cor.csv 0.2 2
#python ml_vcf_rf_v1.py ann_matrix_cor.csv 0.2 2 warm_model.sav scaler.sav
#python ../ml_vcf_rf_v1.py TPF-18-112_KAM56A7_S3_001_twicefiltered.ann_matrix_cor.csv 0.2 0
#logistic regression

import pandas as pd
import numpy as np
import sys
import os
import re
import math
import time
import random
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
import numpy as np
import matplotlib.pyplot as plt

#from sklearn.datasets import make_multilabel_classification
#from sklearn.datasets import make_classification
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

#normalize=1 or 0
#normalize=sys.argv[2]

test_split=float(sys.argv[2])

sample_balance=float(sys.argv[3])

if len(sys.argv) > 4 :
    warm_model = sys.argv[4]
    #scaler_model = sys.argv[5]
    print(warm_model)
    #print(scaler_model)

f1name_re=re.sub('csv','',f1name)
print(f1name_re)

try:
    f1hand=open(f1name)
    print(f1name)
except:
    print('Input file 1 no exist')

df = pd.read_csv(f1name)
df.shape
df.head()
pd.value_counts(df['%INFO/TRUTH'])
df.describe().transpose()

#Training vs test
X = df.iloc[:,7:-1]
#X = X.drop(['P_CONTAM'],axis=1)
y = df['%INFO/TRUTH']

##transform POP_AF and P_CONTAM result to log
#X['POP_AF'] = [math.log(i,10) for i in X['POP_AF']]
#X['P_CONTAM'] = [math.log(i,10) for i in (X['P_CONTAM']+1e-100)]

if sample_balance > 0 :
    print("balance TN:TP ratio:",sample_balance)
    # make balanced  0 and 1 set with random downsize the 0
    X_1 = X.loc[y==1]
    X_0 = X.loc[y==0]
    k_value=int(X_1.shape[0]*sample_balance)
    new_index=random.sample(X_0.index.tolist(),k=k_value)+X_1.index.tolist()
    X_temp=X.loc[new_index]
    y_temp=y[new_index]
else :
    X_temp=X
    y_temp=y

X_train, X_test, y_train, y_test = train_test_split(X_temp,y_temp,test_size=test_split,random_state=327)

print('Training data size: (%i,%i)' % X_train.shape)
print('Testing data size: (%i,%i)' % X_test.shape)
          
X_train_c = np.asarray(X_train)
y_train_c = np.asarray(y_train).reshape(-1)
X_test_c = np.asarray(X_test)
y_test_c = np.asarray(y_test).reshape(-1)

# some time later...
#pipeline = Pipeline([('scaler', StandardScaler()), ('classifier', DecisionTreeClassifier())])
#pipeline.fit(X_train, y_train)
#predictions = pipeline.predict(X_test)
 
# load the model from disk
#loaded_cca = joblib.load(cca_model)
#loaded_model = joblib.load(svc_model)
#X_test_c = loaded_cca.transform(X_test_c)
#y_pred_c = loaded_model.predict(X_test_c)
##print(result)

print('#Random Forest (allow warm start)')
if len(sys.argv) > 4 :
    print('warm start with',warm_model)
    rf_best = joblib.load(warm_model)
    print(rf_best.get_params)
    rf_best.n_estimators += 100
    #rf_best.set_params(warm_start=True)
    print(rf_best.get_params)
    rf_best.fit(X_train_c, y_train_c)
else :
    print("create a new random forest classifier")
    rf = RandomForestClassifier(warm_start=True)
    #rf = RandomForestClassifier(warm_start=False, class_weight="balanced")
    #create a dictionary of all values we want to test for n_estimators
    params_rf = {'n_estimators': [10, 50, 100, 200]}
    #use gridsearch to test all values for n_estimators
    rf_gs = GridSearchCV(rf, params_rf, cv=5)
    #fit model to training data
    rf_gs.fit(X_train_c, y_train_c)
    #save best model
    rf_best = rf_gs.best_estimator_
    #check best n_estimators value
    print(rf_gs.best_params_)

print('rf: {}'.format(rf_best.score(X_test_c, y_test_c))) #0.9304397815464214
y_pred_c = rf_best.predict(X_test_c)
conf_matrix = confusion_matrix(y_test_c,y_pred_c, labels=[0,1])
tn, fp, fn, tp = confusion_matrix(y_test_c,y_pred_c).ravel()
rf_precision = tp/(tp+fp)
rf_sensitivity = tp/(tp+fn)
rf_auc = rf_precision*rf_sensitivity
print(conf_matrix)
print('Precision: {0}; Sensitivity: {1}; AUC: {2}'.format(rf_precision, rf_sensitivity, rf_auc))
time_stamp = time.strftime("%Y-%m-%d-%H-%M")
joblib.dump(rf_best,'vcf_rf_model.'+time_stamp+'.sav')
Y_pred = y_pred_c.tolist()
X_test['Truth']=y_test
X_test['rf']=Y_pred  
rf_conf = pd.DataFrame(conf_matrix)
rf_conf['PPV/TPR']=[rf_precision,rf_sensitivity]
rf_conf.to_csv(f1name_re+'_vcf_rf_'+time_stamp+"_"+'_conf.tsv',sep='\t')

X_test.to_csv(f1name_re+'_vcf_rf_'+time_stamp+'_testpred.csv',sep=',')