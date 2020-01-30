#python ml_vcf_svc.py ann_matrix_cor.csv 0 # non normalized
#python ml_vcf_svc.py ann_matrix_cor.csv normalized? test_split?# normalized
#python ml_vcf_svc.py ann_matrix_cor.csv 1 0.2# normalized
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
normalize=sys.argv[2]

test_split=float(sys.argv[3])

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
pd.value_counts(df['TRUTH'])
df.describe().transpose()

#Training vs test
X = df.iloc[:,7:-1]
#X = X.drop(['P_CONTAM'],axis=1)
y = df['TRUTH']

#transform POP_AF and P_CONTAM result to log
X['POP_AF'] = [math.log(i,10) for i in X['POP_AF']]
X['P_CONTAM'] = [math.log(i,10) for i in (X['P_CONTAM']+1e-100)]

#make balanced  0 and 1 set with random downsize the 0
X_1 = X.loc[y==1]
X_0 = X.loc[y==0]
#random.sample(np.ndarray.tolist(np.linspace(1,10,10)),k=3)
#random.sample(X_1.index.tolist(),k=3)
#X.loc[random.sample(X_0.index.tolist(),k=X_1.shape[0])]
new_index=random.sample(X_0.index.tolist(),k=X_1.shape[0])+X_1.index.tolist()
X_temp=X.loc[new_index]
y_temp=y[new_index]

X_train, X_test, y_train, y_test = train_test_split(X_temp,y_temp,test_size=test_split,random_state=327)

print('Training data size: (%i,%i)' % X_train.shape)
print('Testing data size: (%i,%i)' % X_test.shape)

if normalize is None:
    print("normalize not define, set as default non-normalize")
    name_comm = ''
else:
    if normalize == '0':
        print("do not normalize for SVC")
        name_comm = ''
        X_train_a = X_train
    else:
        print("normalize for SVC, scale to [0,1]")
        scaler = MinMaxScaler()  # Default behavior is to scale to [0,1]
        #scaler = StandardScaler(with_mean=True, with_std=True)
        scaler.fit(X_train)
        X_train_a = scaler.transform(X_train)
        X_test_a = scaler.transform(X_test)
        name_comm = 'scale.'
        joblib.dump(scaler,'scaler.sav')

X_train_c = np.asarray(X_train_a)
y_train_c = np.asarray(y_train).reshape(-1)
X_test_c = np.asarray(X_test_a)
y_test_c = np.asarray(y_test).reshape(-1)
#not use cca or pca
#model = SVC(C=1, kernel='rbf', gamma='auto', decision_function_shape='ovr', class_weight='balanced', cahce_size=1000)
#start = time.time()
#model.fit(X_train_c, y_train_c)
#print(name_comm+"Kernel Fit Time: {0.4f} s".format(time() - start))
#cca_model = CCA(n_components=2)
#cca_model.fit(X_train_c, y_train_c)
#X_train_c = CCA(n_components=2).fit(X_train_c, Y_train_c).transform(X_train_c)
#X_train_c = cca_model.transform(X_train_c)
#X_train_c = PCA(n_components=2).fit_transform(X_train_c)
#classif = OneVsRestClassifier(SVC(kernel='linear'))
#model_c = SVC(C=1, kernel='rbf', gamma='auto', decision_function_shape='ovr', class_weight='balanced')
#model_c.fit(X_train_c, y_train_c)
#zero_class = np.where(Y_train_c[0])
#one_class = np.where(Y_train_c[1])


#intercept = regression_model.intercept_
#coef = pd.DataFrame(regression_model.coef_.transpose(),
#                    index=X.columns,
#                    columns=['Coefficients'])

#from sklearn.metrics import accuracy_score
#X_test_c = np.asarray(X_test_a)
#y_test_c = np.asarray(y_test).reshape(-1)
#start = time.time()
#y_pred_c = model.predict(X_test_c)
#print(name_comm+"Kernel predict Time: {0.4f} s".format(time() - start))

#test_acc = accuracy_score(y_test_c,y_pred_c)*100

#print('The test set accuracy is %4.2f%%' % test_acc)

#from sklearn.metrics import confusion_matrix

#conf_matrix = confusion_matrix(y_test_c,y_pred_c, labels=[0,1])

#tn, fp, fn, tp = confusion_matrix(y_test,y_pred).ravel()
#print(conf_matrix)

#X_test_c = cca_model.transform(X_test_c)
#X_test_c = PCA(n_components=2).fit_transform(X_test_c)
#y_pred_c = model_c.predict(X_test_c)
#test_acc = accuracy_score(y_test_c,y_pred_c)*100

#print('The test set accuracy is %4.2f%%' % test_acc)

#from sklearn.metrics import confusion_matrix

#conf_matrix = confusion_matrix(y_test_c,y_pred_c, labels=[0,1])

#tn, fp, fn, tp = confusion_matrix(y_test,y_pred).ravel()
#print(conf_matrix)

#Y_pred = y_pred_c.tolist()
#X_test['ML']=Y_pred  
#X_test.to_csv('test_svc_cca.csv',sep=',')


# save the model to disk
#filename = 'finalized_SVC_model.sav'
#joblib.dump(model, filename)

#joblib.dump(cca_model,'finalized_CCA_model.sav')

#joblib.dump(model,'finalized_SVC_model.'+name_comm+'sav')
#joblib.dump(model,'finalized_SVC_model_cca.sav')

# some time later...
 
# load the model from disk
#loaded_cca = joblib.load(cca_model)
#loaded_model = joblib.load(svc_model)
#X_test_c = loaded_cca.transform(X_test_c)
#y_pred_c = loaded_model.predict(X_test_c)
##print(result)

print('#k-Nearest Neighbors (k-NN)')
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
#create new a knn model
knn = KNeighborsClassifier(algorithm='auto')
#create a dictionary of all values we want to test for n_neighbors
params_knn = {'n_neighbors': np.arange(1, 25)}
#use gridsearch to test all values for n_neighbors
knn_gs = GridSearchCV(knn, params_knn, cv=5)
#fit model to training data
knn_gs.fit(X_train_c, y_train_c)
#save best model
knn_best = knn_gs.best_estimator_
#check best n_neigbors value
print(knn_gs.best_params_)
print('knn: {}'.format(knn_best.score(X_test_c, y_test_c))) #0.9186547858580052
y_pred_c = knn_best.predict(X_test_c)
conf_matrix = confusion_matrix(y_test_c,y_pred_c, labels=[0,1])
print(conf_matrix)
joblib.dump(knn_best,'knn_model.'+name_comm+'sav')
Y_pred = y_pred_c.tolist()
X_test['Truth']=y_test
X_test['knn']=Y_pred  
pd.DataFrame(conf_matrix).to_csv(f1name_re+'_knn_'+'conf.tsv',sep='\t')



print('#Random Forest (allow warm start)')
from sklearn.ensemble import RandomForestClassifier
#create a new random forest classifier
rf = RandomForestClassifier(warm_start=False, class_weight="balanced")
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
print(conf_matrix)
joblib.dump(rf_best,'rf_model.'+name_comm+'sav')
Y_pred = y_pred_c.tolist()
X_test['Truth']=y_test
X_test['rf']=Y_pred  
pd.DataFrame(conf_matrix).to_csv(f1name_re+'_rf_'+'conf.tsv',sep='\t')


#logisticRegression
print('from sklearn.linear_model import LogisticRegression')
#create a new logistic regression model
params_log_reg = {'solver': ['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga']}
log_reg = LogisticRegression(class_weight='balanced')
log_reg_gs = GridSearchCV(log_reg,params_log_reg,cv=5)
#fit the model to the training data
log_reg_gs.fit(X_train, y_train)
log_reg_best = log_reg_gs.best_estimator_
print(log_reg_gs.best_params_)
print('log_reg: {}'.format(log_reg_best.score(X_test_c, y_test_c))) #0.5004311583788444
y_pred_c = log_reg_best.predict(X_test_c)
conf_matrix = confusion_matrix(y_test_c,y_pred_c, labels=[0,1])
print(conf_matrix)
joblib.dump(log_reg_best,'log_reg_model.'+name_comm+'sav')
Y_pred = y_pred_c.tolist()
#X_test['Truth']=y_test
X_test['log_reg']=Y_pred  
pd.DataFrame(conf_matrix).to_csv(f1name_re+'_log_reg_'+'conf.tsv',sep='\t')

#SVC rbf #slow and not good for small training set
print('from sklearn.svm import SVC')
svc_rbf = SVC(C=1, kernel='rbf', gamma='auto', decision_function_shape='ovr', class_weight='balanced', cache_size=1000)
params_svc_rbf = {'C': [0.5, 1, 1.5, 2]}
svc_rbf_gs = GridSearchCV(svc_rbf, params_svc_rbf, cv=5)
svc_rbf_gs.fit(X_train_c, y_train_c)
svc_rbf_best = svc_rbf_gs.best_estimator_
print(svc_rbf_gs.best_estimator_)
print('svc_rbf: {}'.format(svc_rbf_best.score(X_test_c, y_test_c))) #0.6630463858151228
y_pred_c = svc_rbf_best.predict(X_test_c)
conf_matrix = confusion_matrix(y_test_c,y_pred_c, labels=[0,1])
print(conf_matrix)
joblib.dump(svc_rbf_best,'svc_rbf_model.'+name_comm+'sav')
Y_pred = y_pred_c.tolist()
#X_test['Truth']=y_test
X_test['svc_rbf']=Y_pred  
pd.DataFrame(conf_matrix).to_csv(f1name_re+'_svc_rbf_'+'conf.tsv',sep='\t')

#test the three models with the test data and print their accuracy scores

print('#Voting Classigfier')
from sklearn.ensemble import VotingClassifier
#create a dictionary of our models
#estimators=[('knn', knn_best), ('rf', rf_best), ('log_reg', log_reg_best)]
estimators=[('knn', knn_best), ('rf', rf_best), ('svc_rbf', svc_rbf_best)]
#create our voting classifier, inputting our models
ensemble = VotingClassifier(estimators, voting='hard')
#fit model to training data
ensemble.fit(X_train_c, y_train_c)

#VotingClassifier(estimators=[('knn', KNeighborsClassifier(algorithm='auto', leaf_size=30, metric='minkowski',
#           metric_params=None, n_jobs=1, n_neighbors=20, p=2,
#           weights='uniform')), ('rf', RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
#            max_depth=None, ma...ty='l2', random_state=None, solver='liblinear', tol=0.0001,
#          verbose=0, warm_start=False))],
#         flatten_transform=None, n_jobs=1, voting='hard', weights=None)
#0.98067758119018889
#test our model on the test data
ensemble.score(X_test_c, y_test_c)
#print('voting: {}'.format(ensemble.score(X_test_c, y_test_c)))
y_pred_c = ensemble.predict(X_test_c)
conf_matrix = confusion_matrix(y_test_c,y_pred_c, labels=[0,1])

#tn, fp, fn, tp = confusion_matrix(y_test,y_pred).ravel()
print(conf_matrix)
y_pred_c = ensemble.predict(X_test_c)
conf_matrix = confusion_matrix(y_test_c,y_pred_c, labels=[0,1])
print(conf_matrix)
joblib.dump(ensemble,'ensemble_model.'+name_comm+'sav')
Y_pred = y_pred_c.tolist()
#X_test['Truth']=y_test
X_test['ensemble']=Y_pred  
pd.DataFrame(conf_matrix).to_csv(f1name_re+'_ensemble_'+'conf.tsv',sep='\t')

X_test.to_csv(f1name_re+'ml_best.csv',sep=',')