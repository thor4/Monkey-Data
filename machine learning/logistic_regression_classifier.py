#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 18:06:00 2017

@author: bconkli4
"""
import h5py
import matplotlib.pyplot as plt
import numpy as np
from sklearn.externals import joblib
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
import time

#load data
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-epochs.h5','r')
X1, y1 = h5f['early-data'][:], h5f['early-response'][:]
X2, y2 = h5f['mid-data'][:], h5f['mid-response'][:]
X3, y3 = h5f['late-data'][:], h5f['late-response'][:]
h5f.close()

#split 80/20 train/test sets, train_test_split shuffles by default and stratify
#to ensure biased class ratio is maintained
X_train1, X_test1, y_train1, y_test1 = train_test_split(X1, y1, test_size=0.2, random_state=42, stratify=y1)
X_train2, X_test2, y_train2, y_test2 = train_test_split(X2, y2, test_size=0.2, random_state=42, stratify=y2)
X_train3, X_test3, y_train3, y_test3 = train_test_split(X3, y3, test_size=0.2, random_state=42, stratify=y3)

#logistic regression classifier
start=time.time()
log_reg = LogisticRegression()
log_reg.fit(X_trainfs, y_train)
print("--- %s seconds ---" % (time.time() - start))

#logistic regression classifier standardized data
start=time.time()
log_reg = LogisticRegression()
log_reg.fit(stdx, y_train)
#param_grid = {'C': [0.001, 0.01, 0.1, 1, 10, 100] }
#log_reg = GridSearchCV(LogisticRegression(penalty='l2'), param_grid)

print("--- %s seconds ---" % (time.time() - start))


#logistic regression classifier epochs
start=time.time()
log_reg1 = LogisticRegression()
log_reg1.fit(X_train6drfs, y_train6dr)
log_reg2 = LogisticRegression()
log_reg2.fit(X_train8adfs, y_train8ad)
log_reg3 = LogisticRegression()
log_reg3.fit(X_train8bfs, y_train8b)
log_reg4 = LogisticRegression()
log_reg4.fit(X_train9lfs, y_train9l)
log_reg5 = LogisticRegression()
log_reg5.fit(X_traindpfcfs, y_traindpfc)
log_reg6 = LogisticRegression()
log_reg6.fit(X_trainlipfs, y_trainlip)
log_reg7 = LogisticRegression()
log_reg7.fit(X_trainpefs, y_trainpe)
log_reg8 = LogisticRegression()
log_reg8.fit(X_trainpecfs, y_trainpec)
log_reg9 = LogisticRegression()
log_reg9.fit(X_trainpgfs, y_trainpg)
print("--- %s seconds ---" % (time.time() - start))

#obtain accuracy
start=time.time()
accuracyLogRfs = cross_val_score(log_reg, X_trainfs, y_train, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))

#obtain accuracy standardized data
start=time.time()
accuracyLogRstd = cross_val_score(log_reg, stdx, y_train, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))

#accuracy areas
start=time.time()
accuracyLR6drfs = cross_val_score(log_reg1, X_train6drfs, y_train6dr, cv=5, scoring="accuracy")
accuracyLR8adfs = cross_val_score(log_reg2, X_train8adfs, y_train8ad, cv=5, scoring="accuracy")
accuracyLR8bfs = cross_val_score(log_reg3, X_train8bfs, y_train8b, cv=5, scoring="accuracy")
accuracyLR9lfs = cross_val_score(log_reg4, X_train9lfs, y_train9l, cv=5, scoring="accuracy")
accuracyLRdpfcfs = cross_val_score(log_reg1, X_traindpfcfs, y_traindpfc, cv=5, scoring="accuracy")
accuracyLRlipfs = cross_val_score(log_reg2, X_trainlipfs, y_trainlip, cv=5, scoring="accuracy")
accuracyLRpefs = cross_val_score(log_reg3, X_trainpefs, y_trainpe, cv=5, scoring="accuracy")
accuracyLRpecfs = cross_val_score(log_reg4, X_trainpecfs, y_trainpec, cv=5, scoring="accuracy")
accuracyLRpgfs = cross_val_score(log_reg5, X_trainpgfs, y_trainpg, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))

#hyperparameter tuning
start=time.time()
searchCV = LogisticRegressionCV(
    Cs=list(np.power(10.0, np.arange(-10, 10)))
    ,penalty='l2'
    ,scoring='accuracy'
    ,random_state=42
    ,fit_intercept=True
    ,solver='sag'
    ,max_iter=10000
)
searchCV.fit(X_trainpefs, y_trainpe)
print("--- %s seconds ---" % (time.time() - start))

print ('Max accuracy score (training):', searchCV.scores_[1].max())
accuracyLRpefsTest = searchCV.score(X_testpefs, y_testpe)
print("[INFO] randomized search accuracy: {:.2f}%".format(accuracyLRpefsTest * 100))
LRpeC = searchCV.C_; #best C value parameter

#save model
joblib.dump(searchCV, '/mnt/ceph/home/bconkli4/Documents/data/ml/models/peLR.pkl') 

#load model
clf = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/lipLR.pkl') 