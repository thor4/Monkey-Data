#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 18:06:00 2017

@author: bconkli4
"""
import h5py
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
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
log_reg.fit(X_train, y_train)
#param_grid = {'C': [0.001, 0.01, 0.1, 1, 10, 100] }
#log_reg = GridSearchCV(LogisticRegression(penalty='l2'), param_grid)

print("--- %s seconds ---" % (time.time() - start))

#obtain accuracy
start=time.time()
accuracyLogR = cross_val_score(log_reg, X_train, y_train, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))