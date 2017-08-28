#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 18:58:26 2017

@author: bconkli4
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
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

#1 nearest neighbors classifier
start = time.time()
nbrs_clf = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf.fit(X_train, y_train)
print("--- %s seconds ---" % (time.time() - start))

#1 nearest neighbors classifier epochs
start = time.time()
nbrs_clf1 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf1.fit(X_train1, y_train1)
nbrs_clf2 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf2.fit(X_train2, y_train2)
nbrs_clf3 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf3.fit(X_train3, y_train3)
nbrs_clf4 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf4.fit(X_train3fs, y_train3)
print("--- %s seconds ---" % (time.time() - start))

#obtain accuracy
start = time.time()
accuracyNN = cross_val_score(nbrs_clf, X_train, y_train, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))
pred = nbrs_clf.predict(X_test3)
accnn = accuracy_score(y_test3, pred)

#accuracy epochs
start = time.time()
accuracy1NN = cross_val_score(nbrs_clf1, X_train1, y_train1, cv=5, scoring="accuracy")
accuracy2NN = cross_val_score(nbrs_clf2, X_train2, y_train2, cv=5, scoring="accuracy")
accuracy3NN = cross_val_score(nbrs_clf3, X_train3, y_train3, cv=5, scoring="accuracy")
accuracy3fsNN = cross_val_score(nbrs_clf4, X_train3fs, y_train3, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))