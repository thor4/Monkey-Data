#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 18:58:26 2017

@author: bconkli4
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
from sklearn.grid_search import RandomizedSearchCV
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
nbrs_clf1.fit(X_train6drfs, y_train6dr)
nbrs_clf2 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf2.fit(X_train8adfs, y_train8ad)
nbrs_clf3 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf3.fit(X_train8bfs, y_train8b)
nbrs_clf4 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf4.fit(X_train9lfs, y_train9l)
nbrs_clf5 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf5.fit(X_traindpfcfs, y_traindpfc)
nbrs_clf6 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf6.fit(X_trainlipfs, y_trainlip)
nbrs_clf7 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf7.fit(X_trainpecfs, y_trainpec)
nbrs_clf8 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf8.fit(X_trainpefs, y_trainpe)
nbrs_clf9 = KNeighborsClassifier(n_neighbors=1, algorithm='auto')
nbrs_clf9.fit(X_trainpgfs, y_trainpg)
print("--- %s seconds ---" % (time.time() - start))

#obtain accuracy
start = time.time()
accuracyNN = cross_val_score(nbrs_clf, X_train, y_train, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))
pred = nbrs_clf.predict(X_test3)
accnn = accuracy_score(y_test3, pred)

#accuracy epochs
start = time.time()
accuracyNN6drfs = cross_val_score(nbrs_clf1, X_train6drfs, y_train6dr, cv=5, scoring="accuracy")
accuracyNN8adfs = cross_val_score(nbrs_clf2, X_train8adfs, y_train8ad, cv=5, scoring="accuracy")
accuracyNN8bfs = cross_val_score(nbrs_clf3, X_train8bfs, y_train8b, cv=5, scoring="accuracy")
accuracyNN9lfs = cross_val_score(nbrs_clf4, X_train9lfs, y_train9l, cv=5, scoring="accuracy")
accuracyNNdpfcfs = cross_val_score(nbrs_clf1, X_traindpfcfs, y_traindpfc, cv=5, scoring="accuracy")
accuracyNNlipfs = cross_val_score(nbrs_clf2, X_trainlipfs, y_trainlip, cv=5, scoring="accuracy")
accuracyNNpecfs = cross_val_score(nbrs_clf3, X_trainpecfs, y_trainpec, cv=5, scoring="accuracy")
accuracyNNpefs = cross_val_score(nbrs_clf4, X_trainpefs, y_trainpe, cv=5, scoring="accuracy")
accuracyNNpgfs = cross_val_score(nbrs_clf5, X_trainpgfs, y_trainpg, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))

#hyperparameter tuning
params = {"n_neighbors": np.arange(1, 32, 2),
          "metric": ["euclidean", "minkowski", "chebyshev"]}
searchCVpeNN = RandomizedSearchCV(nbrs_clf8, params)
start = time.time()
searchCVpeNN.fit(X_trainpefs, y_trainpe)
# evaluate the best randomized searched model on the testing
# data
print("[INFO] randomized search took {:.2f} seconds".format(
	time.time() - start))
accuracyNNpefsTuned = searchCVpeNN.score(X_trainpefs, y_trainpe)
accuracyNNpefsTest = searchCVpeNN.score(X_testpefs, y_testpe)
print("[INFO] randomized search accuracy: {:.2f}%".format(accuracyNNpefsTest * 100))
print("[INFO] randomized search best parameters: {}".format(
	searchCVpeNN.best_params_))

#save model
joblib.dump(searchCVpeNN, '/mnt/ceph/home/bconkli4/Documents/data/ml/models/peNN.pkl') 