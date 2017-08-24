#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 23:38:12 2017

@author: bconkli4
"""

import h5py
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
from sklearn.externals import joblib
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, roc_curve
from sklearn.preprocessing import StandardScaler

h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-epochs.h5','r')
X1, y1 = h5f['early-data'][:], h5f['early-response'][:]
X2, y2 = h5f['mid-data'][:], h5f['mid-response'][:]
X3, y3 = h5f['late-data'][:], h5f['late-response'][:]
h5f.close()

#standardize data to improve model
scaler = StandardScaler()
scaler = scaler.fit(X_train3)
# standardization the dataset and print the first 5 rows
normalized = scaler.transform(X_train3)
print('Mean: %d, StandardDeviation: %d' % (np.mean(normalized), sqrt(np.var(normalized))))
# inverse transform brings back original data
inversed = scaler.inverse_transform(normalized)

#split 80/20 train/test sets, train_test_split shuffles by default and stratify
#to ensure biased class ratio is maintained
X_train1, X_test1, y_train1, y_test1 = train_test_split(X1, y1, test_size=0.2, random_state=42, stratify=y1)
X_train2, X_test2, y_train2, y_test2 = train_test_split(X2, y2, test_size=0.2, random_state=42, stratify=y2)
X_train3, X_test3, y_train3, y_test3 = train_test_split(X3, y3, test_size=0.2, random_state=42, stratify=y3)

#manual split test
X_train3, X_test3, y_train3, y_test3

#shuffle to ensure all cv folds will be similar. also some learning algorithms 
#are sensitive to order of training instances. they perform poorly if they get 
#many similar instances in a row

#train SGD classifier
sgd_clf1 = SGDClassifier(random_state=42)
sgd_clf1.fit(X_train1, y_train1)
sgd_clf2 = SGDClassifier(random_state=42)
sgd_clf2.fit(X_train2, y_train2)
sgd_clf3 = SGDClassifier(random_state=42)
sgd_clf3.fit(X3, y3)

accuracy1 = cross_val_score(sgd_clf1, X_train1, y_train1, cv=5, scoring="accuracy")
accuracy2 = cross_val_score(sgd_clf2, X_train2, y_train2, cv=5, scoring="accuracy")
accuracy3n = cross_val_score(sgd_clf3, normalized, y_train3, cv=5, scoring="accuracy")
accuracy3full = cross_val_score(sgd_clf3, X3, y3, cv=5, scoring="accuracy")

#normalize example

# Normalize time series data
from pandas import Series
from sklearn.preprocessing import MinMaxScaler
# load the dataset and print the first 5 rows
series = Series.from_csv('daily-minimum-temperatures-in-me.csv', header=0)
print(series.head())
# prepare data for normalization
values = series.values
values = values.reshape((len(values), 1))
# train the normalization
scaler = MinMaxScaler(feature_range=(0, 1))
scaler = scaler.fit(values)
print('Min: %f, Max: %f' % (scaler.data_min_, scaler.data_max_))
# normalize the dataset and print the first 5 rows
normalized = scaler.transform(values)
for i in range(5):
	print(normalized[i])
# inverse transform and print the first 5 rows
inversed = scaler.inverse_transform(normalized)
for i in range(5):
	print(inversed[i])