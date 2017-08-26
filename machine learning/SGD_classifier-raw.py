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
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, roc_curve
from sklearn.preprocessing import StandardScaler

h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-epochs-0_20_filtered-base_norm_subsample.h5','r')
X1, y1 = h5f['early-d'][:], h5f['early-r'][:]
X2, y2 = h5f['mid-d'][:], h5f['mid-r'][:]
X3, y3 = h5f['late-d'][:], h5f['late-r'][:]
h5f.close()

#epochs stats
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-epochs-stats-0_20_filtered-base_norm_subsample.h5','r')
X1, y1 = h5f['early-d'][:], h5f['early-r'][:]
X2, y2 = h5f['mid-d'][:], h5f['mid-r'][:]
X3, y3 = h5f['late-d'][:], h5f['late-r'][:]
h5f.close()

#chronux
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-chronux.h5','r')
X, y = h5f['data'][:], h5f['response'][:]
h5f.close()

#raw sub sampled data
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-epochs-0_20_filtered-base_norm_subsample.h5','r')
X, y = h5f['data'][:], h5f['response'][:]
h5f.close()

#standardize data to improve model
scalerx = StandardScaler()
scalerxfs = StandardScaler()
scalerx = scalerx.fit(X_train3)
scalerxfs = scalerxfs.fit(X_train3fs)
# standardization the dataset and print the first 5 rows
stdx = scalerx.transform(X_train3)
stdxfs = scalerxfs.transform(X_train3fs)
print('Mean: %d, StandardDeviation: %d' % (np.mean(stdx), sqrt(np.var(stdx))))
# inverse transform brings back original data
inversed = scalerx.inverse_transform(stdx)

#split 80/20 train/test sets, train_test_split shuffles by default and stratify
#to ensure biased class ratio is maintained
X_train1, X_test1, y_train1, y_test1 = train_test_split(X1, y1, test_size=0.2, random_state=42, stratify=y1)
X_train2, X_test2, y_train2, y_test2 = train_test_split(X2, y2, test_size=0.2, random_state=42, stratify=y2)
X_train3, X_test3, y_train3, y_test3 = train_test_split(X3, y3, test_size=0.2, random_state=42, stratify=y3)

#chronux split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

#feature selection
#maximal information coefficient
selector = SelectKBest(mutual_info_classif, k=2).fit(X_train3, y_train3)
selector.scores_
scores = -np.log10(selector.scores_)
X_train3fs = selector.transform(X_train3)

mi1 = mutual_info_classif(X1, y1)
mi2 = mutual_info_classif(X2, y2)
mi3 = mutual_info_classif(X3, y3)

#shuffle to ensure all cv folds will be similar. also some learning algorithms 
#are sensitive to order of training instances. they perform poorly if they get 
#many similar instances in a row

#train SGD classifier
sgd_clf1 = SGDClassifier(random_state=42)
sgd_clf1.fit(X_train1, y_train1)
sgd_clf2 = SGDClassifier(random_state=42)
sgd_clf2.fit(X_train2, y_train2)
sgd_clf3 = SGDClassifier(random_state=42)
sgd_clf3.fit(X_train3, y_train3)
sgd_clf4 = SGDClassifier(random_state=42)
sgd_clf4.fit(stdx, y_train3)
sgd_clf5 = SGDClassifier(random_state=42)
sgd_clf5.fit(X_train3fs, y_train3)
sgd_clf6 = SGDClassifier(random_state=42)
sgd_clf6.fit(stdxfs, y_train3)

#train SGD classifier chronux
sgd_clf = SGDClassifier(random_state=42)
sgd_clf.fit(X_train, y_train)
sgd_clf.fit(normalx, y_train)

accuracy = cross_val_score(sgd_clf, X_train, y_train, cv=5, scoring="accuracy")
accuracyn = cross_val_score(sgd_clf, normalx, y_train, cv=5, scoring="accuracy")
accuracySGD1 = cross_val_score(sgd_clf1, X_train1, y_train1, cv=5, scoring="accuracy")
accuracySGD2 = cross_val_score(sgd_clf2, X_train2, y_train2, cv=5, scoring="accuracy")
accuracySGD3 = cross_val_score(sgd_clf3, X_train3, y_train3, cv=5, scoring="accuracy")
accuracySGD3std = cross_val_score(sgd_clf4, stdx, y_train3, cv=5, scoring="accuracy")
accuracySGD3fs = cross_val_score(sgd_clf5, X_train3fs, y_train3, cv=5, scoring="accuracy")
accuracySGD3fsstd = cross_val_score(sgd_clf6, stdxfs, y_train3, cv=5, scoring="accuracy")
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