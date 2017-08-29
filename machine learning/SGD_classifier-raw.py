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
import pandas as pd
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
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-spec-fft-0_20_filtered.h5','r')
X, y = h5f['data'][:], h5f['response'][:]
h5f.close()

#filtered fft power by region
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-spec-fft-0_20_filtered_region.h5','r')
Xf, yf = h5f['dataF'][:], h5f['responseF'][:]
Xp, yp = h5f['dataP'][:], h5f['responseP'][:]
h5f.close()

#filtered raw data by area
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-0_20_filtered_base_norm-area-dpfc_lip.h5','r')
Xdpfc, ydpfc = h5f['datadpfc'][:], h5f['responsedpfc'][:]
Xlip, ylip = h5f['datalip'][:], h5f['responselip'][:]
h5f.close()

#standardize data to improve model
scalerf = StandardScaler()
scalerp = StandardScaler()
scalerf = scalerf.fit(X_trainf)
scalerp = scalerp.fit(X_trainp)
X_trainfstd = scalerf.transform(X_trainf)
X_trainpstd = scalerp.transform(X_trainp)
stdx = scalerx.transform(X_train)
stdxfs = scalerxfs.transform(X_train3fs)
print('Mean: %d, StandardDeviation: %d' % (np.mean(X_trainfstd), sqrt(np.var(X_trainpstd))))
# inverse transform brings back original data
inversed = scalerx.inverse_transform(stdx)

#split 80/20 train/test sets, train_test_split shuffles by default and stratify
#to ensure biased class ratio is maintained
X_train1, X_test1, y_train1, y_test1 = train_test_split(X1, y1, test_size=0.2, random_state=42, stratify=y1)
X_train2, X_test2, y_train2, y_test2 = train_test_split(X2, y2, test_size=0.2, random_state=42, stratify=y2)
X_train3, X_test3, y_train3, y_test3 = train_test_split(X3, y3, test_size=0.2, random_state=42, stratify=y3)

#chronux split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

#filtered fft power by region split
X_trainf, X_testf, y_trainf, y_testf = train_test_split(Xf, yf, test_size=0.2, random_state=42, stratify=yf)
X_trainp, X_testp, y_trainp, y_testp = train_test_split(Xp, yp, test_size=0.2, random_state=42, stratify=yp)

#filtered raw data by area split
X_traindpfc, X_testdpfc, y_traindpfc, y_testdpfc = train_test_split(Xdpfc, ydpfc, test_size=0.2, random_state=42, stratify=ydpfc)
X_trainlip, X_testlip, y_trainlip, y_testlip = train_test_split(Xlip, ylip, test_size=0.2, random_state=42, stratify=ylip)

#feature selection
#maximal information coefficient, standardization of training data yields the same scores
selector = SelectKBest(mutual_info_classif, k=250).fit(X_trainlip, y_trainlip)
scores = selector.scores_
#extract & plot which of the 250ms timepoints are most important
df = pd.DataFrame(scores);
top250 = df.nlargest(250,0);
topidx = np.sort(top250.index.values);
plt.figure()
plt.hist(topidx, bins=[0,100,200,300,400,500,600,700,810], facecolor='green', alpha=0.75)
plt.xlabel('Time (ms)')
plt.title('Histogram of Mutual Information Score for Each Time Point')
plt.grid(True)
plt.axis('tight')
#scores = -np.log10(selector.scores_)
X_trainlipfs = selector.transform(X_trainlip)

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

#train filtered fft by region classifier
sgd_clf1 = SGDClassifier(random_state=42)
sgd_clf1.fit(X_trainfstd, y_trainf)
sgd_clf2 = SGDClassifier(random_state=42)
sgd_clf2.fit(X_trainpstd, y_trainp)

#train filtered raw data by area classifier
sgd_clf1 = SGDClassifier(random_state=42)
sgd_clf1.fit(X_traindpfc, y_traindpfc)
sgd_clf2 = SGDClassifier(random_state=42)
sgd_clf2.fit(X_trainlipfs, y_trainlip)

#train SGD classifier chronux
sgd_clf = SGDClassifier(random_state=42)
sgd_clf.fit(X_train, y_train)

#train SGD classifier standardized
sgd_clf = SGDClassifier(random_state=42)
sgd_clf.fit(stdx, y_train)

accuracySGD = cross_val_score(sgd_clf, X_train, y_train, cv=5, scoring="accuracy")
accuracySGDstd = cross_val_score(sgd_clf, stdx, y_train, cv=5, scoring="accuracy")
accuracySGD1 = cross_val_score(sgd_clf1, X_train1, y_train1, cv=5, scoring="accuracy")
accuracySGD2 = cross_val_score(sgd_clf2, X_train2, y_train2, cv=5, scoring="accuracy")
accuracySGD3 = cross_val_score(sgd_clf3, X_train3, y_train3, cv=5, scoring="accuracy")

#accuracy filtered fft by region
accuracySGDfstd = cross_val_score(sgd_clf1, X_trainfstd, y_trainf, cv=5, scoring="accuracy")
accuracySGDpstd = cross_val_score(sgd_clf2, X_trainpstd, y_trainp, cv=5, scoring="accuracy")

#accuracy filtered raw data by area
accuracySGDdpfc = cross_val_score(sgd_clf1, X_traindpfc, y_traindpfc, cv=5, scoring="accuracy")
accuracySGDlipfs = cross_val_score(sgd_clf2, X_trainlipfs, y_trainlip, cv=5, scoring="accuracy")


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