#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 04:10:36 2017

@author: bconkli4
"""

from sklearn import svm
import time

#SVC classifier
start = time.time()
svm_clf = svm.SVC()
svm_clf.fit(X_train, y_train)
print("--- %s seconds ---" % (time.time() - start))

#SVC classifier epochs
start = time.time()
svm_clf1 = svm.SVC()
svm_clf1.fit(X_train1, y_train1)
svm_clf2 = svm.SVC()
svm_clf2.fit(X_train2, y_train2)
svm_clf3 = svm.SVC()
svm_clf3.fit(X_train3, y_train3)
svm_clf4 = svm.SVC()
svm_clf4.fit(X_train3fs, y_train3)
print("--- %s seconds ---" % (time.time() - start))

start = time.time()
p = svm_clf.score(X_test, y_test)
print("--- %s seconds ---" % (time.time() - start))

#accuracy
start=time.time()
accuracySVM = cross_val_score(svm_clf, X_train, y_train, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))

#accuracy
start=time.time()
accuracySVM1 = cross_val_score(svm_clf, X_train1, y_train1, cv=5, scoring="accuracy")
accuracysvm2 = cross_val_score(svm_clf2, X_train2, y_train2, cv=5, scoring="accuracy")
accuracysvm3 = cross_val_score(svm_clf3, X_train3, y_train3, cv=5, scoring="accuracy")
accuracysvm3fs = cross_val_score(svm_clf4, X_train3fs, y_train3, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))