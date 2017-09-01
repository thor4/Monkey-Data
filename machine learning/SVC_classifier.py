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
svm_clf1.fit(X_train6drfs, y_train6dr)
svm_clf2 = svm.SVC()
svm_clf2.fit(X_train8adfs, y_train8ad)
svm_clf3 = svm.SVC()
svm_clf3.fit(X_train8bfs, y_train8b)
svm_clf4 = svm.SVC()
svm_clf4.fit(X_train9lfs, y_train9l)
svm_clf5 = svm.SVC()
svm_clf5.fit(X_traindpfcfs, y_traindpfc)
svm_clf6 = svm.SVC()
svm_clf6.fit(X_trainlipfs, y_trainlip)
svm_clf7 = svm.SVC()
svm_clf7.fit(X_trainpefs, y_trainpe)
svm_clf8 = svm.SVC()
svm_clf8.fit(X_trainpecfs, y_trainpec)
svm_clf9 = svm.SVC()
svm_clf9.fit(X_trainpgfs, y_trainpg)
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
accuracySVM6drfs = cross_val_score(svm_clf1, X_train6drfs, y_train6dr, cv=5, scoring="accuracy")
accuracySVM8adfs = cross_val_score(svm_clf2, X_train8adfs, y_train8ad, cv=5, scoring="accuracy")
accuracySVM8bfs = cross_val_score(svm_clf3, X_train8bfs, y_train8b, cv=5, scoring="accuracy")
accuracySVM9lfs = cross_val_score(svm_clf4, X_train9lfs, y_train9l, cv=5, scoring="accuracy")
accuracySVMdpfcfs = cross_val_score(svm_clf5, X_traindpfcfs, y_traindpfc, cv=5, scoring="accuracy")
accuracySVMlipfs = cross_val_score(svm_clf6, X_trainlipfs, y_trainlip, cv=5, scoring="accuracy")
accuracySVMpefs = cross_val_score(svm_clf7, X_trainpefs, y_trainpe, cv=5, scoring="accuracy")
accuracySVMpecfs = cross_val_score(svm_clf8, X_trainpecfs, y_trainpec, cv=5, scoring="accuracy")
accuracySVMpgfs = cross_val_score(svm_clf9, X_trainpgfs, y_trainpg, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))