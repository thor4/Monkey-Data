#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 04:10:36 2017

@author: bconkli4
"""

from sklearn import svm

#SVC classifier
start = time.time()
svm_clf = svm.SVC()
svm_clf.fit(X_train3, y_train3)
print("--- %s seconds ---" % (time.time() - start))