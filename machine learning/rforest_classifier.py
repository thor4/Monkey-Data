#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 02:58:53 2017

@author: bconkli4
"""

from sklearn.ensemble import RandomForestClassifier
import time

forest_clf = RandomForestClassifier(random_state=42)

start = time.time()
forest_clf.fit(X_train, y_train)
print("--- %s seconds ---" % (time.time() - start))

#obtain accuracy
start=time.time()
accuracyrf = cross_val_score(forest_clf, X_train, y_train, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))