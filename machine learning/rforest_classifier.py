#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 02:58:53 2017

@author: bconkli4
"""

import pandas as pd
from scipy.stats import randint as sp_randint
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
import time

forest_clf = RandomForestClassifier(n_estimators=20, random_state=42)
#epochs stats
start = time.time()
forest_clf1 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf1.fit(X_train1, y_train1)
forest_clf2 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf2.fit(X_train2, y_train2)
forest_clf3 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf3.fit(X_train3, y_train3)
print("--- %s seconds ---" % (time.time() - start))

start = time.time()
forest_clf.fit(X_train, y_train)
print("--- %s seconds ---" % (time.time() - start))

#obtain accuracy
start=time.time()
accuracyrf20 = cross_val_score(forest_clf, X_train, y_train, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))

#accuracy epochs stats
start=time.time()
accuracy1rf = cross_val_score(forest_clf1, X_train1, y_train1, cv=5, scoring="accuracy")
accuracy2rf = cross_val_score(forest_clf2, X_train2, y_train2, cv=5, scoring="accuracy")
accuracy3rf = cross_val_score(forest_clf3, X_train3, y_train3, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))

#hyperparameter tuning
# Utility function to report best scores
def report(results, n_top=3):
    for i in range(1, n_top + 1):
        candidates = np.flatnonzero(results['rank_test_score'] == i)
        for candidate in candidates:
            print("Model with rank: {0}".format(i))
            print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
                  results['mean_test_score'][candidate],
                  results['std_test_score'][candidate]))
            print("Parameters: {0}".format(results['params'][candidate]))
            print("")
# specify parameters and distributions to sample from
param_dist = {"max_depth": [3, None],
              "max_features": sp_randint(1, 8),
              "min_samples_split": sp_randint(2, 11),
              "min_samples_leaf": sp_randint(1, 11),
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}
# run randomized search
n_iter_search = 1000
random_search1 = RandomizedSearchCV(forest_clf1, param_distributions=param_dist,
                                   n_iter=n_iter_search)
random_search2 = RandomizedSearchCV(forest_clf2, param_distributions=param_dist,
                                   n_iter=n_iter_search)
random_search3 = RandomizedSearchCV(forest_clf3, param_distributions=param_dist,
                                   n_iter=n_iter_search)

start = time.time()
random_search1.fit(X_train1, y_train1)
print("RandomizedSearchCV took %.2f seconds for %d candidates"
      " parameter settings." % ((time.time() - start), n_iter_search))
random_search2.fit(X_train2, y_train2)
print("RandomizedSearchCV took %.2f seconds for %d candidates"
      " parameter settings." % ((time.time() - start), n_iter_search))
random_search3.fit(X_train3, y_train3)
print("RandomizedSearchCV took %.2f seconds for %d candidates"
      " parameter settings." % ((time.time() - start), n_iter_search))

#accuracy
report(random_search1.cv_results_)
random_search1.best_estimator_
report(random_search2.cv_results_)
random_search2.best_estimator_
report(random_search3.cv_results_)
random_search3.best_estimator_

late_results = pd.DataFrame(random_search3.cv_results_)

