#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 02:58:53 2017

@author: bconkli4
"""

from scipy.stats import randint as sp_randint
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
from sklearn.model_selection import RandomizedSearchCV
import time

forest_clf = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf.fit(X_trainfs, y_train)
#fit standardized
forest_clf = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf.fit(stdx, y_train)
#fit standardized fft power by region
start = time.time()
forest_clf1 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf1.fit(X_trainfstd, y_trainf)
forest_clf2 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf2.fit(X_trainpstd, y_trainp)
print("--- %s seconds ---" % (time.time() - start))
#fit filtered raw data by area
start = time.time()
forest_clf1 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf1.fit(X_train6drfs, y_train6dr)
forest_clf2 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf2.fit(X_train8adfs, y_train8ad)
forest_clf3 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf3.fit(X_train8bfs, y_train8b)
forest_clf4 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf4.fit(X_train9lfs, y_train9l)
forest_clf5 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf5.fit(X_traindpfcfs, y_traindpfc)
forest_clf6 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf6.fit(X_trainlipfs, y_trainlip)
forest_clf7 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf7.fit(X_trainpecfs, y_trainpec)
forest_clf8 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf8.fit(X_trainpefs, y_trainpe)
forest_clf9 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf9.fit(X_trainpgfs, y_trainpg)
print("--- %s seconds ---" % (time.time() - start))
#epochs stats
forest_clf1 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf2 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf3 = RandomForestClassifier(n_estimators=20, random_state=42)

#filtered 0-20Hz fft fit
start = time.time()
forest_clf1 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf1.fit(X_trainfs, y_train)
forest_clf2 = RandomForestClassifier(n_estimators=20, random_state=42)
forest_clf2.fit(stdx, y_train)
print("--- %s seconds ---" % (time.time() - start))

#raw epochs stats
start = time.time()
forest_clf1.fit(X_train1, y_train1)
forest_clf2.fit(X_train2, y_train2)
forest_clf3.fit(X_train3, y_train3)
print("--- %s seconds ---" % (time.time() - start))

#obtain accuracy
start=time.time()
accuracyRFfs = cross_val_score(forest_clf, X_trainfs, y_train, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))

#accuracy for standadized data
start=time.time()
accuracyRFstd= cross_val_score(forest_clf, stdx, y_train, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))

#obtain accuracy filtered fft by region
start=time.time()
accuracyRFfstd = cross_val_score(forest_clf1, X_trainfstd, y_trainf, cv=5, scoring="accuracy")
accuracyRFpstd = cross_val_score(forest_clf2, X_trainpstd, y_trainp, cv=5, scoring="accuracy")
print("--- %s seconds ---" % (time.time() - start))

#obtain accuracy raw data by area
start=time.time()
accuracyRF6drfs = cross_val_score(forest_clf1, X_train6drfs, y_train6dr, cv=5, scoring="accuracy")
accuracyRF8adfs = cross_val_score(forest_clf2, X_train8adfs, y_train8ad, cv=5, scoring="accuracy")
accuracyRF8bfs = cross_val_score(forest_clf3, X_train8bfs, y_train8b, cv=5, scoring="accuracy")
accuracyRF9lfs = cross_val_score(forest_clf4, X_train9lfs, y_train9l, cv=5, scoring="accuracy")
accuracyRFdpfcfs = cross_val_score(forest_clf1, X_traindpfcfs, y_traindpfc, cv=5, scoring="accuracy")
accuracyRFlipfs = cross_val_score(forest_clf2, X_trainlipfs, y_trainlip, cv=5, scoring="accuracy")
accuracyRFpefs = cross_val_score(forest_clf3, X_trainpefs, y_trainpe, cv=5, scoring="accuracy")
accuracyRFpecfs = cross_val_score(forest_clf4, X_trainpecfs, y_trainpec, cv=5, scoring="accuracy")
accuracyRFpgfs = cross_val_score(forest_clf5, X_trainpgfs, y_trainpg, cv=5, scoring="accuracy")
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
              "max_features": sp_randint(1, 11),
              "min_samples_split": sp_randint(2, 11),
              "min_samples_leaf": sp_randint(1, 11),
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}
# run randomized search
n_iter_search = 1000
random_search = RandomizedSearchCV(forest_clf8, param_distributions=param_dist,
                                   n_iter=n_iter_search)

start = time.time()
random_search.fit(X_trainpefs, y_trainpe)
print("RandomizedSearchCV took %.2f seconds for %d candidates"
      " parameter settings." % ((time.time() - start), n_iter_search))
accuracyRFpefsTuned = random_search.score(X_trainpefs, y_trainpe)
accuracyRFpefsTest = random_search.score(X_testpefs, y_testpe)
print("[INFO] randomized search accuracy: {:.2f}%".format(accuracyRFpefsTest * 100))
print("[INFO] randomized search best parameters: {}".format(
	random_search.best_params_))
#random_search is 6dr name

#accuracy
report(random_search.cv_results_)
random_search.best_estimator_

#save model
joblib.dump(random_search, '/mnt/ceph/home/bconkli4/Documents/data/ml/models/peRF.pkl') 

#load model
clf = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/lipRF.pkl') 