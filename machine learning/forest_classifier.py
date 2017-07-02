#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 17:01:03 2017

@author: bconkli4
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, roc_curve

def plot_precision_recall_vs_threshold(precisions, recalls, thresholds):
    plt.plot(thresholds, precisions[:-1], "b--", label="Precision")
    plt.plot(thresholds, recalls[:-1], "g-", label="Recall")
    plt.xlabel("Threshold")
    plt.legend(loc="upper left")
    plt.ylim([0,1])
    
def plot_precision_vs_recall(precisions, recalls):
    plt.plot(recalls, precisions, "b-")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.ylim([0,1])
    plt.xlim([0,1])

def plot_roc_curve(fpr, tpr, label=None):
    plt.plot(fpr, tpr, linewidth=2, label=label)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.axis([0, 1, 0, 1])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

h5f = h5py.File('/home/bconkli4/Documents/Python/ml/datasets/spec500-resp.h5','r')
X, y = h5f['data'][:], h5f['response'][:]
h5f.close()

#split 80/20 train/test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

#time for a random forest clf
forest_clf = RandomForestClassifier(random_state=42)
y_probas_forest = cross_val_predict(forest_clf, X_train, y_train, cv=5, 
                                    method="predict_proba")
y_train_pred_forest = cross_val_predict(forest_clf, X_train, y_train, 
                                        cv=5)
forest_clf_precision = precision_score(y_train, y_train_pred_forest)
forest_clf_recall = recall_score(y_train, y_train_pred_forest)
forest_clf_f1 = f1_score(y_train, y_train_pred_forest)
y_scores_forest = y_probas_forest[:,1]
fpr_forest, tpr_forest, thresholds_forest = roc_curve(y_train, y_scores_forest)
forest_roc_auc_score = roc_auc_score(y_train, y_scores_forest)

#test scores
y_probas_forest = cross_val_predict(forest_clf, X_test, y_test, cv=5, 
                                    method="predict_proba")
y_test_pred_forest = cross_val_predict(forest_clf, X_test, y_test, 
                                        cv=5)
forest_clf_precision = precision_score(y_test, y_test_pred_forest)
forest_clf_recall = recall_score(y_test, y_test_pred_forest)
forest_clf_f1 = f1_score(y_test, y_test_pred_forest)
y_scores_forest = y_probas_forest[:,1]
fpr_forest, tpr_forest, thresholds_forest = roc_curve(y_test, y_scores_forest)
forest_roc_auc_score = roc_auc_score(y_test, y_scores_forest)

#save model
joblib.dump(forest_clf, "model-forest-spec500.pkl")
#forest_clf = joblib.load("model-forest-spec500.pkl") #load model