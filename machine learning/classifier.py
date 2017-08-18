# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 17:00:37 2017

@author: bryan
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import SGDClassifier
from sklearn.base import BaseEstimator
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, roc_curve

class CorrectClassifier(BaseEstimator):
    def fit(self, X, y=None):
        pass
    def predict(self, X):
        return np.ones((len(X), 1), dtype=bool)

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

load h5py
h5f = h5py.File('/home/bconkli4/Documents/Python/ml/datasets/response.h5','r')
X, y = h5f['data'][:], h5f['response'][:]
h5f.close()

#split 80/20 train/test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
#shuffle to ensure all cv folds will be similar. also some learning algorithms 
#are sensitive to order of training instances. they perform poorly if they get 
#many similar instances in a row
shuffle_index = np.random.permutation(392036)
X_train, y_train = X_train[shuffle_index], y_train[shuffle_index]
y_train = y_train.astype(bool)
#X_train = joblib.load('X_train.pkl')
#X_test = np.load('X_test.npy')
#y_train = np.load('y_train.npy')
#y_test = np.load('y_test.npy')


#classify using SGD
sgd_clf = SGDClassifier(random_state=42)
sgd_clf.fit(X_train, y_train)

#performance measures of classifiers
#accuracy
#cross_val_score(sgd_clf, X_train, y_train, cv=5, scoring="accuracy")
##check what would happen if prediction was always "correct": right 72% of the time- this is better than the SGD classifier accuracy
#correct_clf = CorrectClassifier()
#cross_val_score(correct_clf, X_train, y_train, cv=5, 
#                scoring="accuracy")
#this is why we don't use accuracy when dealing with skewed datasets
#get precision, recall, f1 scores of SGD clf
y_train_pred = cross_val_predict(sgd_clf, X_train, y_train, cv=5)
confusion_matrix(y_train, y_train_pred)
sgd_clf_precision = precision_score(y_train, y_train_pred)
sgd_clf_recall = recall_score(y_train, y_train_pred)
sgd_clf_f1 = f1_score(y_train, y_train_pred)
#plot precision_recall curve as function of threshold
y_scores = cross_val_predict(sgd_clf, X_train, y_train, cv=5, 
                             method="decision_function")
precisions, recalls, thresholds = precision_recall_curve(y_train, 
                                                         y_scores)
plot_precision_recall_vs_threshold(precisions, recalls, thresholds)
plt.show()

plt.figure()
plot_precision_vs_recall(precisions, recalls)
#these plots are showing that this SGD classifier doesn't work
#let's check on the ROC curve
fpr, tpr, thresholds = roc_curve(y_train, y_scores)
plt.figure()
plot_roc_curve(fpr, tpr)
plt.show()
sgd_roc_auc_score = roc_auc_score(y_train, y_scores)
#time for a random forest clf
forest_clf = RandomForestClassifier(random_state=42)
y_probas_forest = cross_val_predict(forest_clf, X_train, y_train, cv=5, 
                                    method="predict_proba")
joblib.dump(y_probas_forest, 'y_probas_forest.pkl') #save RF class probs
#y_probas_forest = joblib.load('y_probas_forest.pkl') #load RF class probs
y_train_pred_forest = cross_val_predict(forest_clf, X_train, y_train, 
                                        cv=5)
joblib.dump(y_train_pred_forest, 'y_train_pred_forest.pkl') #save RF predictions
#y_train_pred_forest = joblib.load('y_train_pred_forest.pkl') #load RF class probs
forest_clf_precision = precision_score(y_train, y_train_pred_forest)
forest_clf_recall = recall_score(y_train, y_train_pred_forest)
forest_clf_f1 = f1_score(y_train, y_train_pred_forest)
y_scores_forest = y_probas_forest[:,1]
fpr_forest, tpr_forest, thresholds_forest = roc_curve(y_train, y_scores_forest)
plt.figure()
plt.plot(fpr, tpr, "b:", label="SGD")
plot_roc_curve(fpr_forest, tpr_forest, "Random Forest")
plt.legend(loc="lower right")
plt.show()
forest_roc_auc_score = roc_auc_score(y_train, y_scores_forest)