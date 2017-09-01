#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 19:06:47 2017

@author: bconkli4
"""

from sklearn.externals import joblib
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

#load model LIP
RF6dr = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/6drRF.pkl') 
NN6dr = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/6drNN.pkl') 
LR6dr = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/6drLR.pkl')
RFlip = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/lipRF.pkl') 
NNlip = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/lipNN.pkl') 
LRlip = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/lipLR.pkl')
RFpe = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/peRF.pkl') 
NNpe = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/peNN.pkl') 
LRpe = joblib.load('/mnt/ceph/home/bconkli4/Documents/data/ml/models/peLR.pkl')

#build ROC curve

#6DR
y_probRF6dr = cross_val_predict(RF6dr, X_test6drfs, y_test6dr, cv=5,
method="predict_proba")
y_scoresRF6dr = y_probRF6dr[:, 1] # score = proba of positive class
fpr_RF6dr, tpr_RF6dr, thresholds_RF6dr = roc_curve(y_test6dr,y_scoresRF6dr)
y_probNN6dr = cross_val_predict(NN6dr, X_test6drfs, y_test6dr, cv=5,
method="predict_proba")
y_scoresNN6dr = y_probNN6dr[:, 1] # score = proba of positive class
fpr_NN6dr, tpr_NN6dr, thresholds_NN6dr = roc_curve(y_test6dr,y_scoresNN6dr)
y_probLR6dr = cross_val_predict(LR6dr, X_test6drfs, y_test6dr, cv=5,
method="predict_proba")
y_scoresLR6dr = y_probLR6dr[:, 1] # score = proba of positive class
fpr_LR6dr, tpr_LR6dr, thresholds_LR6dr = roc_curve(y_test6dr,y_scoresLR6dr)

#LIP
y_probRFlip = cross_val_predict(RFlip, X_testlipfs, y_testlip, cv=5,
method="predict_proba")
y_scoresRFlip = y_probRFlip[:, 1] # score = proba of positive class
fpr_RFlip, tpr_RFlip, thresholds_RFlip = roc_curve(y_testlip,y_scoresRFlip)
y_probNNlip = cross_val_predict(NNlip, X_testlipfs, y_testlip, cv=5,
method="predict_proba")
y_scoresNNlip = y_probNNlip[:, 1] # score = proba of positive class
fpr_NNlip, tpr_NNlip, thresholds_NNlip = roc_curve(y_testlip,y_scoresNNlip)
y_probLRlip = cross_val_predict(LRlip, X_testlipfs, y_testlip, cv=5,
method="predict_proba")
y_scoresLRlip = y_probLRlip[:, 1] # score = proba of positive class
fpr_LRlip, tpr_LRlip, thresholds_LRlip = roc_curve(y_testlip,y_scoresLRlip)

#PE
y_probRFpe = cross_val_predict(RFpe, X_testpefs, y_testpe, cv=5,
method="predict_proba")
y_scoresRFpe = y_probRFpe[:, 1] # score = proba of positive class
fpr_RFpe, tpr_RFpe, thresholds_RFpe = roc_curve(y_testpe,y_scoresRFpe)
y_probNNpe = cross_val_predict(NNpe, X_testpefs, y_testpe, cv=5,
method="predict_proba")
y_scoresNNpe = y_probNNpe[:, 1] # score = proba of positive class
fpr_NNpe, tpr_NNpe, thresholds_NNpe = roc_curve(y_testpe,y_scoresNNpe)
y_probLRpe = cross_val_predict(LRpe, X_testpefs, y_testpe, cv=5,
method="predict_proba")
y_scoresLRpe = y_probLRpe[:, 1] # score = proba of positive class
fpr_LRpe, tpr_LRpe, thresholds_LRpe = roc_curve(y_testpe,y_scoresLRpe)


#plot ROC curve
#6DR
plt.figure()
plt.subplot(1, 3, 1)
plt.title('6DR', fontsize=20)
plt.plot(fpr_RF6dr, tpr_RF6dr, color='#4F81BD', linewidth=4, label="Random Forest")
plt.plot(fpr_NN6dr, tpr_NN6dr, color='#FA7921', linewidth=4, label="Nearest Neighbor")
plt.plot(fpr_LR6dr, tpr_LR6dr, color='#8ED081', linewidth=4, label="Logistic Regression")
plt.axis([0, 1, 0, 1], 'tight')
plt.plot([0, 1], [0, 1], 'k--')
plt.tick_params(labelsize=16)
plt.xlabel('False Positive Rate', fontsize=16)
plt.ylabel('True Positive Rate', fontsize=16)

#LIP
plt.subplot(1, 3, 2)
plt.title('LIP', fontsize=20)
plt.plot(fpr_RFlip, tpr_RFlip, linewidth=4, color='#4F81BD', label="Random Forest")
plt.plot(fpr_NNlip, tpr_NNlip, linewidth=4, color='#FA7921', label="Nearest Neighbor")
plt.plot(fpr_LRlip, tpr_LRlip, linewidth=4, color='#8ED081', label="Logistic Regression")
plt.axis([0, 1, 0, 1], 'tight')
plt.plot([0, 1], [0, 1], 'k--')
plt.tick_params(labelsize=16)
plt.xlabel('False Positive Rate', fontsize=16)
plt.ylabel('True Positive Rate', fontsize=16)

#PE
plt.subplot(1, 3, 3)
plt.title('PE', fontsize=20)
plt.plot(fpr_RFpe, tpr_RFpe, linewidth=4, color='#4F81BD', label="Random Forest")
plt.plot(fpr_NNpe, tpr_NNpe, linewidth=4, color='#FA7921', label="Nearest Neighbor")
plt.plot(fpr_LRpe, tpr_LRpe, linewidth=4, color='#8ED081', label="Logistic Regression")
plt.axis([0, 1, 0, 1], 'tight')
plt.plot([0, 1], [0, 1], 'k--')
plt.tick_params(labelsize=16)
plt.xlabel('False Positive Rate', fontsize=16)
plt.ylabel('True Positive Rate', fontsize=16)

plt.legend(loc="lower right", fontsize=16)
plt.show()

#AUC score
#6DR
aucRF6dr = roc_auc_score(y_test6dr, y_scoresRF6dr);
aucNN6dr = roc_auc_score(y_test6dr, y_scoresNN6dr);
aucLR6dr = roc_auc_score(y_test6dr, y_scoresLR6dr);
#LIP
aucRFlip = roc_auc_score(y_testlip, y_scoresRFlip);
aucNNlip = roc_auc_score(y_testlip, y_scoresNNlip);
aucLRlip = roc_auc_score(y_testlip, y_scoresLRlip);
#PE
aucRFpe = roc_auc_score(y_testpe, y_scoresRFpe);
aucNNpe = roc_auc_score(y_testpe, y_scoresNNpe);
aucLRpe = roc_auc_score(y_testpe, y_scoresLRpe);