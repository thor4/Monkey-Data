#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 17:01:03 2017

@author: bconkli4
"""

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