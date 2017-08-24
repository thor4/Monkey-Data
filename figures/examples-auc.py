#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 23:18:16 2017

@author: bconkli4
"""
#example 1
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt

plt.figure(0).clf()

pred = np.random.rand(1000)
label = np.random.randint(2, size=1000)
fpr, tpr, thresh = metrics.roc_curve(label, pred)
auc = metrics.roc_auc_score(label, pred)
plt.plot(fpr,tpr,label="data 1, auc="+str(auc))

pred = np.random.rand(1000)
label = np.random.randint(2, size=1000)
fpr, tpr, thresh = metrics.roc_curve(label, pred)
auc = metrics.roc_auc_score(label, pred)
plt.plot(fpr,tpr,label="data 2, auc="+str(auc))

plt.legend(loc=0)

#example 2
y1=[1, 0, 0, 0, 1];
y2=[1, 1, 0, 0, 1];
y3=[0, 0, 1, 0, 1];
from sklearn.metrics import *
import pandas as pd
import matplotlib.pyplot as plt

plt.figure(1)
fpr, tpr, thresholds = roc_curve(y1, y2)
auc1 = auc(fpr,tpr)

plt.plot(fpr, tpr,label="AUC Y2:{0}".format(auc1),color='red', linewidth=2)

fpr, tpr, thresholds = roc_curve(y1, y3)
auc1 = auc(fpr,tpr)

plt.plot(fpr, tpr,label="AUC Y3:{0}".format(auc1),color='blue', linewidth=2)

plt.plot([0, 1], [0, 1], 'k--', lw=1) 
plt.xlim([0.0, 1.0]) 
plt.ylim([0.0, 1.05])

plt.xlabel('False Positive Rate')  
plt.ylabel('True Positive Rate') 
plt.title('ROC') 
plt.grid(True)
plt.legend(loc="lower right")
plt.show()