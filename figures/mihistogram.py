#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 03:20:07 2017

@author: bconkli4
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


scoreslipF = np.load('/mnt/ceph/home/bconkli4/Documents/data/ml/temp/scoreslipF.npy')
scorespeF = np.load('/mnt/ceph/home/bconkli4/Documents/data/ml/temp/scorespeF.npy')

y_scoresLR6dr = np.load('/mnt/ceph/home/bconkli4/Documents/data/ml/temp/y_scoresLR6dr.npy')
y_scoresNN6dr = np.load('/mnt/ceph/home/bconkli4/Documents/data/ml/temp/y_scoresNN6dr.npy')
y_scoresRF6dr = np.load('/mnt/ceph/home/bconkli4/Documents/data/ml/temp/y_scoresRF6dr.npy')
y_test6dr = np.load('/mnt/ceph/home/bconkli4/Documents/data/ml/temp/y_test6dr.npy')



plt.figure()
#6DR
plt.subplot(1, 3, 1)
df6dr = pd.DataFrame(scores6dr);
top2506dr = df6dr.nlargest(250,0);
topidx6dr = np.sort(top2506dr.index.values);
plt.hist(topidx6dr, bins=[0,100,200,300,400,500,600,700,810], facecolor="#8ED081", alpha=0.75, edgecolor="black", linewidth=1.0)
plt.title('6DR', fontsize=20)
plt.xlabel('Time (ms)', fontsize=16)
plt.grid(False)
plt.axis([0, 810, 0, 120], 'tight')
plt.tick_params(labelsize=16)
#LIP
plt.subplot(1, 3, 2)
dflip = pd.DataFrame(scoreslip);
top250lip = dflip.nlargest(250,0);
topidxlip = np.sort(top250lip.index.values);
plt.hist(topidxlip, bins=[0,100,200,300,400,500,600,700,810], facecolor="#8ED081", alpha=0.75, edgecolor="black", linewidth=1.0)
plt.title('LIP', fontsize=20)
plt.xlabel('Time (ms)', fontsize=16)
plt.grid(False)
plt.axis([0, 810, 0, 120], 'tight')
plt.tick_params(labelsize=16)
#PE
plt.subplot(1, 3, 3)
dfpe = pd.DataFrame(scorespe);
top250pe = dfpe.nlargest(250,0);
topidxpe = np.sort(top250pe.index.values);
plt.hist(topidxpe, bins=[0,100,200,300,400,500,600,700,810], facecolor="#8ED081", alpha=0.75, edgecolor="black", linewidth=1.0)
plt.title('PE', fontsize=20)
plt.xlabel('Time (ms)', fontsize=16)
plt.grid(False)
plt.axis([0, 810, 0, 120], 'tight')


plt.tick_params(labelsize=16)