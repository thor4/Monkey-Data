# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 10:43:40 2018

@author: bryan
"""

import h5py
import numpy as np

f = h5py.File('mGoodStableRule1PingRej-split_by_Day_BehResp_and_Chan.mat','r')

g = f['monkey/day']

data = f.get('monkey/day')

type(g) #dataset

data.keys()

data.shape
data.size
data.dtype #object

data[:]

monkey1 = data.get[1,0]


