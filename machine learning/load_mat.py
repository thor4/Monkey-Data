# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 13:57:12 2018

@author: bryan
"""

import scipy.io
import numpy as np


#%reset clears variable workspace
#prep chronux spectral data
m1 = scipy.io.loadmat('m1GoodStableRule1PingRej-split_by_Day_BehResp_and_Chan.mat')
monkey1 = m1['monkey1'] #pull out monkey1 struct
day = monkey1['day'] #pull out day object
days = day[0,0] #pull out day struct
days.dtype #check datatype for keys to index into 
correct_days = days['correct']; incorrect_days = days['incorrect']
correct_days.dtype; del day; del days #cleanup
cor_day1 = correct_days[0,0] #pull out correct day 1
cor_day1.dtype #datatype shows channels recorded from as keys
m1d1CorChan1PEC_ = cor_day1['chan1PEC'] 
m1d1CorChan1PEC = m1d1CorChan1PEC_[0,0]; del m1d1CorChan1PEC_;

m1d1CorChan3dPFC_ = cor_day1['chan3dPFC'] 
m1d1CorChan3dPFC = m1d1CorChan3dPFC_[0,0]; del m1d1CorChan3dPFC_;