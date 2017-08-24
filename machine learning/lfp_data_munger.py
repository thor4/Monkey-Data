# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 11:44:37 2017

@author: bryan
"""

import scipy.io
import h5py
import numpy as np
#load data- on koko or local and whether to load data by location (True) in
#addition to response, or just response (False)
workspace = 'koko'
location = False

local_root = "E:/OneDrive/Documents/PhD @ FAU/research/High Frequency FP Activity in VWM/data/"
koko_root = "/home/bconkli4/Documents/MATLAB/data/"
correct_path = "Correct-areas"
incorrect_path = "Incorrect-areas"

if (workspace == 'koko') & (location):
    correct_f_koko = koko_root + correct_path + "/correctFrontal.mat"
    incorrect_f_koko = koko_root + incorrect_path + "/incorrectFrontal.mat"
    correct_p_koko = koko_root + correct_path + "/correctParietal.mat"
    incorrect_p_koko = koko_root + incorrect_path + "/incorrectParietal.mat"
    correctFrontal = scipy.io.loadmat(correct_f_koko)
    correctFrontal = correctFrontal['correctFrontal'].T
    incorrectFrontal = scipy.io.loadmat(incorrect_f_koko)
    incorrectFrontal = incorrectFrontal['incorrectFrontal'].T
    correctParietal = scipy.io.loadmat(correct_p_koko)
    correctParietal = correctParietal['correctParietal'].T
    incorrectParietal = scipy.io.loadmat(incorrect_p_koko)
    incorrectParietal = incorrectParietal['incorrectParietal'].T
elif (workspace == 'koko'):
    correct_koko = koko_root + "correct.mat"
    incorrect_koko = koko_root + "incorrect.mat"
    correct = h5py.File(correct_koko, 'r') #too big for scipy, need h5py
    correct = correct['correct'].value
    incorrect = scipy.io.loadmat(incorrect_koko)
    incorrect = incorrect['incorrect'].T
    del correct_koko; del correct_path; del incorrect_koko; del incorrect_path;
    del koko_root; del local_root; del location; del workspace;
elif (location):
    correct_f_local = local_root + correct_path + "/correctFrontal.mat"
    incorrect_f_local = local_root + incorrect_path + "/incorrectFrontal.mat"
    correct_p_local = local_root + correct_path + "/correctParietal.mat"
    incorrect_p_local = local_root + incorrect_path + "/incorrectParietal.mat"
    correctFrontal = scipy.io.loadmat(correct_f_local)
    correctFrontal = correctFrontal['correctFrontal'].T
    incorrectFrontal = scipy.io.loadmat(incorrect_f_local)
    incorrectFrontal = incorrectFrontal['incorrectFrontal'].T
    correctParietal = scipy.io.loadmat(correct_p_local)
    correctParietal = correctParietal['correctParietal'].T
    incorrectParietal = scipy.io.loadmat(incorrect_p_local)
    incorrectParietal = incorrectParietal['incorrectParietal'].T
else:
    correct_local = local_root + "correct.mat"
    incorrect_local = local_root + "incorrect.mat"
    correct = h5py.File(correct_local, 'r')
    correct = correct['correct'].value
    incorrect = scipy.io.loadmat(incorrect_local)
    incorrect = incorrect['incorrect'].T
    del correct_local; del correct_path; del incorrect_local; 
    del incorrect_path; del local_root; del location; del workspace;

#%reset clears variable workspace
#prep chronux spectral data
c_i = scipy.io.loadmat('/home/bconkli4/Documents/MATLAB/data/spec-500-5-c_i.mat')
c_i = scipy.io.loadmat('spec-500-5-c_i.mat')
correct = c_i['Pc']
# load the AR spectral power data
correct = scipy.io.loadmat('correct_ARpower.mat')
correct = correct['arpower_c']
incorrect = scipy.io.loadmat('incorrect_ARpower.mat')
incorrect = incorrect['arpower_i']
incorrect = c_i['Pi']
time = c_i['Tc']
freq = c_i['Fc']
cor = np.transpose(correct)

#load raw data epochs
raw = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/raw-epochs.mat')
correct1 = np.transpose(raw['correct1'])
correct2 = np.transpose(raw['correct2'])
correct3= np.transpose(raw['correct3'])
incorrect1 = np.transpose(raw['incorrect1'])
incorrect2 = np.transpose(raw['incorrect2'])
incorrect3= np.transpose(raw['incorrect3'])

#prep raw data
#labels
c1 = np.ones(len(correct1))
c2 = np.ones(len(correct2))
c3 = np.ones(len(correct3))
i1 = np.zeros(len(incorrect1))
i2 = np.zeros(len(incorrect2))
i3 = np.zeros(len(incorrect3))

#combine data and labels
spec500_data = np.vstack((correct,incorrect))
response = np.vstack((c[:,None],i[:,None]))
response = response.reshape(len(spec500_data),)
##for AR spectral power data
arpower_data = np.vstack((correct,incorrect))
response = np.vstack((c[:,None],i[:,None]))
response = response.reshape(len(arpower_data),)
###for raw epoch data
raw_early_data = np.vstack((correct1,incorrect1))
early_response = np.vstack((c1[:,None],i1[:,None]))
early_response = early_response.reshape(len(raw_early_data),)
raw_mid_data = np.vstack((correct2,incorrect2))
mid_response = np.vstack((c2[:,None],i2[:,None]))
mid_response = mid_response.reshape(len(raw_mid_data),)
raw_late_data = np.vstack((correct3,incorrect3))
late_response = np.vstack((c3[:,None],i3[:,None]))
late_response = late_response.reshape(len(raw_late_data),)

#save h5py
h5f = h5py.File('spec500-resp.h5', 'w')
h5f.create_dataset('data', data=spec500_data)
h5f.create_dataset('response', data=response)
h5f.close()
##for AR spectral power data
h5f = h5py.File('arpower-resp.h5', 'w')
h5f.create_dataset('data', data=arpower_data)
h5f.create_dataset('response', data=response)
h5f.close()

###for raw epoch data
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/raw-epochs.h5', 'w')
h5f.create_dataset('early-data', data=raw_early_data)
h5f.create_dataset('early-response', data=early_response)
h5f.create_dataset('mid-data', data=raw_mid_data)
h5f.create_dataset('mid-response', data=mid_response)
h5f.create_dataset('late-data', data=raw_late_data)
h5f.create_dataset('late-response', data=late_response)
h5f.close()

#load h5py
h5f = h5py.File('lfp_data','r')
lfp_data = h5f['lfp_data'][:]
h5f.close()