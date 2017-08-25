# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 11:44:37 2017

@author: bryan
"""

import scipy.io
import h5py
import numpy as np
import pandas as pd

#%reset clears variable workspace
#prep chronux spectral data
c_i = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/spec-chronux-spectrum-base-norm.mat')
correct = np.transpose(c_i['ScNorm'])
incorrect = np.transpose(c_i['SiNorm'])
# load the AR spectral power data
correct = scipy.io.loadmat('correct_ARpower.mat')
correct = correct['arpower_c']
incorrect = scipy.io.loadmat('incorrect_ARpower.mat')
incorrect = incorrect['arpower_i']
incorrect = c_i['Pi']
time = c_i['Tc']
freq = c_i['Fc']
cor = np.transpose(correct)

#load raw data subsamples
raw = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/raw-0_20_filtered-base_norm_subsample.mat')
correct = raw['cSub']
incorrect = raw['iSub']

#load raw data epochs
raw = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/raw-epochs.mat')
correct1 = np.transpose(raw['correct1'])
correct2 = np.transpose(raw['correct2'])
correct3= np.transpose(raw['correct3'])
incorrect1 = np.transpose(raw['incorrect1'])
incorrect2 = np.transpose(raw['incorrect2'])
incorrect3= np.transpose(raw['incorrect3'])

#load raw data stats epochs
raw = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/raw-epochs-stats.mat')
c1stat = raw['c1stat']
c2stat = raw['c2stat']
c3stat = raw['c3stat']
i1stat = raw['i1stat']
i2stat = raw['i2stat']
i3stat = raw['i3stat']

#prep raw data subsample
#labels
c = np.ones(len(correct))
i = np.zeros(len(incorrect))

#prep raw data
#labels
c = np.ones(len(correct))
c2 = np.ones(len(c2stat))
c3 = np.ones(len(c3stat))
i = np.zeros(len(incorrect))
i2 = np.zeros(len(i2stat))
i3 = np.zeros(len(i3stat))

#combine data and labels
chronux_data = np.vstack((correct,incorrect))
response = np.vstack((c[:,None],i[:,None]))
response = response.reshape(len(chronux_data),)
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

####for raw epoch data stats
raw_early_data = np.vstack((c1stat,i1stat))
early_response = np.vstack((c1[:,None],i1[:,None]))
early_response = early_response.reshape(len(raw_early_data),)
raw_mid_data = np.vstack((c2stat,i2stat))
mid_response = np.vstack((c2[:,None],i2[:,None]))
mid_response = mid_response.reshape(len(raw_mid_data),)
raw_late_data = np.vstack((c3stat,i3stat))
late_response = np.vstack((c3[:,None],i3[:,None]))
late_response = late_response.reshape(len(raw_late_data),)

#####for raw subsampled data
data = np.vstack((correct,incorrect))
response = np.vstack((c[:,None],i[:,None]))
response = response.reshape(len(data),)

#save h5py
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/chronux.h5', 'w')
h5f.create_dataset('data', data=chronux_data)
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

###for raw subsampled data
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-0_20_filtered-base_norm_subsample.mat.h5', 'w')
h5f.create_dataset('data', data=data)
h5f.create_dataset('response', data=response)
h5f.close()

#load h5py
h5f = h5py.File('lfp_data','r')
lfp_data = h5f['lfp_data'][:]
h5f.close()