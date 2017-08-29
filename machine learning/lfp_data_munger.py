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

#load filtered 0-20 Hz phase data
c_i = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/phase-fft-0_20_filtered-base_norm_subsample.mat')
correct = c_i['cFiltPhaseSub']
incorrect = c_i['iFiltPhaseSub']

#load filtered 0-20 Hz fft data
c_i = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/spec-fft-0_20_filtered.mat')
correct = c_i['fftCfilt']
incorrect = c_i['fftIfilt']

#load filtered 0-20 Hz fft data subsampled
c_i = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/spec-fft-0_20_filtered-base_norm_subsample.mat')
correct = c_i['cFiltPowerSub']
incorrect = c_i['iFiltPowerSub']

#load filtered 0-20 Hz fft data by region
c_i = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/spec-fft-0_20_filtered_region.mat')
correctF = c_i['cfNorm']
correctP = c_i['cpNorm']
incorrectF = c_i['ifNorm']
incorrectP = c_i['ipNorm']

#load raw data subsamples
raw = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/raw-0_20_filtered-base_norm_subsample.mat')
correct = raw['cSub']
incorrect = raw['iSub']

#load raw data epochs
raw = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/raw-epochs-0_20_filtered-base_norm_subsample.mat')
correct1 = raw['cSub1']
correct2 = raw['cSub2']
correct3= raw['cSub3']
incorrect1 = raw['iSub1']
incorrect2 = raw['iSub2']
incorrect3= raw['iSub3']

#load raw data stats epochs
raw = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/raw-epochs-stats-0_20_filtered-base_norm_subsample.mat')
c1stat = raw['c1stat']
c2stat = raw['c2stat']
c3stat = raw['c3stat']
i1stat = raw['i1stat']
i2stat = raw['i2stat']
i3stat = raw['i3stat']

#load raw data by broadmann area
c_i = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/raw-0_20_filtered_base_norm-electode_sites.mat')
correctdpfc = np.transpose(c_i['CNormFiltdpfc'])
correctlip = np.transpose(c_i['CNormFiltlip'])
incorrectdpfc = np.transpose(c_i['INormFiltdpfc'])
incorrectlip = np.transpose(c_i['INormFiltlip'])

#prep raw data subsample
#labels
c = np.ones(len(correct))
i = np.zeros(len(incorrect))


#prep fft power data by region
#labels
cF = np.ones(len(correctF))
cP = np.ones(len(correctP))
iF = np.zeros(len(incorrectF))
iP = np.zeros(len(incorrectP))

#prep filtered norm raw data by area
#labels
cdpfc = np.ones(len(correctdpfc))
clip = np.ones(len(correctlip))
idpfc = np.zeros(len(incorrectdpfc))
ilip = np.zeros(len(incorrectlip))

#prep raw data
#labels
c1 = np.ones(len(correct1))
c2 = np.ones(len(correct2))
c3 = np.ones(len(correct3))
i1 = np.zeros(len(incorrect1))
i2 = np.zeros(len(incorrect2))
i3 = np.zeros(len(incorrect3))

#prep raw data stats
c1 = np.ones(len(c1stat))
c2 = np.ones(len(c2stat))
c3 = np.ones(len(c3stat))
i1 = np.zeros(len(i1stat))
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

#combine data and labels fft power by region
f_data = np.vstack((correctF,incorrectF))
f_response = np.vstack((cF[:,None],iF[:,None]))
f_response = f_response.reshape(len(f_data),)
p_data = np.vstack((correctP,incorrectP))
p_response = np.vstack((cP[:,None],iP[:,None]))
p_response = p_response.reshape(len(p_data),)

#combine data and labels raw filtered norm by area
dpfc_data = np.vstack((correctdpfc,incorrectdpfc))
dpfc_response = np.vstack((cdpfc[:,None],idpfc[:,None]))
dpfc_response = dpfc_response.reshape(len(dpfc_data),)
lip_data = np.vstack((correctlip,incorrectlip))
lip_response = np.vstack((clip[:,None],ilip[:,None]))
lip_response = lip_response.reshape(len(lip_data),)

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
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-epochs-0_20_filtered-base_norm_subsample.h5', 'w')
h5f.create_dataset('early-d', data=raw_early_data)
h5f.create_dataset('early-r', data=early_response)
h5f.create_dataset('mid-d', data=raw_mid_data)
h5f.create_dataset('mid-r', data=mid_response)
h5f.create_dataset('late-d', data=raw_late_data)
h5f.create_dataset('late-r', data=late_response)
h5f.close()

###for raw epoch data stats subsampled
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-epochs-stats-0_20_filtered-base_norm_subsample.h5', 'w')
h5f.create_dataset('early-d', data=raw_early_data)
h5f.create_dataset('early-r', data=early_response)
h5f.create_dataset('mid-d', data=raw_mid_data)
h5f.create_dataset('mid-r', data=mid_response)
h5f.create_dataset('late-d', data=raw_late_data)
h5f.create_dataset('late-r', data=late_response)
h5f.close()

###for raw subsampled data
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-0_20_filtered-base_norm_subsample.mat.h5', 'w')
h5f.create_dataset('data', data=data)
h5f.create_dataset('response', data=response)
h5f.close()

###for filtered 0-20Hz fft data
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-spec-fft-0_20_filtered.h5', 'w')
h5f.create_dataset('data', data=data)
h5f.create_dataset('response', data=response)
h5f.close()

###for filtered 0-20Hz phase data
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-phase-fft-0_20_filtered-base_norm_subsample.h5', 'w')
h5f.create_dataset('data', data=data)
h5f.create_dataset('response', data=response)
h5f.close()

###for filtered 0-20Hz fft data subsampled
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-spec-fft-0_20_filtered_subsample.h5', 'w')
h5f.create_dataset('data', data=data)
h5f.create_dataset('response', data=response)
h5f.close()

###for filtered 0-20Hz fft data by region
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-spec-fft-0_20_filtered_region.h5', 'w')
h5f.create_dataset('dataF', data=f_data)
h5f.create_dataset('responseF', data=f_response)
h5f.create_dataset('dataP', data=p_data)
h5f.create_dataset('responseP', data=p_response)
h5f.close()

###for filtered 0-20Hz raw data by area
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-0_20_filtered_base_norm-area-dpfc_lip.h5', 'w')
h5f.create_dataset('datadpfc', data=dpfc_data)
h5f.create_dataset('responsedpfc', data=dpfc_response)
h5f.create_dataset('datalip', data=lip_data)
h5f.create_dataset('responselip', data=lip_response)
h5f.close()

#load h5py
h5f = h5py.File('lfp_data','r')
lfp_data = h5f['lfp_data'][:]
h5f.close()