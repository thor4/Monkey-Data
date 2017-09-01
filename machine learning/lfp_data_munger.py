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
c_i = scipy.io.loadmat('/mnt/ceph/home/bconkli4/Documents/data/raw-0_20_filtered_base_norm-9_areas.mat')
correct6dr = np.transpose(c_i['CNormFilt6dr'])
correct8ad = np.transpose(c_i['CNormFilt8ad'])
correct8b = np.transpose(c_i['CNormFilt8b'])
correct9l = np.transpose(c_i['CNormFilt9l'])
correctpe = np.transpose(c_i['CNormFiltpe'])
correctpec = np.transpose(c_i['CNormFiltpec'])
correctpg = np.transpose(c_i['CNormFiltpg'])
incorrect6dr = np.transpose(c_i['INormFilt6dr'])
incorrect8ad = np.transpose(c_i['INormFilt8ad'])
incorrect8b = np.transpose(c_i['INormFilt8b'])
incorrect9l = np.transpose(c_i['INormFilt9l'])
incorrectpe = np.transpose(c_i['INormFiltpe'])
incorrectpec = np.transpose(c_i['INormFiltpec'])
incorrectpg = np.transpose(c_i['INormFiltpg'])
#sample from correct to even-out the correct/incorrect ration
#32,161 for dpfc, 3,958 for lip
correct6dr = correct6dr[np.random.choice(correct6dr.shape[0], len(incorrect6dr), replace=False),:];
correct8ad = correct8ad[np.random.choice(correct8ad.shape[0], len(incorrect8ad), replace=False),:];
correct8b = correct8b[np.random.choice(correct8b.shape[0], len(incorrect8b), replace=False),:];
incorrect9l = incorrect9l[np.random.choice(incorrect9l.shape[0], len(correct9l), replace=False),:];
correctpe = correctpe[np.random.choice(correctpe.shape[0], len(incorrectpe), replace=False),:];
correctpec = correctpec[np.random.choice(correctpec.shape[0], len(incorrectpec), replace=False),:];
correctpg = correctpg[np.random.choice(correctpg.shape[0], len(incorrectpg), replace=False),:];

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
c6dr = np.ones(len(correct6dr))
c8ad = np.ones(len(correct8ad))
c8b = np.ones(len(correct8b))
c9l = np.ones(len(correct9l))
cpe = np.ones(len(correctpe))
cpec = np.ones(len(correctpec))
cpg = np.ones(len(correctpg))
i6dr = np.zeros(len(incorrect6dr))
i8ad = np.zeros(len(incorrect8ad))
i8b = np.zeros(len(incorrect8b))
i9l = np.zeros(len(incorrect9l))
ipe = np.zeros(len(incorrectpe))
ipec = np.zeros(len(incorrectpec))
ipg = np.zeros(len(incorrectpg))

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
data_6dr = np.vstack((correct6dr,incorrect6dr))
response_6dr = np.vstack((c6dr[:,None],i6dr[:,None]))
response_6dr = response_6dr.reshape(len(data_6dr),)
data_8ad = np.vstack((correct8ad,incorrect8ad))
response_8ad = np.vstack((c8ad[:,None],i8ad[:,None]))
response_8ad = response_8ad.reshape(len(data_8ad),)
data_8b = np.vstack((correct8b,incorrect8b))
response_8b = np.vstack((c8b[:,None],i8b[:,None]))
response_8b = response_8b.reshape(len(data_8b),)
data_9l = np.vstack((correct9l,incorrect9l))
response_9l = np.vstack((c9l[:,None],i9l[:,None]))
response_9l = response_9l.reshape(len(data_9l),)
data_pe = np.vstack((correctpe,incorrectpe))
response_pe = np.vstack((cpe[:,None],ipe[:,None]))
response_pe = response_pe.reshape(len(data_pe),)
data_pec = np.vstack((correctpec,incorrectpec))
response_pec = np.vstack((cpec[:,None],ipec[:,None]))
response_pec = response_pec.reshape(len(data_pec),)
data_pg = np.vstack((correctpg,incorrectpg))
response_pg = np.vstack((cpg[:,None],ipg[:,None]))
response_pg = response_pg.reshape(len(data_pg),)

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

###for filtered 0-20Hz raw data by area even
h5f = h5py.File('/mnt/ceph/home/bconkli4/Documents/data/ml/input-raw-0_20_filtered_base_norm-area-6dr_8ad_8b_9l_pe_pec_pg-even.h5', 'w')
h5f.create_dataset('data6dr', data=data_6dr)
h5f.create_dataset('response6dr', data=response_6dr)
h5f.create_dataset('data8ad', data=data_8ad)
h5f.create_dataset('response8ad', data=response_8ad)
h5f.create_dataset('data8b', data=data_8b)
h5f.create_dataset('response8b', data=response_8b)
h5f.create_dataset('data9l', data=data_9l)
h5f.create_dataset('response9l', data=response_9l)
h5f.create_dataset('datape', data=data_pe)
h5f.create_dataset('responsepe', data=response_pe)
h5f.create_dataset('datapec', data=data_pec)
h5f.create_dataset('responsepec', data=response_pec)
h5f.create_dataset('datapg', data=data_pg)
h5f.create_dataset('responsepg', data=response_pg)
h5f.close()

#load h5py
h5f = h5py.File('lfp_data','r')
lfp_data = h5f['lfp_data'][:]
h5f.close()