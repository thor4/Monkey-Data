%load betty
load('betty_6DR_norm_cor_rule1.mat');; load('betty_6DR_norm_inc_rule1.mat');; 
load('betty_8AD_norm_cor_rule1.mat');; load('betty_8AD_norm_inc_rule1.mat');
load('betty_8B_norm_cor_rule1.mat'); load('betty_8B_norm_inc_rule1.mat');
load('betty_dPFC_norm_cor_rule1.mat'); load('betty_dPFC_norm_inc_rule1.mat');
load('betty_LIP_norm_cor_rule1.mat'); load('betty_LIP_norm_inc_rule1.mat');
load('betty_PE_norm_cor_rule1.mat'); load('betty_PE_norm_inc_rule1.mat');
load('betty_PEC_norm_cor_rule1.mat'); load('betty_PEC_norm_inc_rule1.mat');
load('betty_PG_norm_cor_rule1.mat'); load('betty_PG_norm_inc_rule1.mat');
%load clark
load('clark_9L_norm_cor_rule1.mat'); load('clark_9L_norm_inc_rule1.mat');
load('clark_MIP_norm_cor_rule1.mat'); load('clark_MIP_norm_inc_rule1.mat');
load('clark_8B_norm_cor_rule1.mat'); load('clark_8B_norm_inc_rule1.mat');
load('clark_dPFC_norm_cor_rule1.mat'); load('clark_dPFC_norm_inc_rule1.mat');
load('clark_LIP_norm_cor_rule1.mat'); load('clark_LIP_norm_inc_rule1.mat');
load('clark_vPFC_norm_cor_rule1.mat'); load('clark_vPFC_norm_inc_rule1.mat');
load('clark_PEC_norm_cor_rule1.mat'); load('clark_PEC_norm_inc_rule1.mat');
load('clark_PG_norm_cor_rule1.mat'); load('clark_PG_norm_inc_rule1.mat');

%get power using fftPow func with sr=1000, time window is 810ms
%data should be in #trials x nframes
s = 1000; %sample rate
nframes = 810; %time window
%betty
[b6DRpowCorR1, b6DRfreqCorR1] = fftPow(s,b6DRnormCorR1',nframes);
[b6DRpowIncR1, b6DRfreqIncR1] = fftPow(s,b6DRnormIncR1',nframes);
[b8ADpowCorR1, b8ADfreqCorR1] = fftPow(s,b8ADnormCorR1',nframes);
[b8ADpowIncR1, b8ADfreqIncR1] = fftPow(s,b8ADnormIncR1',nframes);
[b8BpowCorR1, b8BfreqCorR1] = fftPow(s,b8BnormCorR1',nframes);
[b8BpowIncR1, b8BfreqIncR1] = fftPow(s,b8BnormIncR1',nframes);
[bdPFCpowCorR1, bdPFCfreqCorR1] = fftPow(s,bdPFCnormCorR1',nframes);
[bdPFCpowIncR1, bdPFCfreqIncR1] = fftPow(s,bdPFCnormIncR1',nframes);
[bLIPpowCorR1, bLIPfreqCorR1] = fftPow(s,bLIPnormCorR1',nframes);
[bLIPpowIncR1, bLIPfreqIncR1] = fftPow(s,bLIPnormIncR1',nframes);
[bPECpowCorR1, bPECfreqCorR1] = fftPow(s,bPECnormCorR1',nframes);
[bPECpowIncR1, bPECfreqIncR1] = fftPow(s,bPECnormIncR1',nframes);
[bPEpowCorR1, bPEfreqCorR1] = fftPow(s,bPEnormCorR1',nframes);
[bPEpowIncR1, bPEfreqIncR1] = fftPow(s,bPEnormIncR1',nframes);
[bPGpowCorR1, bPGfreqCorR1] = fftPow(s,bPGnormCorR1',nframes);
[bPGpowIncR1, bPGfreqIncR1] = fftPow(s,bPGnormIncR1',nframes);
clearvars b6DRnormCorR1 b6DRnormIncR1 b8ADnormCorR1 b8ADnormIncR1 ...
    b8BnormCorR1 b8BnormIncR1 bdPFCnormCorR1 bdPFCnormIncR1 bLIPnormCorR1 ...
    bLIPnormIncR1 bPECnormCorR1 bPECnormIncR1 bPEnormCorR1 bPEnormIncR1 ...
    bPGnormCorR1 bPGnormIncR1

%clark
[c9LpowCorR1, c9LfreqCorR1] = fftPow(s,c9LnormCorR1',nframes);
[c9LpowIncR1, c9LfreqIncR1] = fftPow(s,c9LnormIncR1',nframes);
[cMIPpowCorR1, cMIPfreqCorR1] = fftPow(s,cMIPnormCorR1',nframes);
[cMIPpowIncR1, cMIPfreqIncR1] = fftPow(s,cMIPnormIncR1',nframes);
[c8BpowCorR1, c8BfreqCorR1] = fftPow(s,c8BnormCorR1',nframes);
[c8BpowIncR1, c8BfreqIncR1] = fftPow(s,c8BnormIncR1',nframes);
[cdPFCpowCorR1, cdPFCfreqCorR1] = fftPow(s,cdPFCnormCorR1',nframes);
[cdPFCpowIncR1, cdPFCfreqIncR1] = fftPow(s,cdPFCnormIncR1',nframes);
[cLIPpowCorR1, cLIPfreqCorR1] = fftPow(s,cLIPnormCorR1',nframes);
[cLIPpowIncR1, cLIPfreqIncR1] = fftPow(s,cLIPnormIncR1',nframes);
[cPECpowCorR1, cPECfreqCorR1] = fftPow(s,cPECnormCorR1',nframes);
[cPECpowIncR1, cPECfreqIncR1] = fftPow(s,cPECnormIncR1',nframes);
[cvPFCpowCorR1, cvPFCfreqCorR1] = fftPow(s,cvPFCnormCorR1',nframes);
[cvPFCpowIncR1, cvPFCfreqIncR1] = fftPow(s,cvPFCnormIncR1',nframes);
[cPGpowCorR1, cPGfreqCorR1] = fftPow(s,cPGnormCorR1',nframes);
[cPGpowIncR1, cPGfreqIncR1] = fftPow(s,cPGnormIncR1',nframes);
clearvars c9LnormCorR1 c9LnormIncR1 cMIPnormCorR1 cMIPnormIncR1 ...
    c8BnormCorR1 c8BnormIncR1 cdPFCnormCorR1 cdPFCnormIncR1 cLIPnormCorR1 ...
    cLIPnormIncR1 cPECnormCorR1 cPECnormIncR1 cvPFCnormCorR1 cvPFCnormIncR1 ...
    cPGnormCorR1 cPGnormIncR1