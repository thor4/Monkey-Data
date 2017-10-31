clear
%betty
%delay
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_6DR_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_6DR_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_8AD_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_8AD_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_8B_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_8B_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_dPFC_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_dPFC_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_LIP_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_LIP_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_PE_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_PE_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_PEC_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_PEC_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_PG_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_PG_delay_inc_rule1.mat')
%baseline
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_6DR_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_6DR_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_8AD_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_8AD_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_8B_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_8B_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_dPFC_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_dPFC_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_LIP_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_LIP_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_PE_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_PE_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_PEC_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_PEC_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_PG_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\baseline\betty_PG_baseline_inc_rule1.mat')

%calc single avg baseline signal across trials then across timepoints and
%subtract from delay period

b6DRnormCorR1 = b6DRdelayCorR1 - mean(mean(b6DRbaseCorR1'));
b6DRnormIncR1 = b6DRdelayIncR1 - mean(mean(b6DRbaseIncR1'));
b8ADnormCorR1 = b8ADdelayCorR1 - mean(mean(b8ADbaseCorR1'));
b8ADnormIncR1 = b8ADdelayIncR1 - mean(mean(b8ADbaseIncR1'));
b8BnormCorR1 = b8BdelayCorR1 - mean(mean(b8BbaseCorR1'));
b8BnormIncR1 = b8BdelayIncR1 - mean(mean(b8BbaseIncR1'));
bdPFCnormCorR1 = bdPFCdelayCorR1 - mean(mean(bdPFCbaseCorR1'));
bdPFCnormIncR1 = bdPFCdelayIncR1 - mean(mean(bdPFCbaseIncR1'));
bLIPnormCorR1 = bLIPdelayCorR1 - mean(mean(bLIPbaseCorR1'));
bLIPnormIncR1 = bLIPdelayIncR1 - mean(mean(bLIPbaseIncR1'));
bPEnormCorR1 = bPEdelayCorR1 - mean(mean(bPEbaseCorR1'));
bPEnormIncR1 = bPEdelayIncR1 - mean(mean(bPEbaseIncR1'));
bPECnormCorR1 = bPECdelayCorR1 - mean(mean(bPECbaseCorR1'));
bPECnormIncR1 = bPECdelayIncR1 - mean(mean(bPECbaseIncR1'));
bPGnormCorR1 = bPGdelayCorR1 - mean(mean(bPGbaseCorR1'));
bPGnormIncR1 = bPGdelayIncR1 - mean(mean(bPGbaseIncR1'));
clearvars b6DRdelayCorR1 b6DRbaseCorR1 b6DRdelayIncR1 b6DRbaseIncR1 ...
    b8ADdelayCorR1 b8ADbaseCorR1 b8ADdelayIncR1 b8ADbaseIncR1 ...
    b8BdelayCorR1 b8BbaseCorR1 b8BdelayIncR1 b8BbaseIncR1 ...
    bdPFCdelayCorR1 bdPFCbaseCorR1 bdPFCdelayIncR1 bdPFCbaseIncR1 ...
    bLIPdelayCorR1 bLIPbaseCorR1 bLIPdelayIncR1 bLIPbaseIncR1 ...
    bPEdelayCorR1 bPEbaseCorR1 bPEdelayIncR1 bPEbaseIncR1 ...
    bPECdelayCorR1 bPECbaseCorR1 bPECdelayIncR1 bPECbaseIncR1 ...
    bPGdelayCorR1 bPGbaseCorR1 bPGdelayIncR1 bPGbaseIncR1

clear
%clark
%delay
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_9L_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_9L_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_MIP_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_MIP_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_8B_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_8B_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_dPFC_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_dPFC_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_LIP_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_LIP_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_vPFC_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_vPFC_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_PEC_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_PEC_delay_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_PG_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\delay\clark_PG_delay_inc_rule1.mat')
%baseline
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_9L_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_9L_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_MIP_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_MIP_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_8B_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_8B_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_dPFC_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_dPFC_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_LIP_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_LIP_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_vPFC_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_vPFC_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_PEC_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_PEC_baseline_inc_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_PG_baseline_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\clark\baseline\clark_PG_baseline_inc_rule1.mat')

%calc single avg baseline signal across trials then across timepoints and
%subtract from delay period

c9LnormCorR1 = c9LdelayCorR1 - mean(mean(c9LbaseCorR1'));
c9LnormIncR1 = c9LdelayIncR1 - mean(mean(c9LbaseIncR1'));
cMIPnormCorR1 = cMIPdelayCorR1 - mean(mean(cMIPbaseCorR1'));
cMIPnormIncR1 = cMIPdelayIncR1 - mean(mean(cMIPbaseIncR1'));
c8BnormCorR1 = c8BdelayCorR1 - mean(mean(c8BbaseCorR1'));
c8BnormIncR1 = c8BdelayIncR1 - mean(mean(c8BbaseIncR1'));
cdPFCnormCorR1 = cdPFCdelayCorR1 - mean(mean(cdPFCbaseCorR1'));
cdPFCnormIncR1 = cdPFCdelayIncR1 - mean(mean(cdPFCbaseIncR1'));
cLIPnormCorR1 = cLIPdelayCorR1 - mean(mean(cLIPbaseCorR1'));
cLIPnormIncR1 = cLIPdelayIncR1 - mean(mean(cLIPbaseIncR1'));
cvPFCnormCorR1 = cvPFCdelayCorR1 - mean(mean(cvPFCbaseCorR1'));
cvPFCnormIncR1 = cvPFCdelayIncR1 - mean(mean(cvPFCbaseIncR1'));
cPECnormCorR1 = cPECdelayCorR1 - mean(mean(cPECbaseCorR1'));
cPECnormIncR1 = cPECdelayIncR1 - mean(mean(cPECbaseIncR1'));
cPGnormCorR1 = cPGdelayCorR1 - mean(mean(cPGbaseCorR1'));
cPGnormIncR1 = cPGdelayIncR1 - mean(mean(cPGbaseIncR1'));
clearvars c9LdelayCorR1 c9LbaseCorR1 c9LdelayIncR1 c9LbaseIncR1 ...
    cMIPdelayCorR1 cMIPbaseCorR1 cMIPdelayIncR1 cMIPbaseIncR1 ...
    c8BdelayCorR1 c8BbaseCorR1 c8BdelayIncR1 c8BbaseIncR1 ...
    cdPFCdelayCorR1 cdPFCbaseCorR1 cdPFCdelayIncR1 cdPFCbaseIncR1 ...
    cLIPdelayCorR1 cLIPbaseCorR1 cLIPdelayIncR1 cLIPbaseIncR1 ...
    cvPFCdelayCorR1 cvPFCbaseCorR1 cvPFCdelayIncR1 cvPFCbaseIncR1 ...
    cPECdelayCorR1 cPECbaseCorR1 cPECdelayIncR1 cPECbaseIncR1 ...
    cPGdelayCorR1 cPGbaseCorR1 cPGdelayIncR1 cPGbaseIncR1