% New analysis 2/10/21
% Work on raw time-frequency power from wavelets_power.m script
% Goal: save analytic signal for every trial
% homepc: G:\monkey_data


%% Step 3: Create analytic signal from lfp per day

% FIRST extract_and_save_LFPs.m
% lfp signal is chan x time x trials
% ex: 'lfp' from clark06121303.all.mat

% SECOND decompose_all_trial_lfps_to_as.m

%setup loops to crawl through each day, load the LFPs and create the
%analytic signals then save them to disk
