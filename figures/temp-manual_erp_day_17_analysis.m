%% Segment the trials

% navigate to monkey B day 17 (suspicious):
% OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\betty\090917\session01

load('trial_info.mat')

idx = 0; %init counter
for k=1:trial_info.numTrials %parse all trials for day
%stable not specified, give stable perf and transition trials
    if (trial_info.good_trials(k) == 1) && ...%artifacts/none 0/1
            (trial_info.BehResp(k) == 1) && ... %correct(1)/incorrect(0)
            (trial_info.rule(k) == 1) %identity(1)/location(2)
        idx = idx + 1;  
        trial_lfp = sprintf('betty09091701.%04d.mat',k);
        load(trial_lfp,'lfp_data');
        base = lfp_data(:,floor(trial_info.CueOnset(k))-504:floor(trial_info.CueOnset(k))-1);
        sample = lfp_data(:,floor(trial_info.CueOnset(k)):floor(trial_info.CueOnset(k))+504);
        delay = lfp_data(:,floor(trial_info.CueOffset(k)):floor(trial_info.CueOffset(k))+810);
        match = lfp_data(:,floor(trial_info.MatchOnset(k)):floor(trial_info.MatchOnset(k))+199);
        lfp(:,:,idx) = cat(2,base,sample,delay,match);
    end
end

%504 samples in baseline
%505 samples in cue
%811 samples in delay
%200 samples in match
%2020 total samples across all chans
%gives all trials cut up and stitched together, ready for ERPing

%% ERP 

