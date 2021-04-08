%% Create analytic signal from lfp per day

% FIRST extract_and_save_LFPs.m [this file]
% lfp signal is chan x time x trials
% ex: 'lfp' from clark06121303.all.mat

% SECOND decompose_all_trial_lfps_to_as.m 

%setup loops to crawl through each day, load the LFPs and create the
%analytic signals then save them to disk

%% Extract all LFPs from both monkeys, all trials and save as separate file in each session

% init vars, run this twice (once for each monkey):
%   homepc:
path = 'G:\\monkey_data\\';
%   labpc:
%   path = 'C:\\Users\\bryan\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\'
monk = 1; %1 = clark, 2 = betty
days_betty = { '090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001' };
days_clark = { '060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221' };
if monk==1
    monkey='clark';
    alldays = days_clark;
    days = { days_clark, "session02", "session03" };
else
    monkey='betty';
    alldays = days_betty;
    days = { days_betty, "session01" };
end
trial_info_path = strcat(path,'%s\\%s\\%s\\trial_info.mat'); %build trial_info path
recording_info_path = strcat(path,'%s\\%s\\%s\\recording_info.mat'); %build recording_info path
lfp_path = strcat(path,'%s\\%s\\%s\\%s%s%s.%04d.mat'); %build lfp raw data path
for day=alldays %cycle through all days
    for j=2:3
        if (j==3) && (monkey=="betty") %only one session for betty
            continue %skip rest of loop
        end
        trial_infoN = sprintf(trial_info_path,monkey,day{:},days{j}); %create full path to trial_info.mat
        load(trial_infoN,'trial_info'); %load trial_info for day's trials
        recording_infoN = sprintf(recording_info_path,monkey,day{:},days{j}); %create full path to recording_info.mat
        load(recording_infoN,'recording_info'); %load recording_info for day's trials
        areas = recording_info.area;
        for k=1:trial_info.numTrials %parse all trials for day                       
            trial_lfp = sprintf(lfp_path,monkey,day{:},days{j},monkey,day{:},days{j}{1}(8:9),k);
            load(trial_lfp,'lfp_data'); 
            %pull out epochs:
            %take 504ms before sample onset up until 1ms before sample onset 
            base = lfp_data(:,floor(trial_info.CueOnset(k))-504:floor(trial_info.CueOnset(k))-1);
            %take first ms sample is turned on up until 504ms of sample
            sample = lfp_data(:,floor(trial_info.CueOnset(k)):floor(trial_info.CueOnset(k))+504);
            %take first ms sample goes away up until 810ms after
            delay = lfp_data(:,floor(trial_info.CueOffset(k)):floor(trial_info.CueOffset(k))+810);
            %take first ms match is turned on up until 258ms of match
            match = lfp_data(:,floor(trial_info.MatchOnset(k)):floor(trial_info.MatchOnset(k))+258);
            lfp(:,:,k) = cat(2,base,sample,delay,match) .* 1e6; %patch epochs together
            %and convert to µV (1V = 10^6µV = 1,000,000µV) 
        end
        lfp_all_path = strcat(path,'%s\\%s\\%s\\%s%s%s.all.mat'); %build lfp all raw data path
        lfp_all = sprintf(lfp_all_path,monkey,day{:},days{j},monkey,day{:},days{j}{1}(8:9));
        save(lfp_all,'lfp') %save chan x time (2079) x trials lfp matrix
        clear lfp 
    end
end %all trial lfp's are now saved in a single file per session