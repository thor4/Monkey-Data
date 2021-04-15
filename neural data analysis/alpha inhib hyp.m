% New analysis on the Alpha Inhibition Hypothesis 4/8/21
% 
% Question: does the power during different stages of the task depend on the location of the stimulus?
% 
% use analytic signal created from decompose_all_trial_lfps_to_as.m
%
% identify for every trial the rule and where the stimulus was presented, 
% which corner of the triangle. the triangle was flipped per day
% [0,1,2] is one orientation of the triangle and [3,4,5] is the other
%
% output: a single figure per channel with 3 subplots showing the avg power
% for each stimulus location for a specific rule across all correct trials



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
lfp_path = strcat(path,'%s\\%s\\%s\\%s%s%s.ansig.mat'); %build lfp raw data path

srate = 1000; % 1,000Hz
min_freq = 4; %in Hz (need several cycles in an epoch, these epochs are 500ms min so 4Hz = 2 cycles)
max_freq = 100; %nothing above 100
num_frex = 35; %50 for 200Hz, 35 for 100Hz, better for statistics mult comp corr, less smooth spectrogram
min_fwhm = .400; % in seconds (350ms)
max_fwhm = .100; % in seconds (75ms)
%there are N/2+1 frequencies between 0 and srate/2:
frex = logspace(log10(min_freq),log10(max_freq),num_frex); %total num of freq's

%504 samples in baseline
%505 samples in cue
%811 samples in delay
%274 samples in match (check this)
%2094 total samples across all chans, 35 total frequencies

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
        day_ansig = sprintf(lfp_path,monkey,day{:},days{j},monkey,day{:},days{j}{1}(8:9));
        %load analytic signal for all trials across all chans in a single variable 'ansig':
        load(day_ansig,'ansig'); %chan (recording_info.numChannels) x frex (35) x time (2079) x trial (trial_info.numTrials)
        %size(ansig) %verify above
              
        for k=1:recording_info.numChannels %parse all channels
            r1_locs = []; r2_locs = []; %init location vars for ea rule
            %first, compile trials of interest
            for l=1:trial_info.numTrials
                if (trial_info.BehResp(l)==1) && (trial_info.rule(l)==1) %cor + rule1
                    if (trial_info.CueObj(l)==trial_info.MatchObj1(l)) %match identity for rule 1
                        r1_locs(end+1,:) = [l trial_info.MatchPos1(l)]; %save trial # & position
                    elseif (trial_info.CueObj(l)==trial_info.MatchObj2(l))
                        r1_locs(end+1,:) = [l trial_info.MatchPos2(l)];
                    end
                elseif (trial_info.BehResp(l)==1) && (trial_info.rule(l)==2) %cor + rule2
                    if (trial_info.CueLoc(l)==trial_info.MatchPos1(l)) %match location for rule 2
                        r2_locs(end+1,:) = [l trial_info.MatchPos1(l)];
                    elseif (trial_info.CueLoc(l)==trial_info.MatchPos2(l))
                        r2_locs(end+1,:) = [l trial_info.MatchPos2(l)];
                    end
                end
            end
            if ~isempty(r1_locs) %ensures r1_locs isn't empty
                r1_locs_split=[]; %init var to hold trial numbers split by location
                all_locs = unique(r1_locs(:,2))'; %row vec of all locations
                for location=all_locs %cycle through locations (needs to be row vec)
                    temp_loc = sprintf('l%d',location);
                    r1_locs_split.(temp_loc) = r1_locs(find(r1_locs(:,2)==location)); %returns trials at loc specified by 'location'
                end
            end
            if ~isempty(r2_locs) %ensures r2_locs isn't empty
                r2_locs_split=[]; %init var to hold trial numbers split by location
                all_locs = unique(r2_locs(:,2))'; %row vec of all locations
                for location=all_locs %cycle through locations (needs to be row vec)
                    temp_loc = sprintf('l%d',location);
                    r2_locs_split.(temp_loc) = r2_locs(find(r2_locs(:,2)==location)); %returns trials at loc specified by 'location'
                end
            end
        %STOPPED HERE
        %PULLED OUT SOME TRIALS OF INTEREST: COR + RULE 1/2 (sep var)
        %NOW CONTINUE TO COMPUTE AVG POWER ACROSS TRIALS PER STIMULUS LOCATION
            
            
            
            
            for fi=1:length(frex)
                % & save avg power per freq component avg across trials
                pow(fi) = mean( abs( ansig(k,fi,:,l) as.^2 ),4 );
                %need to identify stim location based on match loc & rule
                % store dB-norm'd down-sampled power for each frequency in freq
                % x time x trials
%                 dbpow(fi,:,:) = 10*log10( abs( as(times2saveidx,:) ) .^2 ./basePow); 
                % save raw power for ea. freqidx x down-sampled time x
                % trial
%                 pow(fi,:,:) = abs( as(times2saveidx,:) ) .^2;
                % mean( abs( as_ ).^2, 2);
%                     clear as % start anew with these var's ea. loop
                if 
                    continue
                elseif (trial_info.rule(l)==2) %rule2 (location)
                    continue
                end
            end
                    
                
            %save downsampled power as chan x freqidx x time x trials
%             data.(monkey{:})(i).(dday{:}).power(chan,:,:,:) = pow; 
            %save avg baseline power for all frex as chan x freqidx
            data.(monkey{:})(i).(dday{:}).basepow(chan,:) = basePow; % as is now a time x trial complex matrix
                ansig(k,fi,:,:) = as;
        end
        lfp_as_path = strcat(path,'%s\\%s\\%s\\%s%s%s.ansig.mat'); %build lfp all analytic signals path
        lfp_as = sprintf(lfp_as_path,monkey,day{:},days{j},monkey,day{:},days{j}{1}(8:9));
        save(lfp_as,'ansig','-v7.3') %save chan x freq (35) x time (2079) x trials ansig matrix
        clear ansig
%         % testing:
%         find(ansig(:,:,:,:)==0);
%         size(ansig)
    end
end




