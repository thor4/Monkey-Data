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

% define trial timeline
signalt = -.504:1/srate:1.589; %505 (baseline) + 504 (nonzero sample) + 811 (delay) + 274 (match)=1589+505=2094ms
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1466; % in ms, 1466 = 505 (sample) + 811 (delay) + 150 (match)
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');

tic
for day=alldays(1) %cycle through all days
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
            r1=[]; %init var to hold trial numbers split by match location
            all_locs = unique(r1_locs(:,2))'; %row vec of all locations
            for location=all_locs %cycle through locations (needs to be row vec)
                temp_loc = sprintf('l%d_tri',location);
                r1.(temp_loc) = r1_locs(find(r1_locs(:,2)==location)); %returns trials at loc specified by 'location'
            end
            for location=all_locs %cycle through locations (needs to be row vec)
                temp_loc = sprintf('l%d_tri',location);
                temp_loc_pow = sprintf('l%d_pow',location);
                %init temp power var (chan x down-sampled time x frequency):
                temp_pow = zeros(length(areas),length(times2save),length(frex));
                for k=1:length(areas) %parse all channels
                    for fi=1:length(frex)
                        % save avg raw power for ea. chan x freq x down-sampled time avg across trials
                        temp_pow(k,:,fi) = squeeze( mean( abs( ansig(k,fi,times2saveidx,r1.(temp_loc)) ) .^2 ,4) );
                        % end up with down-sampled time vector of avg pow vals (187x1)
                    end
                end
                r1.(temp_loc_pow) = temp_pow; %save chan x time x frex avg pow over tri per loc
            end
        end
        if ~isempty(r2_locs) %ensures r2_locs isn't empty
            r2=[]; %init var to hold trial numbers split by match location
            all_locs = unique(r2_locs(:,2))'; %row vec of all locations
            for location=all_locs %cycle through locations (needs to be row vec)
                temp_loc = sprintf('l%d_tri',location);
                r2.(temp_loc) = r2_locs(find(r2_locs(:,2)==location)); %returns trials at loc specified by 'location'
            end
            for location=all_locs %cycle through locations (needs to be row vec)
                temp_loc = sprintf('l%d_tri',location);
                temp_loc_pow = sprintf('l%d_pow',location);
                %init temp power var (chan x down-sampled time x frequency):
                temp_pow = zeros(length(areas),length(times2save),length(frex));
                for k=1:length(areas) %parse all channels
                    for fi=1:length(frex)
                        % save avg raw power for ea. chan x freq x down-sampled time avg across trials
                        temp_pow(k,:,fi) = squeeze( mean( abs( ansig(k,fi,times2saveidx,r2.(temp_loc)) ) .^2 ,4) );
                        % end up with down-sampled time vector of avg pow vals (187x1)
                    end
                end
                r2.(temp_loc_pow) = temp_pow; %save chan x time x frex avg pow over tri per loc
            end
        end
        powloc_path = strcat(path,'%s\\%s\\%s\\%s%s%s.powloc.mat'); %build power x location data path
        pow_loc = sprintf(powloc_path,monkey,day{:},days{j},monkey,day{:},days{j}{1}(8:9));
        if exist('r1','var') && ~exist('r2','var') %only r1 in this session
            save(pow_loc,'r1','-v7.3') %save tri & pow x location struct
            clear r1
        elseif exist('r2','var') && ~exist('r1','var') %only r2 in this session
            save(pow_loc,'r2','-v7.3') %save tri & pow x location struct
            clear r2
        else %both r1 & r2 in this session
            save(pow_loc,'r1','r2','-v7.3') %save tri & pow x location struct
            clear r1 r2
        end
    end
end
%all correct trials and power split by stim location are now saved in a single file per session
toc

% 1423.38 sec for clark


% output: a single figure per channel with 3 subplots showing the avg power
% for each stimulus location for a specific rule across all correct trials

% work on the figures for each monkey
day_rules=[];

for day=alldays %cycle through all days
    dday = append('d',day{:}); %setup day for indexing
    for j=2:3
        if (j==3) && (monkey=="betty") %only one session for betty
            continue %skip rest of loop
        end
        trial_infoN = sprintf(trial_info_path,monkey,day{:},days{j}); %create full path to trial_info.mat
        load(trial_infoN,'trial_info'); %load trial_info for day's trials
        recording_infoN = sprintf(recording_info_path,monkey,day{:},days{j}); %create full path to recording_info.mat
        load(recording_infoN,'recording_info'); %load recording_info for day's trials
        areas = recording_info.area;
        powloc_path = strcat(path,'%s\\%s\\%s\\%s%s%s.powloc.mat'); %build power x location data path
        pow_loc = sprintf(powloc_path,monkey,day{:},days{j},monkey,day{:},days{j}{1}(8:9));
        load(pow_loc); %r1 or r2 or r1 & r2
        if exist('r1') %proceed with rule 1 trials
            %continue
            %day_rules.(dday).(days{j}) = 'rule1';
            pow_fields_r1 = fieldnames(r1)';
            pow_fields_r1{4:6}
            for locN=4:6 % STOPPED HERE, WORK ON FIGURE WITHIN EACH RULE
                r1.(pow_fields_r1{locN}) % avg pow over tri in chan x time x frex
            clear r1
        elseif exist('r2') %proceed with rule 2 trials
            %day_rules.(dday).(days{j}) = 'rule2';
            pow_fields_r2 = fieldnames(r2);
            clear r2
        else %proceed with both rule 1 and rule 2 trials
            %day_rules.(dday).(days{j}) = 'rules 1 & 2';
            pow_fields_r1 = fieldnames(r1);
            pow_fields_r2 = fieldnames(r2);
            clear r1 r2
        end
        %STOPPED HERE. NEED TO PARSE THROUGH THE r1/r2 structs AND PULL OUT
        %POWER AND PLOT IN FIGURES FOR CLARK. AFTERWARDS, START OVER WITH
        %BETTY BY MAKING NEW POWLOC FILES AND DOING FIGURES AGAIN
        for k=1:length(areas) %parse all channels
            for fi=1:length(frex)
                % save avg raw power for ea. chan x freq x down-sampled time avg across trials
                temp_pow(k,:,fi) = squeeze( mean( abs( ansig(k,fi,times2saveidx,r1.(temp_loc)) ) .^2 ,4) );
                % end up with down-sampled time vector of avg pow vals (187x1)
                %REVIEW THE BELOW FIGURE CODE
                % plotting the raw difference data
                
                %pow is saved as chan x time x frex in ie: r1.
                areaN = 2; %plot which area?
                clim = [0 90];

                pull out correct & incorrect power avg over trials
                cor = mean( monkey(monkeyN).(responses{1}).(areas{areaN}),3 ); 
                inc = mean( monkey(monkeyN).(responses{2}).(areas{areaN}),3 );
                compute the difference in power between the two conditions
                diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));
                diffmap = cor - inc;

                figure(9), clf
                subplot(221)
                imagesc(times2save,[],cor)
                hold on %vertical line
                yL = get(gca,'YLim'); line([0 500;0 500],yL,'Color','k','LineStyle',':'); 
                set(gca,'clim',clim,'ydir','n')
                set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
                xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
                lim = get(cbar,'Limits'); cbar.Ticks=lim;
                cbar.Label.String = 'Raw Power'; pos = cbar.Label.Position; 
                cbar.Label.Position=[pos(1)-1 pos(2)];
                title(sprintf('Monkey %d, Area %s, Resp Correct',monkeyN,areas{areaN}(2:end)));

                subplot(222)
                imagesc(times2save,[],inc)
                hold on %vertical line
                yL = get(gca,'YLim'); line([0 500;0 500],yL,'Color','k','LineStyle',':'); 
                set(gca,'clim',clim,'ydir','n')
                set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
                xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
                lim = get(cbar,'Limits'); cbar.Ticks=lim;
                cbar.Label.String = 'Raw Power'; pos = cbar.Label.Position; 
                cbar.Label.Position=[pos(1)-1 pos(2)];
                title(sprintf('Monkey %d, Area %s, Resp Incorrect',monkeyN,areas{areaN}(2:end)));

                subplot(223)
                imagesc(times2save,[], diffmap)
                hold on %vertical line
                yL = get(gca,'YLim'); line([0 500;0 500],yL,'Color','k','LineStyle',':'); 
                set(gca,'clim',[-mean(clim)/5 mean(clim)/5],'ydir','n')
                set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
                xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
                lim = get(cbar,'Limits'); cbar.Ticks=lim;
                cbar.Label.String = 'Raw Power'; pos = cbar.Label.Position; 
                cbar.Label.Position=[pos(1)-1 pos(2)];
                title(sprintf('Monkey %d, Area %s, Correct - Incorrect',monkeyN,areas{areaN}(2:end)));
            end
        end
    end
end