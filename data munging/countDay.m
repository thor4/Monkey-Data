function [lengths, idx] = countDay(path,monkey,day,good,stable,behResp,rule,epoch)
%%%% Crawler for raw data %%%%
%Will crawl specific file location based on PC used: koko, home or lab

%doc: https://www.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html

% monkey: [ 'betty', 'clark' ]
% day:(betty) [   '090615', '090616', '090617', '090618', '090622', 
%                 '090625', '090626', '090629', '090701', '090702', 
%                 '090706', '090708', '090709', '090901', '090903', 
%                 '090916', '090917', '090921', '090923', '090924', 
%                 '090928', '090929', '090930', '091001']
%     (clark) [   '060328', '060406', '060411', '060414', '060426', 
%                 '060427', '060428', '060502', '060503', '060509', 
%                 '060511', '060531', '060601', '060602', '060824', 
%                 '060825', '060831', '060907', '061212', '061213', 
%                 '061214', '061215', '061221']
% good: [ 0(artifacts), 1(no artifacts) ]
% stable: [ 0(transition), 1(stable performance), 2(both) ]
% behResp: [0(incorrect), 1(correct) ]
% rule: [ 1(identity), 2(location) ]
% epoch: [ 'base', 'sample', 'delay', 'match', 'all' ]

%   homepc:
%   path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\'
%   labpc:
%   path = 'C:\\Users\\bryan\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\'

% Version 1.0   (2020 January)


    p = inputParser; %create inputParser object to check inputs
    %define default optional parameter values
%   test  monkey = "betty" 
%   test day = '090615'
%   test good = 1;
%   test stable = 2;
%   test behResp = 1;
%   test rule = 1;
%   test epoch = 'delay';
    monkeys = { 'betty', 'clark' };
    %validate monkey exists and is accurate
    checkMonkey = @(x) any(strcmp(x,monkeys));
    days_betty = { '090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001' };
    days_clark = { '060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221' };
    %validate session exists or if user wants all sessions
    validDay = [ days_betty days_clark 'all' ];
    checkDay = @(x) any(validatestring(x,validDay)); 
    checkGood = @(x) ismember(x,[0,1]); %0 bad trials, 1 good/no artifacts
    checkBehResp = @(x) ismember(x,[0,1]); %0 incorrect trials, 1 correct
    checkRule = @(x) ismember(x,[1,2]); %1 identity, 2 location
    validEpoch = { 'base','sample','delay','match','all' }; 
    checkEpoch = @(x) any(strcmp(x,validEpoch));
    checkStable = @(x) ismember(x,[0,1,2]); %0 transition, 1 stable perf, 2 all
    %add required input + optional parameter values & verify datatype
    addRequired(p,'path',@ischar);
    addRequired(p,'monkey',checkMonkey)
    addRequired(p,'day',checkDay)
    addRequired(p,'good',checkGood)
    addRequired(p,'stable',checkStable)
    addRequired(p,'behResp',checkBehResp)
    addRequired(p,'rule',checkRule)
    addRequired(p,'epoch',checkEpoch)
    parse(p,path,monkey,day,good,stable,behResp,rule,epoch) %parse all the inputs
    %report parsed inputs back to user for confirmation
    disp(['Monkey: ',p.Results.monkey]), disp(['Day: ',p.Results.day])
    disp(['Good: ',num2str(p.Results.good)]), disp(['Stable: ',num2str(p.Results.stable)])
    disp(['Response: ',num2str(p.Results.behResp)]), disp(['Rule: ',num2str(p.Results.rule)])
    disp(['Epoch: ',p.Results.epoch])
    if ~isempty(p.UsingDefaults) %check if using any default values
        disp('Using defaults: ')
        disp(p.UsingDefaults)
    end
    trial_info_path = strcat(path,'%s\\%s\\%s\\trial_info.mat'); %build trial_info path
%     recording_info_path = strcat(path,'%s\\%s\\%s\\recording_info.mat'); %build recording_info path
    lfp_path = strcat(path,'%s\\%s\\%s\\%s%s%s.%04d.mat'); %build lfp raw data path

    if monkey=="betty"
        days = { days_betty, "session01" };
    else
        days = { days_clark, "session02", "session03" };
    end
        
    idx(1,1:2) = 0; %init counter
    lengths = zeros(2000,2); %init length mat
    if ~(string(day) == "all") %single day
        for j=2:3
            if (j==3) && (monkey=="betty") %only one session for betty
                continue %skip rest of loop
            end
            trial_infoN = sprintf(trial_info_path, monkey,day,days{j}); %create full path to trial_info.mat
            load(trial_infoN,'trial_info'); %load trial_info for day's trials
            for k=1:trial_info.numTrials %parse all trials for day
                if ~(stable==2) && ...%stable is specified
                        (trial_info.good_trials(k) == good) && ...%artifacts/none
                        (trial_info.stable_trials(k) == stable) && ...%stabile/transition
                        (trial_info.BehResp(k) == behResp) && ... %correct/incorrect
                        (trial_info.rule(k) == rule) %identify/location
                    idx(j-1) = idx(j-1) + 1;                        
                    trial_lfp = sprintf(lfp_path,monkey,day,days{j},monkey,day,days{j}{1}(8:9),k);
                    load(trial_lfp,'lfp_data');
                    switch epoch
                        case 'base'
                            %take 501ms before sample onset up until 1ms 
                            %before sample onset
                            lengths(idx(j-1),j-1) = length(lfp_data(:,1:floor(trial_info.CueOnset(k))-1)); %len of baseline in ms
                        case 'sample'
                            %take first ms sample is turned on up until 
                            %509ms of sample. turned off around 510 or so 
                            lengths(idx(j-1),j-1) = length(lfp_data(:,floor(trial_info.CueOnset(k)):floor(trial_info.CueOffset(k))-1)); %len of sample in ms
                        case 'delay'
                            %take one second after sample goes away up 
                            %through 810ms afterwards. 
                            lengths(idx(j-1),j-1) = length(lfp_data(:,floor(trial_info.CueOffset(k)):floor(trial_info.MatchOnset(k))-1)); %len of delay in ms
                        case 'match'
                            %take first ms match is turned on up until 
                            %200ms of match. 
                            lengths(idx(j-1),j-1) = length(lfp_data(:,floor(trial_info.MatchOnset(k)):end)); %len of match in ms
                        case 'all'
                            lengths(idx(j-1),j-1) = length(lfp_data); %len of entire trial in ms
                        otherwise
                            warning('no such epoch exists')
                    end
                elseif (stable==2) && ...%stable not specified, give stable perf and transition trials
                        (trial_info.good_trials(k) == good) && ...%artifacts/none
                        (trial_info.BehResp(k) == behResp) && ... %correct/incorrect
                        (trial_info.rule(k) == rule) %identify/location
                    idx(j-1) = idx(j-1) + 1;  
                    trial_lfp = sprintf(lfp_path,monkey,day,days{j},monkey,day,days{j}{1}(8:9),k);
                    load(trial_lfp,'lfp_data');
                    switch epoch
                        case 'base'
                            %take 501ms before sample onset up until 1ms 
                            %before sample onset
                            lengths(idx(j-1),j-1) = length(lfp_data(:,1:floor(trial_info.CueOnset(k))-1)); %len of baseline in ms
                        case 'sample'
                            %take first ms sample is turned on up until 
                            %509ms of sample. turned off around 510 or so 
                            lengths(idx(j-1),j-1) = length(lfp_data(:,floor(trial_info.CueOnset(k)):floor(trial_info.CueOffset(k))-1)); %len of sample in ms
                        case 'delay'
                            %take one second after sample goes away up 
                            %through 810ms afterwards. 
                            lengths(idx(j-1),j-1) = length(lfp_data(:,floor(trial_info.CueOffset(k)):floor(trial_info.MatchOnset(k))-1)); %len of delay in ms
                        case 'match'
                            %take first ms match is turned on up until 
                            %200ms of match. 
                            lengths(idx(j-1),j-1) = length(lfp_data(:,floor(trial_info.MatchOnset(k)):end)); %len of match in ms
                        case 'all'
                            lengths(idx(j-1),j-1) = length(lfp_data); %len of entire trial in ms
                        otherwise
                            warning('no such epoch exists')
                    end
                end
            end
        end
    else
        warning('can only provide single day data')
    end
%     else % go through all days, need to update 'lengths' to store all days and sessions
%         for i=1:length(days{1}) %all days
%             for j=2:3
%                 if (j==3) && (monkey=="betty") %only one session for betty
%                     continue
%                 end
%                 trial_infoN = sprintf(trial_info_path, monkey,days{1}{i},days{j}); %create full path to trial_info.mat
%                 load(trial_infoN,'trial_info'); %load trial_info for day's trials
% %                 [trial_info, ~] = loadinfo(days{1}{i},days{j}); %load trial & rec info for day
%                 idx = 0; %init counter
%                 for k=1:trial_info.numTrials %parse all trials for day
%                     if ~(stable==2) %stable is specified
%                         if (trial_info.good_trials(k) == good) && ...%artifacts/none
%                                 (trial_info.stable_trials(k) == stable) && ...%stabile/transition
%                                 (trial_info.BehResp(k) == behResp) && ... %correct/incorrect
%                                 (trial_info.rule(k) == rule) %identify/location
% %                             idx = idx + length(recording_info.area); %inc by # of channels
%                             trial_lfp = sprintf(lfp_path,monkey,days{1}{i},days{j},monkey,days{1}{i},days{j}{1}(8:9),k);
%                             load(trial_lfp,'lfp_data');
%                             %%verify the timing below with the ERP analysis
%                             %done before. Ensure same time periods are
%                             %taken
%                             switch epoch
%                                 case 'base'
%                                     %take 501ms before sample onset up
%                                     %until 1ms before sample onset =
%                                     %lfp(chan,500,trial)
%                                     %lfp(:,:,idx) = lfp_data(:,trial_info.CueOnset(k)-501:trial_info.CueOnset(k)-1);
%                                     lengths(idx) = length(lfp_data(:,1:floor(trial_info.CueOnset(k))-1)); %len of baseline in ms
%                                 case 'sample'
%                                     %take first ms sample is turned on up
%                                     %until 509ms of sample. turned off
%                                     %around 510 or so = lfp(chan,510,trial)
%                                     %lfp(:,:,idx) = lfp_data(:,trial_info.CueOnset(k):trial_info.CueOnset(k)+509);
%                                     lengths(idx) = length(lfp_data(:,floor(trial_info.CueOnset(k)):floor(trial_info.CueOffset(k))-1)); %len of sample in ms
%                                 case 'delay'
%                                     %take one second after sample goes away
%                                     %up through 810ms afterwards. insert in
%                                     %matrix with trial count as 3rd
%                                     %dimension = lfp(chan,810,trial)
%                                     %lfp(:,:,idx) = lfp_data(:,trial_info.CueOffset(k)+1:trial_info.CueOffset(k)+810);
%                                     lengths(idx) = length(lfp_data(:,floor(trial_info.CueOffset(k)):floor(trial_info.MatchOnset(k))-1)); %len of delay in ms
%                                 case 'match'
%                                     %take first ms match is turned on up
%                                     %until 200ms of match. =
%                                     %lfp(chan,210,trial)
%                                     %lfp(:,:,idx) = lfp_data(:,trial_info.MatchOnset(k):trial_info.MatchOnset(k)+209);
%                                     lengths(idx) = length(lfp_data(:,floor(trial_info.MatchOnset(k)):end)); %len of match in ms
%                                 case 'all'
%                                     %lfp(:,:,idx) = lfp_data;
%                                     lengths(idx) = length(lfp_data); %len of entire trial in ms
%                                 otherwise
%                                     warning('no such epoch exists')
%                             end
%                             idx = idx + 1;
%                             %%save entire trial_info and rec_info per day in struct
%                         end
%                     else %stable not specified, give stable perf and transition trials
%                         if (trial_info.good_trials(k) == good) && ...%artifacts/none
%                                 (trial_info.BehResp(k) == behResp) && ... %correct/incorrect
%                                 (trial_info.rule(k) == rule) %identify/location
%                             idx = idx + 1;
%                             trial_lfp = sprintf(lfp_path,monkey,days{1}{i},days{j},monkey,days{1}{i},days{j}{1}(8:9),k);
%                             load(trial_lfp,'lfp_data');
%                             %%load up trials accordingly
%                             switch epoch
%                                 case 'base'
%                                     %take 501ms before sample onset up
%                                     %until 1ms before sample onset =
%                                     %lfp(chan,500,trial)
%                                     %lfp(:,:,idx) = lfp_data(:,trial_info.CueOnset(k)-501:trial_info.CueOnset(k)-1);
%                                     lengths(idx) = length(lfp_data(:,1:floor(trial_info.CueOnset(k))-1)); %len of baseline in ms
%                                 case 'sample'
%                                     %take first ms sample is turned on up
%                                     %until 509ms of sample. turned off
%                                     %around 510 or so = lfp(chan,510,trial)
%                                     %lfp(:,:,idx) = lfp_data(:,trial_info.CueOnset(k):trial_info.CueOnset(k)+509);
%                                     lengths(idx) = length(lfp_data(:,floor(trial_info.CueOnset(k)):floor(trial_info.CueOffset(k))-1)); %len of sample in ms
%                                 case 'delay'
%                                     %take one second after sample goes away
%                                     %up through 810ms afterwards. insert in
%                                     %matrix with trial count as 3rd
%                                     %dimension = lfp(chan,810,trial)
%                                     %lfp(:,:,idx) = lfp_data(:,trial_info.CueOffset(k)+1:trial_info.CueOffset(k)+810);
%                                     lengths(idx) = length(lfp_data(:,floor(trial_info.CueOffset(k)):floor(trial_info.MatchOnset(k))-1)); %len of delay in ms
%                                 case 'match'
%                                     %take first ms match is turned on up
%                                     %until 200ms of match. =
%                                     %lfp(chan,210,trial)
%                                     %lfp(:,:,idx) = lfp_data(:,trial_info.MatchOnset(k):trial_info.MatchOnset(k)+209);
%                                     lengths(idx) = length(lfp_data(:,floor(trial_info.MatchOnset(k)):end)); %len of match in ms
%                                 case 'all'
%                                     %lfp(:,:,idx) = lfp_data;
%                                     lengths(idx) = length(lfp_data); %len of entire trial in ms
%                                 otherwise
%                                     warning('no such epoch exists')
%                             end
%                         end
%                     end
%                 end
%             end
%         end
end