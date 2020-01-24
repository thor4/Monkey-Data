%%%% Crawler for raw data %%%%
%Will crawl specific file location based on PC used: koko, home or lab

%doc: https://www.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html

function [M1,M2] = craw(path,varargin)
    p = inputParser; %create inputParser object to check inputs
    %define default optional parameter values
%     path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\'
    defaultMonkey = 'betty'; monkeys = [ "betty", "clark" ];
    %validate monkey exists and is accurate
    checkMonkey = @(x) any(validatestring(x,monkeys));
    defaultSession = 'all'; 
    days_betty = { '090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001' };
    days_clark = { '060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221' };
    %validate session exists or if user wants all sessions
    validSession = [ days_betty days_clark 'all' ];
    checkSession = @(x) any(validatestring(x,validSession)); 
    defaultGood x= 'yes'; validGood = { 'yes','no' }; %yes has no artifacts
    checkGood = @(x) any(validatestring(x,validGood));
    defaultStable = 'no'; %don't account for stability
    validStable = { 'yes','no','transition' }; %yes stable perf or trial during transition bet rules
    checkStable = @(x) any(validatestring(x,validStable));
    defaultBehResp = 1; %1 correct trials
    checkBehResp = @(x) ismember(x,[0,1]); %0 incorrect trials
    defaultRule = 1; checkRule = @(x) ismember(x,[1,2]); %1 identity, 2 location
    defaultEpoch = 'delay'; validEpoch = { 'base','sample','delay','match' }; 
    checkEpoch = @(x) any(validatestring(x,validEpoch));
    %add required input + optional parameter values & verify datatype
%     addRequired(p,'data',@isstruct);
    addRequired(p,'path',@ischar);
    addOptional(p,'Monkey',defaultMonkey,checkMonkey)
    addOptional(p,'Session',defaultSession,checkSession)
    addOptional(p,'Good',defaultGood,checkGood)
    addOptional(p,'Stable',defaultStable,checkStable)
    addOptional(p,'BehResp',defaultBehResp,checkBehResp)
    addOptional(p,'Rule',defaultRule,checkRule)
    addOptional(p,'Epoch',defaultEpoch,checkEpoch)
    parse(p,varargin{:}) %parse all the inputs
    %report parsed inputs back to user for confirmation
    disp(['Crawling: ',p.Results.data])
    if ~isempty(p.UsingDefaults) %check if using any default values
        disp('Using defaults: ')
        disp(p.UsingDefaults)
    end
    trial_info_path = strcat(path,'%s\\%s\\%s\\trial_info.mat'); %build trial_info path
    recording_info_path = strcat(path,'%s\\%s\\%s\\recording_info.mat'); %build recording_info path
    lfp_path = strcat(path,'%s\\%s\\%s\\%s%s%s.%04d.mat'); %build lfp raw data path
    monkeys = [ "betty", "clark" ];
    betty = { days_betty, "session01" };
    clark = { days_clark, "session02", "session03" };
    z = 1;
    M1=0;
    %betty
    
    idx = 0; % counter

    for i=1:length(betty{1})
        trial_infoN = sprintf(trial_info_path, monkeys(1), betty{1}{i}, betty{2}); %create full path to trial_info.mat
        load(trial_infoN); %load trial_info for day's trials
        recording_infoN = sprintf(recording_info_path, monkeys(1), betty{1}{i}, betty{2}); %create full path to recording_info.mat
        load(recording_infoN); %load recording_info for day's trials
        %% left off here %%
        for j=1:trial_info.numTrials
            if (trial_info.good_trials(j) == 1) && ...%no artifacts
                    (trial_info.rule(j) == 1) && ... %only identity rule
                    (trial_info.BehResp(j) == 1)  %only correct resp
                % extract rule, response and rt in 1x3 matrix
                idx = idx + 1;
                trial_lfp_myfilename = sprintf(lfp_path, monkeys(1), betty{1}{i}, betty{2}, monkeys(1), betty{1}{i}, betty{2}{1}(8:9), j);
                load(trial_lfp_myfilename);
                chan = length(recording_info.area); %get number of channels
    %             unique(recording_info.area);
                %define delay period length
                delay_period = trial_info.MatchOnset(j)-trial_info.CueOffset(j);
                %extract LFP during delay period: take one second after sample 
                %goes away up through 810ms afterwards. insert into delay
                %matrix with trial count as 3rd dimension
                delay(:,:,idx) = lfp_data(:,trial_info.CueOffset(j)+1:trial_info.CueOffset(j)+810);
                obj(idx) = trial_info.CueObj(j);
            end
        end
        m2_obj_delay.day(i).delay = delay;
        m2_obj_delay.day(i).objects = obj;
        clearvars delay obj
        idx=0; %reset counter for new day
    end

    
    
    fieldsBettyDays = fieldnames(bettyGoodStableTrials);
    for i = 1:numel(fieldsBettyDays)
        fieldsBettyTrials = fieldnames(bettyGoodStableTrials.(fieldsBettyDays{i}));
        recordingRegion = bettyGoodStableTrials.(fieldsBettyDays{i}).recording_info.cortex;
        recordingArea = bettyGoodStableTrials.(fieldsBettyDays{i}).recording_info.area;
        for j = 1:numel(fieldsBettyTrials)-2
            if BehResp && Rule && region %provide BehResp,Rule,region slice
                for k = 1:numel(recordingRegion)                
                    if (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == BehResp) && ...
                            (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == Rule) && ...
                            (recordingRegion(k) == region) %only look at trials from requested region
                        M2(z,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                        M2(z,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                        z = z + 1;
                    end
                end
            elseif BehResp && Rule && area %provide BehResp,Rule,area slice
                for k = 1:numel(recordingArea)
                   if (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == BehResp) && ...
                            (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == Rule) && ...
                            (recordingArea(k) == area) %only look at trials from requested area
                        M2(z,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                        M2(z,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                        z = z + 1;
                   end 
                end
            else
                print("Please enter either BehResp,Rule,Region combo or BehResp,Rule,area");
            end
        end
    end
end
    %clark
%     tic
%     fieldsClarkDays = fieldnames(clarkGoodStableTrials);
%     for i = 1:numel(fieldsClarkDays)
%         fieldsClarkSessions = fieldnames(clarkGoodStableTrials.(fieldsClarkDays{i}));
%         for j = 1:numel(fieldsClarkSessions)
%             fieldsClarkTrials = fieldnames(clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}));
%             recordingRegion = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).recording_info.cortex;
%             for k = 1:numel(fieldsClarkTrials)-2
%                 for l = 1:numel(recordingRegion)
%                     if (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
%                             (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
%                             (recordingRegion(l) == 'F') %only look at correct rule 1 frontal trials
%                         correctFrontalClark(zfc,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
%                         correctFrontalClark(zfc,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
%                         zfc = zfc + 1;
%                     elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
%                             (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
%                             (recordingRegion(l) == 'F') %only look at incorrect rule 1 frontal trials
%                         incorrectFrontalClark(zfi,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
%                         incorrectFrontalClark(zfi,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
%                         zfi = zfi + 1;
%                     elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
%                             (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
%                             (recordingRegion(l) == 'P') %only look at correct rule 1 parietal trials
%                         correctParietalClark(zpc,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
%                         correctParietalClark(zpc,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
%                         zpc = zpc + 1;
%                     elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
%                             (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
%                             (recordingRegion(l) == 'P') %only look at incorrect rule 1 parietal trials
%                         incorrectParietalClark(zpi,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
%                         incorrectParietalClark(zpi,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
%                         zpi = zpi + 1;
%                     end
%                 end
%             end
%         end
%     end
% end
% toc