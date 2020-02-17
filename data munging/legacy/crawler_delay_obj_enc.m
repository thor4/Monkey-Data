%extract delay-period data for self-organizing-map modeling object encoding
clear
monkeys = string({'betty', 'clark'});
days_betty = string({'090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001'});
days_clark = string({'060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221'});
betty = { days_betty, string('session01') };
clark = { days_clark, string('session02'), string('session03') };
% total number of trials indexed by monkey, day, session
path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\trial_info.mat';
record_path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\recording_info.mat';
lfp_path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\%s%s%s.%04d.mat';

idx = 0; % counter

for i=1:length(betty{1})
    myfilename = sprintf(path, monkeys(1), betty{1}{i}, betty{2}); %create full path to trial_info.mat
    load(myfilename); %load trial_info for day's trials
    myfilename_record = sprintf(record_path, monkeys(1), betty{1}{i}, betty{2}); %create full path to recording_info.mat
    load(myfilename_record); %load recording_info for day's trials
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

i=1; %day 1
myfilename = sprintf(path, monkeys(1), betty{1}{i}, betty{2}); %create full path to trial_info.mat
    load(myfilename); %load trial_info for day's trials
    myfilename_record = sprintf(record_path, monkeys(1), betty{1}{i}, betty{2}); %create full path to recording_info.mat
    load(myfilename_record); %load recording_info for day's trials
    fchans = {'6DR', '8AD', '8B', 'dPFC'}; pchans = {'LIP', 'PE', 'PEC', 'PG'}; %m2
    areas = {'6DR', '8AD', '8B', 'dPFC', 'LIP', 'PE', 'PEC', 'PG'}; %monkey2
    combos = combnk(recording_info.area,2); %pull out all chan pairs
    areaF=[]; areaP=[]; %initialize areas
    for comboN=1:size(combos,1) %total number of combinations
        for areaN=1:numel(areas) %which area pair
            if (combos{comboN,1}==areas{areaN}) %which chan
                if areaN<=4
                    areaF=areaN;
                else, areaP=areaN;
                end
            end
            if endsWith(combos{comboN,2},areas{areaN}) %which chan
                if areaN<=4
                    areaF=areaN;
                else, areaP=areaN;
                end
            end
        end
    end
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


% save delay, obj, recording_info & trial_info vars in
% m2-day1_delay_obj.mat
            
%find all unique values in a vector
unique(trial_info.CueObj)
                