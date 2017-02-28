clear
monkeys = string({'betty', 'clark'});
days_betty = string({'090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001'});
days_clark = string({'060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221'});
betty = { days_betty, string('session01') };
clark = { days_clark, string('session02'), string('session03') };
% initialize data matrix:
% lfpData-sessionID-day-trialNum-chanLoc-chan-rule-resp [double floating point array]
data = zeros(1000000,1222+3+1+1+2+11+2+2);
days_betty_num = [90615, 90616, 90617, 90618, 90622, 90625, 90626, 90629, 90701, 90702, 90706, 90708, 90709, 90901, 90903, 90916, 90917, 90921, 90923, 90924, 90928, 90929, 90930, 91001];
days_clark_num = [60328, 60406, 60411, 60414, 60426, 60427, 60428, 60502, 60503, 60509, 60511, 60531, 60601, 60602, 60824, 60825, 60831, 60907, 61212, 61213, 61214, 61215, 61221];
trial_path = 'E:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\trial_info.mat';
recording_path = 'E:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\recording_info.mat';
trial_lfp_path = 'E:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\%s%s%s.%04d.mat';

idx = 1; % counter

% betty loop
for i=1:length(betty{1})
    trial_myfilename = sprintf(trial_path, monkeys(1), betty{1}{i}, betty{2});
    load(trial_myfilename);
    recording_myfilename = sprintf(recording_path, monkeys(1), betty{1}{i}, betty{2});
    load(recording_myfilename);
    for j=1:trial_info.numTrials
        % parse through each good, stable trial with a delay period less
        % than 1222ms
        if (trial_info.good_trials(j) == 1) && (trial_info.stable_trials(j) == 1) && (trial_info.MatchOnset(j)-trial_info.CueOffset(j)<1223)
            trial_lfp_myfilename = sprintf(trial_lfp_path, monkeys(1), betty{1}{i}, betty{2}, monkeys(1), betty{1}{i}, betty{2}{1}(8:9), j);
            load(trial_lfp_myfilename);
            [r,c] = size(lfp_data);
            delay_period = trial_info.MatchOnset(j)-trial_info.CueOffset(j);
            for k=1:r % parse through each channel
                % insert voltage values during delay period, take one second after
				% sample goes away up through one second before match appears
                data(idx,1:delay_period-1) = lfp_data(k,trial_info.CueOffset(j)+1:trial_info.MatchOnset(j)-1);
                data(idx,1223:1225) = [1,0,0]; % sessionID 1 for betty
                data(idx,1226) = days_betty_num(i); % day
                data(idx,1227) = j; % trial number
                % channel location, the + changes from logical type to
                % double [F, P]
                data(idx,1228:1229) = +[recording_info.cortex(k)=='F',recording_info.cortex(k)=='P'];
                % channel [9L, 8B, 6DR, 8AD, vPFC, dPFC, LIP, MIP, PE, PG,
                % PEC]
                data(idx,1230:1240) = +[recording_info.area{k}==string('9L') ...
                    recording_info.area{k}==string('8B'), recording_info.area{k}==string('6DR') ...
                    recording_info.area{k}==string('8AD'), recording_info.area{k}==string('vPFC') ...
                    recording_info.area{k}==string('dPFC'), recording_info.area{k}==string('LIP') ...
                    recording_info.area{k}==string('MIP'), recording_info.area{k}==string('PE') ...
                    recording_info.area{k}==string('PG'), recording_info.area{k}==string('PEC')];
                % rule in play [identity, location]
                data(idx,1241:1242) = +[trial_info.rule(j)==1, trial_info.rule(j)==2];
                % correct or incorrect trial [correct, incorrect]
                data(idx,1243:1244) = +[trial_info.BehResp(j)==1, trial_info.BehResp(j)==0];
                idx = idx + 1;
            end 
        end
    end
end;

% clark loop
for i=1:length(clark{1})
    for j=2:3
        trial_myfilename = sprintf(trial_path, monkeys(2), clark{1}{i}, clark{j});
        load(trial_myfilename);
        recording_myfilename = sprintf(recording_path, monkeys(2), clark{1}{i}, clark{j});
        load(recording_myfilename);
        for k=1:trial_info.numTrials
            % parse through each good, stable trial with a delay period less
            % than 1222ms
            if (trial_info.good_trials(k) == 1) && (trial_info.stable_trials(k) == 1) && (trial_info.MatchOnset(k)-trial_info.CueOffset(k)<1223)
                trial_lfp_myfilename = sprintf(trial_lfp_path,  monkeys(2), clark{1}{i}, clark{j}, monkeys(2), clark{1}{i}, clark{j}{1}(8:9), k);
                load(trial_lfp_myfilename);
                [r,c] = size(lfp_data);
                delay_period = trial_info.MatchOnset(k)-trial_info.CueOffset(k);
                for l=1:r % parse through each channel
                    % insert voltage values during delay period, take one second after
					% sample goes away up through one second before match appears
                    data(idx,1:delay_period-1) = lfp_data(k,trial_info.CueOffset(k)+1:trial_info.MatchOnset(k)-1);
                    if j==2
                        data(idx,1223:1225) = [0,1,0]; % sessionID 2 for clark
                    else
                        data(idx,1223:1225) = [0,0,1]; % sessionID 3 for clark
                    end
                    data(idx,1226) = days_clark_num(i); % day
                    data(idx,1227) = k; % trial number
                    % channel location, the + changes from logical type to
                    % double [F, P]
                    data(idx,1228:1229) = +[recording_info.cortex(l)=='F',recording_info.cortex(l)=='P'];
                    % channel [9L, 8B, 6DR, 8AD, vPFC, dPFC, LIP, MIP, PE, PG,
                    % PEC]
                    data(idx,1230:1240) = +[recording_info.area{l}==string('9L') ...
                        recording_info.area{l}==string('8B'), recording_info.area{l}==string('6DR') ...
                        recording_info.area{l}==string('8AD'), recording_info.area{l}==string('vPFC') ...
                        recording_info.area{l}==string('dPFC'), recording_info.area{l}==string('LIP') ...
                        recording_info.area{l}==string('MIP'), recording_info.area{l}==string('PE') ...
                        recording_info.area{l}==string('PG'), recording_info.area{l}==string('PEC')];
                    % rule in play [identity, location]
                    data(idx,1241:1242) = +[trial_info.rule(k)==1, trial_info.rule(k)==2];
                    % correct or incorrect trial [correct, incorrect]
                    data(idx,1243:1244) = +[trial_info.BehResp(k)==1, trial_info.BehResp(k)==0];
                    idx = idx + 1;
                end 
            end
        end
    end
end;


% collapse to single first dimension since it's redundant to have betty in
% first row and clark in second when they are already distinguished by the
% session (third dimension)
% a = sum(numTrials);
% find minimums per session by inflating the zeroes to 9999: 
% min(a+9999*(a==0))
% find maximums per session max(a)