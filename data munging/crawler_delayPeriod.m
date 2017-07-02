clear
monkeys = string({'betty', 'clark'});
days_betty = string({'090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001'});
days_clark = string({'060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221'});
betty = { days_betty, string('session01') };
clark = { days_clark, string('session02'), string('session03') };
% initialize delay matrix: first column is day, second is session, third is
% trial number, fourth is CueOffset, fifth is MatchOnset and sixth is total
% delay period in ms
delay = zeros(61400,6);
days_betty_num = [90615, 90616, 90617, 90618, 90622, 90625, 90626, 90629, 90701, 90702, 90706, 90708, 90709, 90901, 90903, 90916, 90917, 90921, 90923, 90924, 90928, 90929, 90930, 91001];
days_clark_num = [60328, 60406, 60411, 60414, 60426, 60427, 60428, 60502, 60503, 60509, 60511, 60531, 60601, 60602, 60824, 60825, 60831, 60907, 61212, 61213, 61214, 61215, 61221];
trial_path = 'E:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\trial_info.mat';
recording_path = 'E:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\recording_info.mat';
trial_lfp_path = 'E:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\trial_info.mat';

% counter
idx = 1;

for i=1:length(betty{1})
    trial_myfilename = sprintf(trial_path, monkeys(1), betty{1}{i}, betty{2});
    load(trial_myfilename);
    recording_myfilename = sprintf(recording_path, monkeys(1), betty{1}{i}, betty{2});
    load(recording_myfilename);
    for j=1:trial_info.numTrials
        if (trial_info.good_trials(j) == 1) && (trial_info.stable_trials(j) == 1)
            trial_myfilename = sprintf(trial_path, monkeys(1), betty{1}{i}, betty{2});
            load(trial_myfilename);
            delay(idx,1:1222) = days_betty_num(i);
            delay(idx,2) = 1;
            delay(idx,3) = j;
            delay(idx,4) = trial_info.CueOffset(j);
            delay(idx,5) = trial_info.MatchOnset(j);
            delay(idx,6) = delay(idx,5) - delay(idx,4);
            idx = idx + 1;
        end
    end
end;

% clark loop
for i=1:length(clark{1})
    for j=2:3
        trial_myfilename = sprintf(trial_path, monkeys(2), clark{1}{i}, clark{j});
        load(trial_myfilename);
        for k=1:trial_info.numTrials
            if (trial_info.good_trials(k) == 1) && (trial_info.stable_trials(k) == 1)
                delay(idx,1) = days_clark_num(i);
                delay(idx,2) = j;
                delay(idx,3) = k;
                delay(idx,4) = trial_info.CueOffset(k);
                delay(idx,5) = trial_info.MatchOnset(k);
                delay(idx,6) = delay(idx,5) - delay(idx,4);
                idx = idx + 1;
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