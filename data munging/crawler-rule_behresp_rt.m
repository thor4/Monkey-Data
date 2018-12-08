%extract all behavior-related data for modeling rule switching
clear
monkeys = string({'betty', 'clark'});
days_betty = string({'090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001'});
days_clark = string({'060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221'});
betty = { days_betty, string('session01') };
clark = { days_clark, string('session02'), string('session03') };
% total number of trials indexed by monkey, day, session
delay = zeros(61400,4);
path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\trial_info.mat';
record_path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\recording_info.mat';

behavior=[]; %init var to hold 3 behavior var's for all good trials across all sessions
eyesr=[]; %init var to hold all eye sampling rates per day

for i=1:length(betty{1})
    myfilename = sprintf(path, monkeys(1), betty{1}{i}, betty{2}); %create full path to trial_info.mat
    load(myfilename); %load trial_info for day's trials
    myfilename_record = sprintf(record_path, monkeys(1), betty{1}{i}, betty{2}); %create full path to recording_info.mat
    load(myfilename_record); %load recording_info for day's trials
    eyesr = [eyesr recording_info.eye_sampleingrate];
    for j=1:trial_info.numTrials
        if (trial_info.good_trials(j) == 1) && ...%no artifacts
                ((trial_info.BehResp(j) == 0) || ... %only incorrect resp
                (trial_info.BehResp(j) == 1)) %or correct responses
            % extract rule, response and rt in 1x3 matrix
            behavior = [behavior; 
                trial_info.rule(j) trial_info.BehResp(j) trial_info.FirstSac(j)];
        end
    end
end
% 
% % clark loop
% for i=1:length(clark{1})
%     for j=2:3
%         myfilename = sprintf(path, monkeys(2), clark{1}{i}, clark{j});
%         load(myfilename);
%         delay(2,i,j) = trial_info.numTrials;
%     end
% end;
