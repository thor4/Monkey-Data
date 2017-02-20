clear
monkeys = string({'betty', 'clark'});
days_betty = string({'090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001'});
days_clark = string({'060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221'});
betty = { days_betty, string('session01') };
clark = { days_clark, string('session02'), string('session03') };
% total number of trials indexed by monkey, day, session
delay = zeros(61400,4);
path = 'E:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\trial_info.mat';


% load each of the gsc trials and save to specific delay file
% for k = 1:gsc_tot
%     myfilename = sprintf('betty09061701.%04d.mat', gsc_trials(k));
%     load(myfilename);
%     lfp_data_delay=lfp_data(:,trial_info.CueOffset(gsc_trials(k)):(trial_info.MatchOnset(gsc_trials(k)))-1);
%     save(sprintf('betty09061701.d.%04d.mat',gsc_trials(k)), 'lfp_data_delay');
% end

for i=1:length(betty{1})
    myfilename = sprintf(path, monkeys(1), betty{1}{i}, betty{2});
    load(myfilename);
    numTrials = trial_info.numTrials;
    for j=1:trial_info.numTrials
        delay(j,1) = j;
        delay(j,2) = trial_info.CueOffset(j);
        delay(j,3) = trial_info.MatchOnset(j);
        delay(j,4) = delay(j,3) - delay(j,2);
    end
end;

% clark loop
for i=1:length(clark{1})
    for j=2:3
        myfilename = sprintf(path, monkeys(2), clark{1}{i}, clark{j});
        load(myfilename);
        delay(2,i,j) = trial_info.numTrials;
    end
end;

% collapse to single first dimension since it's redundant to have betty in
% first row and clark in second when they are already distinguished by the
% session (third dimension)
% a = sum(numTrials);
% find minimums per session by inflating the zeroes to 9999: 
% min(a+9999*(a==0))
% find maximums per session max(a)