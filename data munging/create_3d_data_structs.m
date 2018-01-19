clear
monkeys = ["betty", "clark"];
days_betty = ["090615", "090616", "090617", "090618", "090622", "090625", "090626", "090629", "090701", "090702", "090706", "090708", "090709", "090901", "090903", "090916", "090917", "090921", "090923", "090924", "090928", "090929", "090930", "091001"];
days_clark = ["060328", "060406", "060411", "060414", "060426", "060427", "060428", "060502", "060503", "060509", "060511", "060531", "060601", "060602", "060824", "060825", "060831", "060907", "061212", "061213", "061214", "061215", "061221"];
betty = { days_betty, "session01" };
clark = { days_clark, "session02", "session03" };
% initialize data matrix:
% data = zeros(25,7000,2000); %delay period total
days_betty_num = [90615, 90616, 90617, 90618, 90622, 90625, 90626, 90629, 90701, 90702, 90706, 90708, 90709, 90901, 90903, 90916, 90917, 90921, 90923, 90924, 90928, 90929, 90930, 91001];
days_clark_num = [60328, 60406, 60411, 60414, 60426, 60427, 60428, 60502, 60503, 60509, 60511, 60531, 60601, 60602, 60824, 60825, 60831, 60907, 61212, 61213, 61214, 61215, 61221];
trial_path = 'C:\\Users\\bryan\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\trial_info.mat';
recording_path = 'C:\\Users\\bryan\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\recording_info.mat';
trial_lfp_path = 'C:\\Users\\bryan\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\%s%s%s.%04d.mat';
% trial_path = '//home//bconkli4//Documents//data//original//%s//%s//%s//trial_info.mat';
% recording_path = '//home//bconkli4//Documents//data//original//%s//%s//%s//recording_info.mat';
% trial_lfp_path = '//home//bconkli4//Documents//data//original//%s//%s//%s//%s%s%s.%04d.mat';
areas = ["9L", "8B", "6DR", "8AD", "vPFC", "dPFC", "LIP", "MIP", "PE", "PG", "PEC"];
% trial_length = 2100; %in ms, found in short-trial search

tic
%betty loop
for i=1:length(betty{1})
    trial_myfilename = sprintf(trial_path, monkeys(1), betty{1}{i}, betty{2});
    load(trial_myfilename);
    recording_myfilename = sprintf(recording_path, monkeys(1), betty{1}{i}, betty{2});
    load(recording_myfilename);
    %exploring trial lengths    
    day = sprintf('d%s',days_betty(i));
    z = 1;
    for j=1:trial_info.numTrials
        if (trial_info.good_trials(j) == 1) && (trial_info.stable_trials(j) == 1)
              trial_lfp_myfilename = sprintf(trial_lfp_path, monkeys(1), betty{1}{i}, betty{2}, monkeys(1), betty{1}{i}, betty{2}{1}(8:9), j);
              load(trial_lfp_myfilename);
              trial = sprintf('t%d',j);
              [r,c] = size(lfp_data);
    %         bTrialLength(j,1) = c;
              lfp_data_mv = lfp_data .* 1000000; % convert to µV (1V = 10^6µV = 1,000,000µV)
    %         data(1:r,1:trial_length,j) = lfp_data_mv(:,1:trial_length); %these columns don't dynamically update, they just grow so there are extra 0's. can't remove all 0's since some are signal
%             bettyData.(day).(trial) = lfp_data_mv;
%             clear lfp_data_mv
              bettyGoodStableTrials.(day).(trial) = lfp_data_mv(:,trial_info.CueOnset(j)-500:end); %store shifted time series with 500ms prestimulus baseline as start
              new_trial_info.TrialLength(z) = size(bettyGoodStableTrials.(day).(trial),2);
              offset = c - new_trial_info.TrialLength(z); %calc offset for revised markers
              new_trial_info.TrialNum(z) = j;
              new_trial_info.CueOnset(z) = trial_info.CueOnset(j) - offset; %0 is now a point, so new CueOnset will be 501 for all shifted trials
              new_trial_info.CueOffset(z) = trial_info.CueOffset(j) - offset; %should be 1001-1051 for 500-550 ms sample window
              new_trial_info.MatchOnset(z) = trial_info.MatchOnset(j) - offset; 
              new_trial_info.BehResp(z) = trial_info.BehResp(j); new_trial_info.Rule(z) = trial_info.rule(j);
              new_trial_info.CueObj(z) = trial_info.CueObj(j); new_trial_info.CueLoc(z) = trial_info.CueLoc(j);
              new_trial_info.MatchObj1(z) = trial_info.MatchObj1(j); new_trial_info.MatchObj2(z) = trial_info.MatchObj2(j);
              new_trial_info.MatchPos1(z) = trial_info.MatchPos1(j); new_trial_info.MatchPos2(z) = trial_info.MatchPos2(j);
              new_trial_info.FirstSac(z) = trial_info.FirstSac(j);
              z = z + 1;
              clear lfp_data_mv
              %need to determine the distribution of sample-size windows, do a
              %histogram of (CueOffset - CueOnset) values and see for one
              %day. ask bressler how you do a spectrogram incl baseline and
              %sample in addition to delay if you have so many variable
              %window sizes, of which sample is just one
        end
    end
    bettyGoodStableTrials.(day).new_trial_info = new_trial_info;
    bettyGoodStableTrials.(day).recording_info = recording_info;
    %creating cell within a struct for each recording session's data
%     bettyTrials.(day) = trial_info.TrialLength;
%     bettyData.(day) = data;
%     clear data
%     trial_info.TrialLength = bTrialLength;
%     save(trial_myfilename, 'trial_info');
%     clear bTrialLength
    clear new_trial_info
    clear recording_info
end
toc
% data = data(:,1:23349);

tic
% clark loop
for i=1:length(clark{1})
    day = sprintf('d%s',days_clark(i));
    for j=2:3
        sesh = string(clark(j));
        z = 1;
        trial_myfilename = sprintf(trial_path, monkeys(2), clark{1}{i}, clark{j});
        load(trial_myfilename);
        recording_myfilename = sprintf(recording_path, monkeys(2), clark{1}{i}, clark{j});
        load(recording_myfilename);
        for k=1:trial_info.numTrials
            if (trial_info.good_trials(k) == 1) && (trial_info.stable_trials(k) == 1)
              trial_lfp_myfilename = sprintf(trial_lfp_path,  monkeys(2), clark{1}{i}, clark{j}, monkeys(2), clark{1}{i}, clark{j}{1}(8:9), k);
              load(trial_lfp_myfilename);
              trial = sprintf('t%d',k);
%             cTrialLength(k,1) = c;
              lfp_data_mv = lfp_data .* 1000000; % convert to µV (1V = 10^6µV = 1,000,000µV)
              [r,c] = size(lfp_data);
%             delay_period = trial_info.MatchOnset(k)-trial_info.CueOffset(k);
%                 for l=1:r % parse through each channel
%                 % insert voltage values during delay period, take one
%                 % second after, up to 810ms after
%                 % Pull out electrode in a specific recording area
%                     if (recording_info.area(l) == area)
%                         data(1:810,idx) = lfp_data_mv(l,floor(trial_info.CueOffset(k))+1:floor(trial_info.CueOffset(k))+810); %delay
%                         %data(1:400,idx) = lfp_data_mv(l,floor(trial_info.CueOnset(k))-450:floor(trial_info.CueOnset(k))-51); %baseline
%                         idx = idx + 1;
%                     end
%                 end
%                 end
              clarkGoodStableTrials.(day).(sesh).(trial) = lfp_data_mv(:,floor(trial_info.CueOnset(k)-500):end); %store shifted time series with 500ms prestimulus baseline as start
              new_trial_info.TrialLength(z) = size(clarkGoodStableTrials.(day).(sesh).(trial),2);
              offset = c - new_trial_info.TrialLength(z); %calc offset for revised markers
              new_trial_info.TrialNum(z) = k;
              new_trial_info.CueOnset(z) = floor(trial_info.CueOnset(k) - offset); %0 is now a point, so new CueOnset will be 501 for all shifted trials
              new_trial_info.CueOffset(z) = floor(trial_info.CueOffset(k) - offset); %should be 1001-1051 for 500-550 ms sample window
              new_trial_info.MatchOnset(z) = floor(trial_info.MatchOnset(k) - offset); 
              new_trial_info.BehResp(z) = trial_info.BehResp(k); new_trial_info.Rule(z) = trial_info.rule(k);
              new_trial_info.CueObj(z) = trial_info.CueObj(k); new_trial_info.CueLoc(z) = trial_info.CueLoc(k);
              new_trial_info.MatchObj1(z) = trial_info.MatchObj1(k); new_trial_info.MatchObj2(z) = trial_info.MatchObj2(k);
              new_trial_info.MatchPos1(z) = trial_info.MatchPos1(k); new_trial_info.MatchPos2(z) = trial_info.MatchPos2(k);
              new_trial_info.FirstSac(z) = trial_info.FirstSac(k);
              z = z + 1;
              clear lfp_data_mv
            end
        end
        clarkGoodStableTrials.(day).(sesh).new_trial_info = new_trial_info;
        clarkGoodStableTrials.(day).(sesh).recording_info = recording_info;
        clear new_trial_info
        clear recording_info
    end
end
toc

%determine ultimate shortest trial
fieldsB = fieldnames(bettyTrials);
fieldsC = fieldnames(clarkTrials);
for i = 1:numel(fieldsB)
  minBetty(i) = min(bettyTrials.(fieldsB{i}));
end
minBetty = minBetty';
for i = 1:numel(fieldsC)
  minClark(i) = min(clarkTrials.(fieldsC{i}));
end
minClark = minClark';
min(minBetty)
%2148: d090615 & d091001
min(minClark)
%2198: d061221session02

%extra to crawl trials
minBetty = zeros(numel(daysB),1);
for i = 1:numel(daysB)
    trialsB = fieldnames(bettyData.(daysB{i}));
    for j = 1:numel(trialsB)
        if j == 1
            minBetty(i) = size(bettyData.(daysB{i}).(trialsB{j}),2);
        elseif (size(bettyData.(daysB{i}).(trialsB{j}),2) < size(bettyData.(daysB{i}).(trialsB{j-1}),2))
            minBetty(i) = size(bettyData.(daysB{i}).(trialsB{j}),2);
        else
            minBetty(i) = size(bettyData.(daysB{i}).(trialsB{j-1}),2);
        end
    end
end
minBetty = minBetty';
for i = 1:numel(fieldsC)
  minClark(i) = min(clarkTrials.(fieldsC{i}));
end
