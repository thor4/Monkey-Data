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
trial_path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\trial_info.mat';
recording_path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\recording_info.mat';
trial_lfp_path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\%s%s%s.%04d.mat';
areas = ["9L", "8B", "6DR", "8AD", "vPFC", "dPFC", "LIP", "MIP", "PE", "PG", "PEC"];

tic
%betty loop
for i=1:length(betty{1})
    trial_myfilename = sprintf(trial_path, monkeys(1), betty{1}{i}, betty{2});
    load(trial_myfilename);
    recording_myfilename = sprintf(recording_path, monkeys(1), betty{1}{i}, betty{2});
    load(recording_myfilename);
    for j=1:trial_info.numTrials
        trial_lfp_myfilename = sprintf(trial_lfp_path, monkeys(1), betty{1}{i}, betty{2}, monkeys(1), betty{1}{i}, betty{2}{1}(8:9), j);
        load(trial_lfp_myfilename);
        [r,c] = size(lfp_data);
        lfp_data_mv = lfp_data .* 1000000; % convert to µV (1V = 10^6µV = 1,000,000µV)
        data(1:r,1:c,j) = lfp_data_mv; %these columns don't dynamically update, they just grow so there are extra 0's. can't remove all 0's since some are signal
    end
    day = sprintf('d%s',days_betty(i));
    bettyData.(day) = data;
    clear data
end
toc

% data = data(:,1:23349);


% clark loop
for i=1:length(clark{1})
    for j=2:3
        trial_myfilename = sprintf(trial_path, monkeys(2), clark{1}{i}, clark{j});
        load(trial_myfilename);
        recording_myfilename = sprintf(recording_path, monkeys(2), clark{1}{i}, clark{j});
        load(recording_myfilename);
        for k=1:trial_info.numTrials
            % parse through each good, stable trial with a certain response
            if (trial_info.good_trials(k) == 1) && (trial_info.stable_trials(k) == 1) && (trial_info.BehResp(k) == response) && (trial_info.rule(k) == rule)
                trial_lfp_myfilename = sprintf(trial_lfp_path,  monkeys(2), clark{1}{i}, clark{j}, monkeys(2), clark{1}{i}, clark{j}{1}(8:9), k);
                load(trial_lfp_myfilename);
                [r,c] = size(lfp_data);
                lfp_data_mv = lfp_data .* 1000000; % convert to µV (1V = 10^6µV = 1,000,000µV)
                delay_period = trial_info.MatchOnset(k)-trial_info.CueOffset(k);
                %baseline_period = (trial_info.CueOnset(k) - 50) - (trial_info.CueOnset(k) - 450);
%                 if (baseline_period ~= 400)
%                     bad_trials(idx,1) = monkeys(2);
%                     bad_trials(idx,2) = clark{1}{i};
%                     bad_trials(idx,3) = clark{j};
%                     bad_trials(idx,4) = k;
%                     idx = idx + 1;
%                 else
                    for l=1:r % parse through each channel
                    % insert voltage values during delay period, take one
                    % second after, up to 810ms after
                    % Pull out electrode in a specific recording area
                        if (recording_info.area(l) == area)
                            data(1:810,idx) = lfp_data_mv(l,floor(trial_info.CueOffset(k))+1:floor(trial_info.CueOffset(k))+810); %delay
                            %data(1:400,idx) = lfp_data_mv(l,floor(trial_info.CueOnset(k))-450:floor(trial_info.CueOnset(k))-51); %baseline
                            idx = idx + 1;
                        end
                    end
%                 end
            end
        end
    end
end
toc