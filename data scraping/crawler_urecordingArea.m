clear
monkeys = string({'betty', 'clark'});
days_betty = string({'090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001'});
days_clark = string({'060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221'});
betty = { days_betty, string('session01') };
clark = { days_clark, string('session02'), string('session03') };
% total number of trials indexed by monkey, day, session
recordingArea = cell(2,24,3);
path = 'E:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\%s\\%s\\%s\\recording_info.mat';


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
    recordingArea(1,i,1) = {recording_info.area};
end;

% clark loop
for i=1:length(clark{1})
    for j=2:3
        myfilename = sprintf(path, monkeys(2), clark{1}{i}, clark{j});
        load(myfilename);
        recordingArea(2,i,j) = {recording_info.area};
    end
end;


% remove and collapse all nonzero cells
% the first 24 are betty days, next 23 are session02 from clark and last 23
% are session03 from clark. m is 70x1 cell
m = recordingArea(~cellfun(@isempty,recordingArea));
% put each session's recording area on a separate row and save only the
% unique areas recorded from (urecordingArea) row1 is session01, row2
% session02 and row3 session03
n = cell(1,24);
urecordingArea = cell(3,24);
k = 1;
for i=1:24
    for j=1:length(m{i})
        n{1,k} = char(m{i}{j});
        k = k + 1;
    end
end
z = unique(n);
for i=1:length(z)
    urecordingArea{1,i} = z(i);
end
n = cell(1,24);
k = 1;
for i=25:47
    for j=1:length(m{i})
        n{1,k} = char(m{i}{j});
        k = k + 1;
    end
end
z = unique(n);
for i=1:length(z)
    urecordingArea{2,i} = z(i);
end
n = cell(1,24);
k = 1;
for i=48:70
    for j=1:length(m{i})
        n{1,k} = char(m{i}{j});
        k = k + 1;
    end
end
z = unique(n);
for i=1:length(z)
    urecordingArea{3,i} = z(i);
end
