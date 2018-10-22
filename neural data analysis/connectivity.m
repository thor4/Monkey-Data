i=1; %day 1
j=1; %resp 1
k=1; %chan 1
l=1; %trial 1
%initialize all variables
resp = fieldnames(monkey(1).day(1));
time_window = 100; %ms
time_step = 5; %ms
srate = 1000; %hz
% convert ms to indices
% (1000/EEG.srate) provides number of ms per cycle(sample) ie: 1 cycle / 3.90625 ms
% timewindow/(1000/EEG.srate) provides number of cycles (samples) over
% entire timewindow of interest. divide by 2 because you want equal number
% of samples before and after time point of interest, which sum up to total
% number of samples in timewindow
timewindowidx = round(time_window/(1000/srate)/2);
monkeyTime = -500:1:1310; %total trial time (ms)
times2save = -450:time_step:1260; %allows room for half of time_window on front & back
%determine which time points are of interest for rolling 100ms time window
%with defined time_step
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(monkeyTime-times2save(i)));
end

%day
d = 'd%i';

%% Step 1: compute optimal number of bins for each variable
tic
%monkey 1
for i=1:numel(monkey(1).day)
    for j=1:numel(resp)
        chan = fieldnames(monkey(1).day(i).(resp{j}));
        % C = combnk(v,k) returns all combinations of the n elements in v taken k at a time.
        C = combnk(chan,2); %all chan combinations
        day = sprintf(d, i);
        for k=1:size(C,1) %total number of combinations
            for timei = 1:length(times2save)
                %pull out channel combination
                data1 = monkey(1).day(i).(resp{j}).(C{k,1});
                monkeyData1 = data1(:,times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx);
                data2 = monkey(1).day(i).(resp{j}).(C{k,2});
                monkeyData2 = data1(:,times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx);
                % compute optimal number of bins for each variable
                numbins(1,timei,k) = numBins(monkeyData1(:));
                numbins(2,timei,k) = numBins(monkeyData2(:));
            end
        end
        m1bins.(day).(resp{j}).numbins = numbins;
        clear numbins
    end
end

%monkey 2
for i=1:numel(monkey(2).day)
    for j=1:numel(resp)
        chan = fieldnames(monkey(2).day(i).(resp{j}));
        %C = combnk(v,k) returns all combinations of the n elements in v taken k at a time.
        C = combnk(chan,2); %all chan combinations
        day = sprintf(d, i);
        for k=1:size(C,1) %total number of combinations
            for timei = 1:length(times2save)
                %pull out channel combination
                data1 = monkey(2).day(i).(resp{j}).(C{k,1});
                monkeyData1 = data1(:,times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx);
                data2 = monkey(2).day(i).(resp{j}).(C{k,2});
                monkeyData2 = data1(:,times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx);
                % compute optimal number of bins for each variable
                numbins(1,timei,k) = numBins(monkeyData1(:));
                numbins(2,timei,k) = numBins(monkeyData2(:));
            end
        end
        m2bins.(day).(resp{j}).numbins = numbins;
        clear numbins 
    end
end
toc
