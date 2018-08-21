i=1; %day 1
j=1; %chan 1
k=1; %trial 1
chan = fieldnames(monkey(1).day(i).correct);
ppc = monkey(1).day(1).correct.chan1PEC;
pfc = monkey(1).day(1).correct.chan3dPFC;

%day
d = 'd%i';
%correct
for i=1:numel(monkey(1).day)
    chan = fieldnames(monkey(1).day(i).correct);
    day = sprintf(d, i);
    for j=1:numel(chan)
        monkeyData = monkey(1).day(i).correct.(chan{j});
        monkeyDataTrial = monkeyData';
        monkeyDataTrial = reshape(monkeyDataTrial,1,1811,364);
        %transform into data structure like mikeXcohen
        chanXtimeXtrial(j,:,:) = monkeyDataTrial;
    end
end

%take only first trial
day1trial1 = chanXtimeXtrial(:,:,1);

%split into baseline, stimulus, delay1,2,3 epochs
%delay is split into 3 equal-sized 269ms windows: early, mid and late
[monkeyCorrBase, monkeyCovBase] = monkeyCov(day1trial1(:,1:500));
[monkeyCorrStimulus, monkeyCovStimulus] = monkeyCov(day1trial1(:,501:1001));
[monkeyCorrDelay1, monkeyCovDelay1] = monkeyCov(day1trial1(:,1002:1271));
[monkeyCorrDelay2, monkeyCovDelay2] = monkeyCov(day1trial1(:,1272:1541));
[monkeyCorrDelay3, monkeyCovDelay3] = monkeyCov(day1trial1(:,1542:1811));

%variances
diag(monkeyCovStimulus)


monkeyDataStimulus = bsxfun(@minus,day1trial1(:,501:1001),mean(day1trial1(:,501:1001),2));
monkeyCovStimulus = monkeyDataStimulus*monkeyDataStimulus'/(size(day1trial1(:,501:1001),2)-1);
imagesc(monkeyCovStimulus)
colorbar


histogram(pfc(:,1343))


rng('default')
x = randn(30,4);
y = randn(30,4);
y(:,4) = sum(x,2); % introduce correlation

[r,p] = corr(pfc,ppc);
%corr(X,Y) returns a p1-by-p2 matrix containing the pairwise correlation 
%coefficient between each pair of columns in the n-by-p1 and n-by-p2 
%matrices X and Y.

sig = (p < 0.001);

imagesc(r)
colorbar