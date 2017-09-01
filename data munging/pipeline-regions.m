%%load data from regions and rename
%correct
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/correctFrontal.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/correctParietal.mat')
Cf = correctFrontal; Cp = correctParietal;
clearvars correctFrontal correctParietal
%incorrect
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/incorrectFrontal.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/incorrectParietal.mat')
If = incorrectFrontal; Ip = incorrectParietal;
clearvars incorrectFrontal incorrectParietal

%baseline correct
%correct
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/correctFrontal.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/correctParietal.mat')
CBf = correctFrontal; CBp = correctParietal;
clearvars correctFrontal correctParietal
%incorrect
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/incorrectFrontal.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/incorrectParietal.mat')
IBf = incorrectFrontal; IBp = incorrectParietal;
clearvars incorrectFrontal incorrectParietal

%calc single avg baseline signal across trials then across timepoints and
%subtract from delay period
CNormf = Cf - mean(mean(CBf')); CNormp = Cp - mean(mean(CBp'));
INormf = If - mean(mean(IBf')); INormp = Ip - mean(mean(IBp'));
clearvars CBf CBp Cf Cp IBf IBp If Ip

%butterworth filter all data between 0-20Hz
filtOrder = 3; % filter order, higher values = steeper roll off
x = 0; %lower bound
y = 20; %upper bound
sr = 1000; %sampling rate
[b,a] = butter(filtOrder,20/(sr/2),'low'); %low pass, get coefficients a,b
% filter your data using the coefficients you computed above and save the 
% output. data should be an #samples x #trials matrix
CNormFiltf = filtfilt(b,a,CNormf); CNormFiltp = filtfilt(b,a,CNormp);
INormFiltf = filtfilt(b,a,INormf); INormFiltp = filtfilt(b,a,INormp);
clearvars CNormf CNormp INormf INormp

%concatenate across regions to create correct vs incorrect
CNormFilt = horzcat(CNormFiltf, CNormFiltp);
INormFilt = horzcat(INormFiltf, INormFiltp);

%plot to see differences between correct and incorrect
CNormFiltAvg = mean(CNormFilt'); INormFiltAvg = mean(INormFilt');
CNormFiltAvgf = mean(CNormFiltf'); INormFiltAvgf = mean(INormFiltf');
CNormFiltAvgp = mean(CNormFiltp'); INormFiltAvgp = mean(INormFiltp');

%start them at zero
CNormFiltAvg = CNormFiltAvg + (-CNormFiltAvg(1)); INormFiltAvg = INormFiltAvg + (-INormFiltAvg(1));
CNormFiltAvgf = CNormFiltAvgf + (-CNormFiltAvgf(1)); INormFiltAvgf = INormFiltAvgf + (-INormFiltAvgf(1));
CNormFiltAvgp = CNormFiltAvgp + (-CNormFiltAvgp(1)); INormFiltAvgp = INormFiltAvgp + (-INormFiltAvgp(1));

%figure
time = 1:size(CNormFiltAvg,2);
figure
subplot(2,2,[1 2])
plot(time,CNormFiltAvg,time,INormFiltAvg,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Both Regions');
xlabel('Time after stimulus off (ms)');
ylabel('Potential (µV)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(2,2,3)
plot(time,CNormFiltAvgf,time,INormFiltAvgf,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Frontal Region');
xlabel('Time after stimulus off (ms)');
ylabel('Potential (µV)');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(2,2,4)
plot(time,CNormFiltAvgp,time,INormFiltAvgp,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Parietal Region');
xlabel('Time after stimulus off (ms)');
ylabel('Potential (µV)');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
set(gcf,'color','white')