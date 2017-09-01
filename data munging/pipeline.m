%%load data from areas and rename
%correct
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/area6DR.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/area8AD.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/area8B.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/area9L.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/areadPFC.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/areaLIP.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/areaPE.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/areaPEC.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/areaPG.mat')
C6dr = area6DR; 
C8ad = area8AD; 
C8B = area8B; 
C9L = area9L; 
Cdpfc = areadPFC; 
Clip = areaLIP;
%Cmip = areaMIP; %not enough total trials (<2,000)
Cpe = areaPE;
Cpec = areaPEC; 
Cpg = areaPG;
%Cvpfc = areavPFC; %not enough total trials (<2,000)
clearvars area6DR area8AD area8B areadPFC areaLIP areaPE areaPEC areaPG
%incorrect
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/area6DR.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/area8AD.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/area8B.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/area9L.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/areadPFC.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/areaLIP.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/areaPE.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/areaPEC.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/areaPG.mat')
I6dr = area6DR;
I8ad = area8AD;
I8B = area8B;
I9l = area9L; 
Idpfc = areadPFC; 
Ilip = areaLIP;
%Imip = areaMIP; %not enough total trials (2,000)
Ipe = areaPE;
Ipec = areaPEC; 
Ipg = areaPG;
%Ivpfc = areavPFC; %not enough total trial (<2,000)
clearvars area6DR area8AD area8B areadPFC areaLIP areaPE areaPEC areaPG

%baseline correct
%correct
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/area6DR.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/area8AD.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/area8B.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/area9L.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/areadPFC.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/areaLIP.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/areaPE.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/areaPEC.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Correct-areas/baseline/areaPG.mat')
CB6dr = area6DR; 
CB8ad = area8AD; 
CB8B = area8B; 
CB9l = area9L;
CBdpfc = areadPFC; 
CBlip = areaLIP;
CBpe = areaPE;
CBpec = areaPEC; 
CBpg = areaPG;
clearvars area6DR area8AD area8B areadPFC areaLIP areaPE areaPEC areaPG
%incorrect
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/area6DR.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/area8AD.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/area8B.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/area9L.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/areadPFC.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/areaLIP.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/areaPE.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/areaPEC.mat')
load('/mnt/ceph/home/bconkli4/Documents/data/Incorrect-areas/baseline/areaPG.mat')
IB6dr = area6DR; 
IB8ad = area8AD; 
IB8B = area8B; 
IB9l = area9L;
IBdpfc = areadPFC; 
IBlip = areaLIP;
IBpe = areaPE;
IBpec = areaPEC; 
IBpg = areaPG;
clearvars area6DR area8AD area8B areadPFC areaLIP areaPE areaPEC areaPG

%calc single avg baseline signal across trials then across timepoints and
%subtract from delay period

CNorm6dr = C6dr - mean(mean(CB6dr'));
CNorm8ad = C8ad - mean(mean(CB8ad'));
CNorm8b = C8B - mean(mean(CB8B'));
CNorm9l = C9L - mean(mean(CB9l'));
CNormdpfc = Cdpfc - mean(mean(CBdpfc')); 
CNormlip = Clip - mean(mean(CBlip'));
CNormpe = Cpe - mean(mean(CBpe'));
CNormpec = Cpec - mean(mean(CBpec'));
CNormpg = Cpg - mean(mean(CBpg'));
INorm6dr = I6dr - mean(mean(IB6dr'));
INorm8ad = I8ad - mean(mean(IB8ad'));
INorm8b = I8B - mean(mean(IB8B'));
INorm9l = I9l - mean(mean(IB9l'));
INormdpfc = Idpfc - mean(mean(IBdpfc')); 
INormlip = Ilip - mean(mean(IBlip'));
INormpe = Ipe - mean(mean(IBpe'));
INormpec = Ipec - mean(mean(IBpec'));
INormpg = Ipg - mean(mean(IBpg'));
clearvars C6dr C8ad C8B C9l Cdpfc Clip Cpe Cpec Cpg
clearvars I6dr I8ad I8B I9l Idpfc Ilip Ipe Ipec Ipg
clearvars CB6dr CB8ad CB8B CB9l CBdpfc CBlip CBpe CBpec CBpg
clearvars IB6dr IB8ad IB8B IB9l IBdpfc IBlip IBpe IBpec IBpg 

%butterworth filter all data between 0-20Hz
filtOrder = 3; % filter order, higher values = steeper roll off
x = 0; %lower bound
y = 20; %upper bound
sr = 1000; %sampling rate
[b,a] = butter(filtOrder,20/(sr/2),'low'); %low pass, get coefficients a,b
% filter your data using the coefficients you computed above and save the 
% output. data should be an #samples x #trials matrix
CNormFilt6dr = filtfilt(b,a,CNorm6dr);
CNormFilt8ad = filtfilt(b,a,CNorm8ad);
CNormFilt8b = filtfilt(b,a,CNorm8b);
CNormFilt9l = filtfilt(b,a,CNorm9l);
CNormFiltdpfc = filtfilt(b,a,CNormdpfc);
CNormFiltlip = filtfilt(b,a,CNormlip);
CNormFiltpe = filtfilt(b,a,CNormpe);
CNormFiltpec = filtfilt(b,a,CNormpec);
CNormFiltpg = filtfilt(b,a,CNormpg);
INormFilt6dr = filtfilt(b,a,INorm6dr);
INormFilt8ad = filtfilt(b,a,INorm8ad);
INormFilt8b = filtfilt(b,a,INorm8b);
INormFilt9l = filtfilt(b,a,INorm9l);
INormFiltdpfc = filtfilt(b,a,INormdpfc);
INormFiltlip = filtfilt(b,a,INormlip);
INormFiltpe = filtfilt(b,a,INormpe);
INormFiltpec = filtfilt(b,a,INormpec);
INormFiltpg = filtfilt(b,a,INormpg);
clearvars CNorm6dr CNorm8ad CNorm8b CNorm9l CNormdpfc CNormlip CNormpe CNormpec CNormpg
clearvars INorm6dr INorm8ad INorm8b INorm9l INormdpfc INormlip INormpe INormpec INormpg 


%plot to see differences between correct and incorrect
CNormFiltAvg6dr = mean(CNormFilt6dr'); INormFiltAvg6dr = mean(INormFilt6dr');
CNormFiltAvg8ad = mean(CNormFilt8ad'); INormFiltAvg8ad = mean(INormFilt8ad');
CNormFiltAvg8b = mean(CNormFilt8b'); INormFiltAvg8b = mean(INormFilt8b');
CNormFiltAvg9l = mean(CNormFilt9l'); INormFiltAvg9l = mean(INormFilt9l');
CNormFiltAvgdpfc = mean(CNormFiltdpfc'); INormFiltAvgdpfc = mean(INormFiltdpfc');
CNormFiltAvglip = mean(CNormFiltlip'); INormFiltAvglip = mean(INormFiltlip');
CNormFiltAvgpe = mean(CNormFiltpe'); INormFiltAvgpe = mean(INormFiltpe');
CNormFiltAvgpec = mean(CNormFiltpec'); INormFiltAvgpec = mean(INormFiltpec');
CNormFiltAvgpg = mean(CNormFiltpg'); INormFiltAvgpg = mean(INormFiltpg');

%start them at zero
CNormFiltAvg6dr = CNormFiltAvg6dr + (-CNormFiltAvg6dr(1)); INormFiltAvg6dr = INormFiltAvg6dr + (-INormFiltAvg6dr(1));
CNormFiltAvg8ad = CNormFiltAvg8ad + (-CNormFiltAvg8ad(1)); INormFiltAvg8ad = INormFiltAvg8ad + (-INormFiltAvg8ad(1));
CNormFiltAvg8b = CNormFiltAvg8b + (-CNormFiltAvg8b(1)); INormFiltAvg8b = INormFiltAvg8b + (-INormFiltAvg8b(1));
CNormFiltAvg9l = CNormFiltAvg9l + (-CNormFiltAvg9l(1)); INormFiltAvg9l = INormFiltAvg9l + (-INormFiltAvg9l(1));
CNormFiltAvgdpfc = CNormFiltAvgdpfc + (-CNormFiltAvgdpfc(1)); INormFiltAvgdpfc = INormFiltAvgdpfc + (-INormFiltAvgdpfc(1));
CNormFiltAvglip = CNormFiltAvglip + (-CNormFiltAvglip(1)); INormFiltAvglip = INormFiltAvglip + (-INormFiltAvglip(1));
CNormFiltAvgpe = CNormFiltAvgpe + (-CNormFiltAvgpe(1)); INormFiltAvgpe = INormFiltAvgpe + (-INormFiltAvgpe(1));
CNormFiltAvgpec = CNormFiltAvgpec + (-CNormFiltAvgpec(1)); INormFiltAvgpec = INormFiltAvgpec + (-INormFiltAvgpec(1));
CNormFiltAvgpg = CNormFiltAvgpg + (-CNormFiltAvgpg(1)); INormFiltAvgpg = INormFiltAvgpg + (-INormFiltAvgpg(1));
%figure
time = 1:size(CNormFiltAvg6dr,2);
figure
subplot(5,2,1)
text(0.02,0.98,'a','Units', 'Normalized', 'VerticalAlignment', 'Top')
plot(time,mean(cNorm'),time,mean(iNorm'),':', 'LineWidth', 2);
title('Trial-averaged Delay-period All Areas');
xlabel('Time (ms)');
ylabel('Normalized Voltage (µV)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(5,2,2)
text(0.02,0.98,'b','Units', 'Normalized', 'VerticalAlignment', 'Top')
plot(time,CNormFiltAvg6dr,time,INormFiltAvg6dr,':', 'LineWidth', 2);
title('Trial-averaged Delay-period 6DR Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized µV');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(5,2,3)
text(0.02,0.98,'c','Units', 'Normalized', 'VerticalAlignment', 'Top')
plot(time,CNormFiltAvg8ad,time,INormFiltAvg8ad,':', 'LineWidth', 2);
title('Trial-averaged Delay-period 8AD Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized µV');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(5,2,4)
text(0.02,0.98,'d','Units', 'Normalized', 'VerticalAlignment', 'Top')
plot(time,CNormFiltAvg8b,time,INormFiltAvg8b,':', 'LineWidth', 2);
title('Trial-averaged Delay-period 8B Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized µV');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(5,2,5)
text(0.02,0.98,'e','Units', 'Normalized', 'VerticalAlignment', 'Top')
plot(time,CNormFiltAvg9l,time,INormFiltAvg9l,':', 'LineWidth', 2);
title('Trial-averaged Delay-period 9L Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized µV');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(5,2,6)
text(0.02,0.98,'f','Units', 'Normalized', 'VerticalAlignment', 'Top')
plot(time,CNormFiltAvgdpfc,time,INormFiltAvgdpfc,':', 'LineWidth', 2);
title('Trial-averaged Delay-period dPFC Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized µV');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(5,2,7)
text(0.02,0.98,'g','Units', 'Normalized', 'VerticalAlignment', 'Top')
plot(time,CNormFiltAvglip,time,INormFiltAvglip,':', 'LineWidth', 2);
title('Trial-averaged Delay-period LIP Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized µV');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(5,2,8)
text(0.02,0.98,'h','Units', 'Normalized', 'VerticalAlignment', 'Top')
plot(time,CNormFiltAvgpe,time,INormFiltAvgpe,':', 'LineWidth', 2);
title('Trial-averaged Delay-period PE Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized µV');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(5,2,9)
text(0.02,0.98,'i','Units', 'Normalized', 'VerticalAlignment', 'Top')
plot(time,CNormFiltAvgpec,time,INormFiltAvgpec,':', 'LineWidth', 2);
title('Trial-averaged Delay-period PEC Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized µV');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(5,2,10)
text(0.02,0.98,'j','Units', 'Normalized', 'VerticalAlignment', 'Top')
plot(time,CNormFiltAvgpg,time,INormFiltAvgpg,':', 'LineWidth', 2);
title('Trial-averaged Delay-period PG Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized µV');
%legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);

