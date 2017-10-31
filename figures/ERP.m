load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_6DR_delay_cor_rule1.mat')
load('D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\betty\delay\betty_6DR_delay_inc_rule1.mat')
%create x variable for gramm (time)
%time = (1:size(b6DRdelayCorR1,1))';
% %create subset of data for gramm (behavioral response)
% cor = ones(size(b6DRdelayCorR1,2),1);
% inc = zeros(size(b6DRdelayIncR1,2),1);
% resp = vertcat(cor,inc);
%create gramm
% g=gramm('x',time,'y',volt,'color',resp);

%boundedline plot
%x axis
time = (1:size(cdPFCnormCorR1,1))';
%need average and standard deviation for ERP
avgC = mean(cdPFCnormCorR1,2);
avgI = mean(cdPFCnormIncR1,2);
stdC = std(cdPFCnormCorR1,0,2);
stdI = std(cdPFCnormIncR1,0,2);



%butterworth filter all data between 0-20Hz
filtOrder = 3; % filter order, higher values = steeper roll off
x = 0; %lower bound
y = 30; %upper bound
sr = 1000; %sampling rate
[b,a] = butter(filtOrder,20/(sr/2),'low'); %low pass, get coefficients a,b
% filter your data using the coefficients you computed above and save the 
% output. data should be an #samples x #trials matrix
cdPFCnormfiltCorR1 = filtfilt(b,a,cdPFCnormCorR1);
cdPFCnormfiltIncR1 = filtfilt(b,a,cdPFCnormIncR1);

avgCfilt = mean(cdPFCnormfiltCorR1,2);
avgIfilt = mean(cdPFCnormfiltIncR1,2);
stdCfilt = std(cdPFCnormfiltCorR1,0,2);
stdIfilt = std(cdPFCnormfiltIncR1,0,2);

figure
subplot(1,2,1)
title('Clark dPFC Rule 1 ERP 0-500Hz');
xlabel('Time (ms)');
ylabel('Normalized Voltage (µV)');
boundedline(time, avgC, stdC, time, avgI, stdI, 'r', 'alpha')
axis tight;
subplot(1,2,2)
xlabel('Time (ms)');
ylabel('Normalized Voltage (µV)');
title('Clark dPFC Rule 1 ERP 0-30Hz Butterworth Filter');
boundedline(time, avgCfilt, stdCfilt, time, avgIfilt, stdIfilt, 'r', 'alpha')
axis tight;