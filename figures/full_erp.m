load('good_trials-betty.mat')
load('good_trials-clark.mat')
load('bettyGoodRule1-split_by_BehResp_and_Region.mat')
load('clarkGoodRule1-split_by_BehResp_and_Region.mat')

Triggers = [0 500]; 
%BehResp x Region per Monkey ERP
%clark (Monkey 1)
correctFrontalClarkAvg = mean(correctFrontalClark);
incorrectFrontalClarkAvg = mean(incorrectFrontalClark);
correctParietalClarkAvg = mean(correctParietalClark);
incorrectParietalClarkAvg = mean(incorrectParietalClark);
%shift all traces to start at 0
correctFrontalClarkAvgShift = correctFrontalClarkAvg + (-correctFrontalClarkAvg(1));
incorrectFrontalClarkAvgShift = incorrectFrontalClarkAvg + (-incorrectFrontalClarkAvg(1));
correctParietalClarkAvgShift = correctParietalClarkAvg - (correctParietalClarkAvg(1));
incorrectParietalClarkAvgShift = incorrectParietalClarkAvg + (-incorrectParietalClarkAvg(1));
%betty (Monkey 2)
correctFrontalBettyAvg = mean(correctFrontalBetty);
incorrectFrontalBettyAvg = mean(incorrectFrontalBetty);
correctParietalBettyAvg = mean(correctParietalBetty);
incorrectParietalBettyAvg = mean(incorrectParietalBetty);
%shift all traces to start at 0
correctFrontalBettyAvgShift = correctFrontalBettyAvg + (-correctFrontalBettyAvg(1));
incorrectFrontalBettyAvgShift = incorrectFrontalBettyAvg + (-incorrectFrontalBettyAvg(1));
correctParietalBettyAvgShift = correctParietalBettyAvg - (correctParietalBettyAvg(1));
incorrectParietalBettyAvgShift = incorrectParietalBettyAvg + (-incorrectParietalBettyAvg(1));
%figure
time  = -500:size(correctFrontalBettyAvg,2)-501; % time, from -500ms baseline
figure
subplot(2,2,1)
traces1 = plot(time,correctFrontalClarkAvgShift,time,incorrectFrontalClarkAvgShift,':', 'LineWidth', 2);
title('Trial-averaged Monkey 1 Frontal Region Time-series'); xlabel('Time (ms)'); ylabel('Voltage (µV)'); 
set(gca,'box','off','Xlim',[time(1);time(end)]);
y1 = get(gca,'ylim'); hold on
triggers1 = plot([Triggers(1) Triggers(1)],y1,'--',[Triggers(2) Triggers(2)],y1,'--'); hold off 
legend(traces1,'Correct Trials','Incorrect Trials');
subplot(2,2,3)
traces2 = plot(time,correctParietalClarkAvgShift,time,incorrectParietalClarkAvgShift,':', 'LineWidth', 2);
title('Trial-averaged Monkey 1 Parietal Region Time-series'); xlabel('Time (ms)'); ylabel('Voltage (µV)');
set(gca,'box','off','Xlim',[time(1);time(end)]);
y2 = get(gca,'ylim'); hold on
triggers2 = plot([Triggers(1) Triggers(1)],y2,'--',[Triggers(2) Triggers(2)],y2,'--'); hold off
legend(traces2,'Correct Trials','Incorrect Trials');
subplot(2,2,2)
traces3 = plot(time,correctFrontalBettyAvgShift,time,incorrectFrontalBettyAvgShift,':', 'LineWidth', 2);
title('Trial-averaged Monkey 2 Frontal Region Time-series'); xlabel('Time (ms)'); ylabel('Voltage (µV)'); 
set(gca,'box','off','Xlim',[time(1);time(end)]);
y1 = get(gca,'ylim'); hold on
triggers3 = plot([Triggers(1) Triggers(1)],y1,'--',[Triggers(2) Triggers(2)],y1,'--'); hold off 
legend(traces3,'Correct Trials','Incorrect Trials');
subplot(2,2,4)
traces4 = plot(time,correctParietalBettyAvgShift,time,incorrectParietalBettyAvgShift,':', 'LineWidth', 2);
title('Trial-averaged Monkey 2 Parietal Region Time-series'); xlabel('Time (ms)'); ylabel('Voltage (µV)');
set(gca,'box','off','Xlim',[time(1);time(end)]);
y2 = get(gca,'ylim'); hold on
triggers4 = plot([Triggers(1) Triggers(1)],y2,'--',[Triggers(2) Triggers(2)],y2,'--'); hold off
legend(traces4,'Correct Trials','Incorrect Trials');