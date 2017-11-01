%correct vs incorrect ERP
corAvg = mean(correct');
incorrAvg = mean(incorrect');
ncorAvg = corAvg + (-corAvg(1));
nincorrAvg = incorrAvg + (-incorrAvg(1));
time = 1:size(ncorAvg,2);
plot(time,ncorAvg,time,nincorrAvg,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Time-series');
xlabel('Time (ms)');
ylabel('Normalized Voltage (�V)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);

%correct frontal vs incorrect frontal & correct parietal vs incorrect parietal ERP
%frontal
corfAvg = mean(correctFrontal');
incorrfAvg = mean(incorrectFrontal');
ncorfAvg = corfAvg + (-corfAvg(1));
nincorrfAvg = incorrfAvg + (-incorrfAvg(1));
%parietal
corpAvg = mean(correctParietal');
incorrpAvg = mean(incorrectParietal');
ncorpAvg = corpAvg - (corpAvg(1));
nincorrpAvg = incorrpAvg + (-incorrpAvg(1));
%figure
time = 1:size(corfAvg,2);
figure
subplot(2,1,1)
plot(time,ncorfAvg,time,nincorrfAvg,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Frontal Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized Voltage (�V)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(2,1,2)
plot(time,ncorpAvg,time,nincorrpAvg,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Parietal Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized Voltage (�V)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);

%correct vs incorrect dPFC vs PEC ERP
%dPFC
corfAvg = mean(filtCFdpfc');
incorrfAvg = mean(filtIFdpfc');
ncorfAvg = corfAvg + (-corfAvg(1));
nincorrfAvg = incorrfAvg + (-incorrfAvg(1));
%PEC
corpAvg = mean(filtCPpec');
incorrpAvg = mean(filtIPpec');
ncorpAvg = corpAvg - (corpAvg(1));
nincorrpAvg = incorrpAvg + (-incorrpAvg(1));
%figure
time = 1:size(corfAvg,2);
figure
subplot(2,1,1)
plot(time,ncorfAvg,time,nincorrfAvg,':', 'LineWidth', 2);
title('Trial-averaged Delay-period dPFC Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized Voltage (�V)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(2,1,2)
plot(time,ncorpAvg,time,nincorrpAvg,':', 'LineWidth', 2);
title('Trial-averaged Delay-period PEC Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized Voltage (�V)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);

%total, frontal and parietal slices 
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