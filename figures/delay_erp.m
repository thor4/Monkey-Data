%correct vs incorrect ERP
corAvg = mean(correct');
incorrAvg = mean(incorrect');
ncorAvg = corAvg + (-corAvg(1));
nincorrAvg = incorrAvg + (-incorrAvg(1));
time = 1:size(ncorAvg,2);
plot(time,ncorAvg,time,nincorrAvg,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Time-series');
xlabel('Time (ms)');
ylabel('Normalized Voltage (µV)');
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
ylabel('Normalized Voltage (µV)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
subplot(2,1,2)
plot(time,ncorpAvg,time,nincorrpAvg,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Parietal Region Time-series');
xlabel('Time (ms)');
ylabel('Normalized Voltage (µV)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);


