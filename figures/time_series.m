corAvg = mean(correct');
incorrAvg = mean(incorrect');
time = 1:size(corAvg,2);
plot(time,corAvg,time,incorrAvg,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Time-series');
xlabel('Time (ms)');
ylabel('Voltage (µV)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);