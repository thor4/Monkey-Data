%calc single avg baseline signal across trials then across timepoints
basec = mean(mean(filtCb'));
basei = mean(mean(filtIb'));

%average signal across trials
noisec = mean(filtC');
noisei = mean(filtI');

%remove baseline signal from delay
c = filtC - basec;
i = filtI - basei;
cAvg = noisec - basec;
iAvg = noisei - basei;

%shift down to start at 0
nc = cAvg + (-cAvg(1));
ni = iAvg + (-iAvg(1));
time = 1:size(nc,2);
plot(time,nc,time,ni,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Time-series');
xlabel('Time (ms)');
ylabel('Normalized Voltage (µV)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);
