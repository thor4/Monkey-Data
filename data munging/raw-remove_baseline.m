%calc single avg baseline signal across trials then across timepoints
basec = mean(mean(filtCb'));
basei = mean(mean(filtIb'));

%average signal across trials
noisec = mean(filtC');
noisei = mean(filtI');
%average signal across trials which were already normalized, now skip to 
%shift down to start at 0
cnaSub = mean(cSub);
inaSub = mean(iSub);

%remove baseline signal from delay
cNorm = filtC - basec;
iNorm = filtI - basei;

%shift down to start at 0
nc = cnaSub + (-cnaSub(1));
ni = inaSub + (-inaSub(1));
time = 1:size(nc,2);
plot(time,nc,time,ni,':', 'LineWidth', 2);
title('Trial-averaged Delay-period Time-series');
xlabel('Time (ms)');
ylabel('Normalized Voltage (�V)');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[time(1);time(end)]);

%reduce dataset
%split into frontal and parietal
cfNorm = cNorm(:, 1:162922);
cpNorm = cNorm(:, 162923:352967);
ifNorm = iNorm(:, 1:64786);
ipNorm = iNorm(:, 64787:137078);

%randomly sample to lower corpus count, produces trials x samples
cfnSub = datasample(cfNorm', 5000);
cpnSub = datasample(cpNorm', 5000);
ifnSub = datasample(ifNorm', 5000);
ipnSub = datasample(ipNorm', 5000);

%concatenate to create correct and incorrect subsamples
cSub = vertcat(cfnSub, cpnSub);
iSub = vertcat(ifnSub, ipnSub);