s = 1000;
[fftpower_cAvgDelay, cfaxd] = fftPow(s, filtC', 810);
fftpower_iAvgDelay = fftPow(s, filtI', 810);
[fftpower_cAvgBase, cfax] = fftPow(s, filtCb', 810);
fftpower_iAvgBase = fftPow(s, filtIb', 810);

%baseline normalization
fftpower_cAvgNorm = fftpower_cAvgDelay ./ fftpower_cAvgBase;
fftpower_iAvgNorm = fftpower_iAvgDelay ./ fftpower_iAvgBase;

%only interested in frequencies 1-250Hz
freqIdx = find(cfax>=0 & cfax<250);
freqfft = cfax(freqIdx(1):freqIdx(end));
fftpower_cAvgNorm = fftpower_cAvgNorm(freqIdx(1):freqIdx(end));
fftpower_iAvgNorm = fftpower_iAvgNorm(freqIdx(1):freqIdx(end));

%reduce dataset of filtered 0-20Hz fft power
%split into frontal and parietal
cfNorm = fftCfilt(1:162922,:);
cpNorm = fftCfilt(162923:352967,:);
ifNorm = fftIfilt(1:64786,:);
ipNorm = fftIfilt(64787:137078,:);

%randomly sample to lower corpus count, produces trials x samples
cfFiltPowerSub = datasample(cfNorm, 5000);
cpFiltPowerSub = datasample(cpNorm, 5000);
ifFiltPowerSub = datasample(ifNorm, 5000);
ipFiltPowerSub = datasample(ipNorm, 5000);

%concatenate to create correct and incorrect subsamples
cFiltPowerSub = vertcat(cfFiltPowerSub, cpFiltPowerSub);
iFiltPowerSub = vertcat(ifFiltPowerSub, ipFiltPowerSub);

% 
figure(1),clf
plot(fax,fftPow),hold on              % correct trials
plot(fax,fftPow_wrong)                      % incorrect trials
set(gca,'box','off','FontSize',24,'Xlim',[fax(1);fax(nyq)])

%change data to trials x samples and setup parameters, take baseline
%normalized raw data (baseline subtracted from delay)
cNorm = cNorm';
iNorm = iNorm';
s = 1000;
nframes = size(cNorm,2);
nyq = sr/2;
frStep = nyq/(nframes/2);
fax = 0:frStep:(nyq-frStep);
%compute power
fftd = fft(cNorm,nframes,2);                       % apply discrete fourier transform
%fftPow = mean(abs(fftd).^2);                  % convert to power estimates and avg across trials
fftPow = abs(fftd).^2;                  % convert to power estimates for all trials
fftPow = fftPow(:,1:size(fftPow,2)/2); %get rid of frequencies over nyquist
%extract frequencies 0-20Hz
freqIdx = find(fax>=0 & fax<20);
freqfft = fax(freqIdx(1):freqIdx(end));
fftCfilt = fftC(:,freqIdx(1):freqIdx(end));
fftIfilt = fftI(:,freqIdx(1):freqIdx(end

%extract phase
phaseI = angle(fftd);
phaseIfilt = phaseI(:,freqIdx(1):freqIdx(end));
phaseC = angle(fftd);
phaseCfilt = phaseC(:,freqIdx(1):freqIdx(end));
figure
%plot(freqfft,mean(phaseCfilt),freqfft,mean(phaseIfilt),':','LineWidth',2);
plot(freqfft,phaseCfilt(1:5,:),freqfft,phaseIfilt(1:5,:),':','LineWidth',2);
grid on
xlabel('Frequency(Hz)');
ylabel('Phase (rad)');
title('Butterworth 0-20Hz Filtered Delay-period Phase');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[freqfft(1);freqfft(end)])

%reduce dataset
%split into frontal and parietal
cfNorm = phaseCfilt(1:162922,:);
cpNorm = phaseCfilt(162923:352967,:);
ifNorm = phaseIfilt(1:64786,:);
ipNorm = phaseIfilt(64787:137078,:);

%randomly sample to lower corpus count, produces trials x samples
cfFiltPhaseSub = datasample(cfNorm, 5000);
cpFiltPhaseSub = datasample(cpNorm, 5000);
ifFiltPhaseSub = datasample(ifNorm, 5000);
ipFiltPhaseSub = datasample(ipNorm, 5000);

%concatenate to create correct and incorrect subsamples
cFiltPhaseSub = vertcat(cfFiltPhaseSub, cpFiltPhaseSub);
iFiltPhaseSub = vertcat(ifFiltPhaseSub, ipFiltPhaseSub);