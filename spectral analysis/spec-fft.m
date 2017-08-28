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

% 
figure(1),clf
plot(fax,fftPow),hold on              % correct trials
plot(fax,fftPow_wrong)                      % incorrect trials
set(gca,'box','off','FontSize',24,'Xlim',[fax(1);fax(nyq)])

%change data to trials x samples and setup parameters
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
subplot(2,1,1)
plot(freqfft,phaseCfilt,'LineWidth',1);
xlabel('Frequency(Hz)');
ylabel('Phase (rad)');
title('Butterworth 0-20Hz Filtered Delay-period Phase for Correct Trials');
legend('Correct Trials','Incorrect Trials');
set(gca,'box','off','Xlim',[freqfft(1);freqfft(end)])
subplot(2,1,2)
plot(freqfft,phaseIfilt,'LineWidth',1);
grid on
title('Butterworth 0-20Hz Filtered Delay-period Phase for Incorrect Trials');
xlabel('Frequency(Hz)');
ylabel('Phase (rad)');
set(gca,'box','off','Xlim',[freqfft(1);freqfft(end)])