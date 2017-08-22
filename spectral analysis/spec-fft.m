s = 1000;
[fftpower_cAvgDelay, cfaxd] = fftPow(s, correct', 810);
fftpower_iAvgDelay = fftPow(s, incorrect');
[fftpower_cAvgBase, cfax] = fftPow(s, correctBase', 810);
fftpower_iAvgBase = fftPow(s, incorrectBase', 810);

%baseline normalization
fftpower_cAvgNorm = fftpower_cAvgDelay ./ fftpower_cAvgBase;
fftpower_iAvgNorm = fftpower_iAvgDelay ./ fftpower_iAvgBase;

%only interested in frequencies 1-250Hz
freqIdx = find(cfax>=0 & cfax<250);
f = cfax(freqIdx(1):freqIdx(end));
fftpower_cAvgNorm = fftpower_cAvgNorm(freqIdx(1):freqIdx(end));
fftpower_iAvgNorm = fftpower_iAvgNorm(freqIdx(1):freqIdx(end));

% 
figure(1),clf
plot(fax,fftPow),hold on              % correct trials
plot(fax,fftPow_wrong)                      % incorrect trials
set(gca,'box','off','FontSize',24,'Xlim',[fax(1);fax(nyq)])