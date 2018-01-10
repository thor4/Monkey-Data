movingwin=[0.5 0.05]; % set the moving window dimensions
params.Fs=1000; % sampling frequency
params.fpass=[0 60]; % frequency of interest
params.tapers=[5 9]; % tapers
params.trialave=1; % average over trials
params.err=0; % no error computation
%
tic
% maxDb=25; %Limit power range to 15 Db
[P,T,F]=mtspecgramc(Correct.area9L,movingwin,params);

figure %not in chapter
subplot(2,6,1)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Frontal 9L'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.area8B,movingwin,params);
subplot(2,6,2)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Frontal 8B'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.area6DR,movingwin,params);
subplot(2,6,3)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Frontal 6DR'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.area8AD,movingwin,params);
subplot(2,6,4)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Frontal 8AD'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areavPFC,movingwin,params);
subplot(2,6,5)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Frontal vPFC'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areadPFC,movingwin,params);
subplot(2,6,6)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Frontal dPFC'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areaLIP,movingwin,params);
subplot(2,6,7)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Parietal LIP'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areaMIP,movingwin,params);
subplot(2,6,8)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Parietal MIP'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areaPE,movingwin,params);
subplot(2,6,9)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Parietal PE'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areaPG,movingwin,params);
subplot(2,6,10)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Parietal PG'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areaPEC,movingwin,params);
subplot(2,6,11)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Parietal PEC'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

toc

plot_matrix(P,T,F); xlabel([]); % plot spectrogram
%caxis([8 28]); 
colorbar;
set(gca,'FontName','Times New Roman','Fontsize', 14);
title({['LFP 1,  W=' num2str(params.tapers(1)/movingwin(1)) 'Hz']; ['moving window = ' num2str(movingwin(1)) 's, step = ' num2str(movingwin(2)) 's']});
ylabel('frequency Hz');

%time frequency coherence (uses mtspecgramc)
%from Salazar supp: To evaluate the time course of the amplitude and synchronization between simultaneously
%recorded fronto-parietal LFPs, we calculated the time-frequency power and coherency spectra (200 msec. 
%sliding window stepped at 50 ms) using multi-taper spectral analysis (Mitra and Pesaran, 1999) implemented 
%in the Chronux toolbox (Bokil et al., 2010) (http://www.chronux.org/chronux/) with three tapers and a 
%time-bandwidth product of two (frequency resolution of 3.12 Hz and frequency smoothing of 12 Hz). 

window = .500; %in s
step = .050; %in s
movingwin=[window step]; % set the moving window dimensions, step size
params.Fs=1000; % sampling frequency
params.fpass=[2 60]; % frequency of interest
params.tapers=[3 5]; % tapers
params.trialave=1; % average over trials
params.err=0; % no error computation
%
tic
% maxDb=25; %Limit power range to 15 Db
[C,phi,S12,S1,S2,t,f]=cohgramc(lfp_data(1,:)',lfp_data(9,:)',movingwin,params);
toc

%spectrogram
figure
subplot(121)
contourf(t,f,S1',40,'linecolor','none')
colorbar
%set(gca,'clim',[-3 3],'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
%title('Logarithmic frequency scaling')

subplot(122)
contourf(t,f,S2',40,'linecolor','none')
colorbar
% subplot(122)
% contourf(EEG.times,frex,eegpower,40,'linecolor','none')
% set(gca,'clim',[-3 3],'xlim',[-200 1000])
% title('Linear frequency scaling')

%baseline normalization
 % Average power over trials (this code performs baseline transform,
 % which you will learn about in chapter 18)
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials
baseidx = dsearchn(EEG.times',[-500 -200]');
temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));


