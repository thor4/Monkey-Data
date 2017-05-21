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
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Frontal 8B'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.area6DR,movingwin,params);
subplot(2,6,3)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Frontal 6DR'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.area8AD,movingwin,params);
subplot(2,6,4)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Frontal 8AD'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areavPFC,movingwin,params);
subplot(2,6,5)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Frontal vPFC'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areadPFC,movingwin,params);
subplot(2,6,6)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Frontal dPFC'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areaLIP,movingwin,params);
subplot(2,6,7)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Parietal LIP'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areaMIP,movingwin,params);
subplot(2,6,8)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Parietal MIP'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areaPE,movingwin,params);
subplot(2,6,9)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Parietal PE'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areaPG,movingwin,params);
subplot(2,6,10)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Parietal PG'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(Correct.areaPEC,movingwin,params);
subplot(2,6,11)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Parietal PEC'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

toc

plot_matrix(P,T,F); xlabel([]); % plot spectrogram
%caxis([8 28]); 
colorbar;
set(gca,'FontName','Times New Roman','Fontsize', 14);
title({['LFP 1,  W=' num2str(params.tapers(1)/movingwin(1)) 'Hz']; ['moving window = ' num2str(movingwin(1)) 's, step = ' num2str(movingwin(2)) 's']});
ylabel('frequency Hz');