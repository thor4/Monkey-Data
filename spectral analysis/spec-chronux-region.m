movingwin=[0.5 0.05]; % set the moving window dimensions, window=500ms step=50ms, w=10Hz
params.Fs=1000; % sampling frequency
params.fpass=[0 60]; % frequency of interest
params.tapers=[5 9]; % tapers
params.trialave=1; % average over trials
params.err=0; % no error computation
%
%min = -5;
%max = 26;
tic
figure 
% maxDb=25; %Limit power range to 15 Db
[P,T,F]=mtspecgramc(correctFrontal,movingwin,params);
subplot(2,1,1)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Correct Frontal'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

[P,T,F]=mtspecgramc(correctParietal,movingwin,params);
subplot(2,1,2)
imagesc(T,F,10*log10(P')) %Plot power in dB
axis xy; title('Correct Parietal'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

% [P,T,F]=mtspecgramc(IncorrectFrontal,movingwin,params);
% subplot(2,4,3)
% imagesc(T,F,10*log10(P')) %Plot power in dB
% axis xy; title('Incorrect Frontal'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;
% 
% [P,T,F]=mtspecgramc(IncorrectParietal,movingwin,params);
% subplot(2,4,4)
% imagesc(T,F,10*log10(P')) %Plot power in dB
% axis xy; title('Incorrect Parietal'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

toc