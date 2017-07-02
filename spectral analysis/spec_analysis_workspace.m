fs = 1000;
figure
% pwelch(lfp_data(1,:),128,127,128,fs); xlabel('Frequency (Hz)')
[pxx,f] = pwelch(lfp_data(1,:),128,127,128,fs); xlabel('Frequency (Hz)')


% Chronux time spectrum
params.Fs = 1000;
params.fpass = [0 100];
params.pad = 2;
[S,f] = mtspectrumc(lfp_data(1,:),params);
figure
plot_vector(S,f)

% Chronux time frequency spectrum
movingwin=[0.5 0.05]; % set the moving window dimensions
params.Fs=1000; % sampling frequency
%params.fpass=[0 60]; % frequency of interest, Default all frequencies 
% between 0 and Fs/2
params.tapers=[5 9]; % tapers
params.trialave=1; % average over trials
params.err=0; % no error computation
[P,T,F]=mtspecgramc(data,movingwin,params);

% Matlab
[s,f,t] = spectrogram(lfp_data(1,:),128,127,128,1000,'yaxis');
ylim([0 300])
yticks([0 4 8 12 30 50 100 150 200 250 300])
title('Spectrogram')