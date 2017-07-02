fs = 1000;
figure
% pwelch(lfp_data(1,:),128,127,128,fs); xlabel('Frequency (Hz)')
[pxx,f] = pwelch(lfp_data(1,:),128,127,128,fs); xlabel('Frequency (Hz)')


% Chronux
params.Fs = 1000;
params.fpass = [0 100];
params.pad = 2;
[S,f] = mtspectrumc(lfp_data(1,:),params);
figure
plot_vector(S,f)

% Matlab
[s,f,t] = spectrogram(lfp_data(1,:),128,127,128,1000,'yaxis');
ylim([0 300])
yticks([0 4 8 12 30 50 100 150 200 250 300])
title('Spectrogram')