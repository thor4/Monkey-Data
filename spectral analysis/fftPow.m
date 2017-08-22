function [fftPow, fax] = fftPow(sr, data, nframes) %sr=sampling rate, data should be in #trials x nframes

%data = incorrectBase'; %update depending on what you want to take the fft of

sr = sr;                               % ephys sampling rate (Hz)
nframes = nframes;                 % number of samples in each epoch you want to fft (Hz)
nyq = sr/2;                             % nyquist (Hz)
frStep = nyq/(nframes/2);               % frequency resolution based on your nyquist and sampling window
fax = 0:frStep:(nyq-frStep);            % index of all frequencies up to the nyquist 

% sub out "data" below for your actual data. Should be in #trials x nframes
% format. 

% also, you'll want to implement notch filters around line noise freq and 
% harmonics (guessing 60, 120, 180 if this was done in N. America).

%data = rand(100,1000);                  % generate fake data
fftd = fft(data,nframes,2);                       % apply discrete fourier transform
fftPow = mean(abs(fftd).^2);                  % convert to power estimates and avg across trials
fftPow = fftPow(1:length(fftPow)/2);    % ditch all frequencies above the nyquist from redundant part of DFT
end