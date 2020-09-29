%% Step 1: Define Parameters

%initialize variables
srate = 1000; % 1,000Hz
% -500:0 baseline, 1:500 sample, 501:1310 delay (all in ms)
wavet = -.5:1/srate:.5; % in seconds 
min_freq = 3; %in Hz
max_freq = 200;
num_frex = 50;
min_fwhm = .350; % in seconds (350ms)
max_fwhm = .050; % in seconds (50ms)
wavpts = length(wavet);
%there are N/2+1 frequencies between 0 and srate/2:
hz = linspace(0,srate/2,floor(wavpts/2)+1); % pos frequencies (not neg) up to Nyquist
frex = logspace(log10(min_freq),log10(max_freq),num_frex); %total num of freq's
% s    = logspace(log10(3),log10(12),num_frex)./(2*pi*frex); %width of gaussian
%fwhm of gaussian windows used to create wavelets logarithmically spaced:
fwhm = logspace(log10(min_fwhm),log10(max_fwhm),length(frex)); % in seconds


%% Step 2: Make wavelets, ensure they taper to 0 at either end in time domain and
% ensure Gaussian in frequency domain

% create wavelets using equation 3 from Cohen, 2018: "A better way to define
% and describe Morlet wavelets for time-frequency analysis"

% now using fwhm-specified in time domain
midp = dsearchn(wavet',0); %identify 0ms on wavelet time axis
% init outputs for empirical fwhm (time domain)
empfwhmT = zeros(length(frex),1);
for fi=1:length(frex) %loop over frequencies
    % create the Gaussian using the FWHM formula (equation 3):
    gwin = exp( (-4*log(2)*wavet.^2) ./ fwhm(fi)^2 ); 
    % measure the empirical fwhm (time domain):
    empfwhmT(fi) = wavet(midp-1+dsearchn(gwin(midp:end)',.5)) - ...
    wavet(dsearchn(gwin(1:midp)',.5));
end 
% empirical approximates specified (empfwhmT ~= fwhm)
% now make the wavelets and visualize them in time and frequency domains in
% addition to computing their fwhm 

empfwhmF = zeros(length(frex),1); %init outputs for emp fwhm (freq domain)

for fi=1:length(frex)
    %create wavelets using eq. 3: sin .* gaussian window
    wavelets(fi,:) = exp(2*1i*pi*frex(fi)*wavet).*exp( (-4*log(2)*wavet.^2) ./ fwhm(fi)^2 );
    wavelets_fft(fi,:) = fft(wavelets(fi,:)); % frequency domain (complex)
    wavelets_fft(fi,:) = wavelets_fft(fi,:)./max(wavelets_fft(fi,:)); % normalize 
    % measure the empirical fwhm (frequency domain)
    magnitude = abs(wavelets_fft(fi,1:length(hz)));
     % find left and right 1/2
    [~,peakx]  = max(wavelets_fft(fi,:)); 
    [~,left5]  = min(abs(wavelets_fft(fi,1:peakx)-.5));
    [~,right5] = min(abs(wavelets_fft(fi,peakx:end)-.5));
    right5 = right5+peakx-1;
    empfwhmF(fi) = hz(right5)-hz(left5);
end

% plot wavelets in time domain- ensure all taper to 0 (or very close)
figure(1), clf
subplot(211)
plot(wavet,real(wavelets))
axis off
% export_fig test2.png -transparent % no background
% plot(time,signal,'k','linew',2)
xlabel('Time (s)'), ylabel('Amplitude (gain)')
title('Time domain')
% plot wavelets in frequency domain- ensure all are symmetric about their
% peak frequency and they taper to 0 on the ends
subplot(212)
plot(hz,abs(wavelets_fft(:,1:length(hz))).^2)
set(gca,'xlim',[0 225])
xlabel('Frequency (Hz)'), ylabel('Normalized Power')
title('Frequency domain')

% plot FWHM in time domain to show how much temporal smoothing occurs
figure(2), clf
subplot(211)
plot(frex,empfwhmT*1000,'o:','markersize',8,'markerfacecolor','w','linew',2)
% semilogx(frex,empfwhmT*1000,'o:','markersize',8,'markerfacecolor','w','linew',2)
set(gca,'xlim',[0 max(frex)*1.05],'ylim',[0 max(empfwhmT)*1000*1.05])
xlabel('Wavelet frequency (Hz)'), ylabel('Empirical FWHM (ms)')
title('Time domain')
% plot FWHM in frequency domain to show how much spectral smoothing occurs
subplot(212)
plot(frex,empfwhmF,'s:','markersize',10,'markerfacecolor','w','linew',2)
% semilogx(frex,empfwhmF,'s:','markersize',10,'markerfacecolor','w','linew',2)
set(gca,'xlim',[0 max(frex)*1.05],'ylim',[0 max(empfwhmF)*1.05])
xlabel('Wavelet frequency (Hz)'), ylabel('Empirical FWHM (Hz)')
title('Frequency domain')

% PREP THESE FIGUERS FOR PAPER