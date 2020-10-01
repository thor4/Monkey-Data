%% Step 1: Define Parameters

%initialize variables
srate = 1000; % 1,000Hz
% -500:0 baseline, 1:500 sample, 501:1310 delay (all in ms)
wavet = -1:1/srate:1; % in seconds (length doesn't matter since multiplying by zeroes won't affect convolution)
min_freq = 4; %in Hz (need several cycles in an epoch, these epochs are 500ms min so 4Hz = 2 cycles)
max_freq = 200; %may be no reason to go this high, try 100 after
num_frex = 50; %try 25-35 when only going to 100Hz, better for statistics mult comp corr, less smooth spectrogram
min_fwhm = .400; % in seconds (350ms)
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
plot(wavet,real(wavelets),'linew',2)
xlabel('Time (s)'), ylabel('Amplitude (gain)')
text(-0.72,1.2,'A','fontsize',35); box off
title('Time domain'); ax=gca; ax.FontSize = 25;
% plot wavelets in frequency domain- ensure all are symmetric about their
% peak frequency and they taper to 0 on the ends
subplot(212)
plot(hz,abs(wavelets_fft(:,1:length(hz))).^2,'linew',2)
set(gca,'xlim',[0 225]); text(-23,1.08,'B','fontsize',35);
xlabel('Frequency (Hz)'), ylabel('Normalized Power')
title('Frequency domain'); ax=gca; ax.FontSize = 25; box off
% export_fig('wavelets','-png','-transparent'); %save transparent pdf in pwd

% plot FWHM in time domain to show how much temporal smoothing occurs
figure(2), clf
subplot(211)
% plot(frex,empfwhmT*1000,'o:','markersize',8,'markerfacecolor','w','linew',2)
% semilogx(frex,empfwhmT*1000,'o:','markersize',8,'markerfacecolor','w','linew',2)
semilogx(frex,empfwhmT*1000,'.','markersize',25,'markerfacecolor','b')
% set(gca,'xlim',[0 max(frex)*1.05],'ylim',[0 max(empfwhmT)*1000*1.05])
set(gca,'ylim',[0 max(empfwhmT)*1000*1.05],'xminortick','off'); 
xticks(round(frex(1:3:end),1)); ylabel('FWHM (ms)')
text(2.2,max(empfwhmT)*1000*1.05+35,'A','fontsize',35);
title('Temporal Resolution log-lin'); ax=gca; ax.FontSize = 25;
ax.XTickLabel=[]; box off

% plot FWHM in frequency domain to show how much spectral smoothing occurs
subplot(212)
% plot(frex,empfwhmF,'s:','markersize',10,'markerfacecolor','w','linew',2)
% semilogx(frex,empfwhmF,'s:','markersize',8,'markerfacecolor','w','linew',2)
semilogx(frex,empfwhmF,'.','markersize',25,'markerfacecolor','b')
set(gca,'ylim',[0 max(empfwhmF)*1.05],'xminortick','off')
text(2.2,max(empfwhmF)*1.05+1,'B','fontsize',35);
xticks(round(frex(1:3:end),1)); box off
% set(gca,'xlim',[0 max(frex)*1.05],'ylim',[0 max(empfwhmF)*1.05],'xtick',frex)
xlabel('Wavelet frequency (Hz)'), ylabel('FWHM (Hz)')
title('Spectral Resolution log-lin'); ax=gca; ax.FontSize = 25;
% export_fig('fwhm smoothing','-png','-transparent'); %save transparent pdf in pwd

% Prepped figs for paper- check after the T/F analysis, may need
% to change parameters of wavelets to extract analytic signal in diff ways

%% Create analytic signal using the complex morlet wavelet family for both 
% monkeys, both conditions, all days, all channels and all trials

% load lfp's from:
% \OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data

load('time_domain-m1.mat') %dataM1goodCorR1 & dataM1goodIncR1
load('time_domain-m2.mat') %dataM2goodCorR1 & dataM2goodIncR1

%504 samples in baseline
%505 samples in cue
%811 samples in delay
%274 samples in match (check this)
%2094 total samples across all chans

% begin definining convolution parameters
n_wavelet = length(wavet);
half_of_wavelet_size = floor(n_wavelet/2)+1;
monkeyN = 2; % which monkey (1 or 2)
mAchans = {'8B', '9L', 'dPFC', 'vPFC', 'LIP', 'MIP', 'PEC', 'PG'};
mAareas = {'a8B', 'a9L', 'adPFC', 'avPFC', 'aLIP', 'aMIP', 'aPEC', 'aPG'};
mBchans = {'6DR', '8AD', '8B', 'dPFC', 'LIP', 'PE', 'PEC', 'PG'};
mBareas = {'a6DR', 'a8AD', 'a8B', 'adPFC', 'aLIP', 'aPE', 'aPEC', 'aPG'};

% define trial timeline
signalt = -.504:1/srate:1.589; %504 (nonzero sample) + 811 (delay) + 274 (match)=1589ms
signalt = -.5:1/srate:1.31; % in seconds
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1489; % in ms
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');
% define baseline time
baset = [-.4 -.1]; % in seconds
baseidx = dsearchn(signalt',baset');

% setup response structs
monkey(monkeyN).correct=[];
monkey(monkeyN).incorrect=[];


