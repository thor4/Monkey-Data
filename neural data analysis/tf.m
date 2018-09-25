%% Step 1: Get Signal

%initialize variables
srate = 1000; % 1,000Hz
% -500:0 baseline, 1:500 sample, 501:1310 delay (all in ms)
wavet = -.5:1/srate:.5; % in seconds 
min_freq =  3;
max_freq = 235;
num_frex = 120;
N     = length(wavet);
hz    = linspace(0,srate/2,floor(N/2)+1); % pos frequencies (not neg) up to Nyquist
hzf = linspace(0,srate,N); % total frequencies from fft
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(12),num_frex)./(2*pi*frex); %width of gaussian
fwhm = logspace(.8,.7,nfrex); %logspace or linspace for fwhm?

%% Make wavelets, ensure they taper to 0 at either end in time domain and
% ensure Gaussian in frequency domain

% Cohen's code from Google group (using n-cycles)
srate = 1000;
frex  = linspace(.5,200,50);
s2 = 2*(linspace(5,25,length(frex))./(2*pi*frex)).^2;

wavet = -.5:1/srate:.5;

for fi=1:length(frex)
    waves(fi,:) = exp(2*1i*pi*frex(fi)*wavet).*exp(-wavet.^2/s2(fi));
end

subplot(211), plot(wavet,real(waves))

% fft all the wavelets
for ffti=1:size(waves,1)
    waves_fft(ffti,:) = fft(waves(ffti,:));
end

subplot(212)
plot(linspace(0,srate,length(wavet)),abs(waves_fft).^2)
set(gca,'xlim',[30 62])
% ----------------

%% now using fwhm-specified in time domain
frex  = logspace(log10(3.5),log10(200),50);
fwhm = logspace(log10(.350),log10(.025),length(frex)); % in seconds
midp = dsearchn(wavet',0);
wavpts = length(wavet);
hz = linspace(0,srate/2,floor(wavpts/2)+1);
% outputs
empfwhmT = zeros(length(frex),1);
% loop over frequencies
for fi=1:length(frex)
    % create the Gaussian using the FWHM formula (equation 3)
    gwin = exp( (-4*log(2)*wavet.^2) ./ fwhm(fi)^2 );
    % measure the empirical fwhm (time domain)
    empfwhmT(fi) = wavet(midp-1+dsearchn(gwin(midp:end)',.5)) - ...
    wavet(dsearchn(gwin(1:midp)',.5));
end 
% empirical is pretty close to specified
% now make the wavelets and visualize them in time and frequency domains in
% addition to computing their fwhm 

empfwhmF = zeros(length(frex),1);

for fi=1:length(frex)
    wavelets(fi,:) = exp(2*1i*pi*frex(fi)*wavet).*exp( (-4*log(2)*wavet.^2) ./ fwhm(fi)^2 );
    wavelets_fft(fi,:) = fft(wavelets(fi,:)); % frequency domain
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
subplot(211), plot(wavet,real(wavelets))
% plot(time,signal,'k','linew',2)
xlabel('Time (s)'), ylabel('Amplitude (gain)')
title('Time domain')
% plot wavelets in frequency domain- ensure all are symmetric about their
% peak frequency and they taper to 0 on the ends
subplot(212)
plot(hz,abs(wavelets_fft(:,1:length(hz))).^2)
set(gca,'xlim',[0 250])
xlabel('Frequency (Hz)'), ylabel('Normalized Power')
title('Frequency domain')

% plot FWHM in time domain to show how much temporal smoothing occurs
figure(2), clf
subplot(211), plot(frex,empfwhmT*1000,'o:','markersize',8,'markerfacecolor','w','linew',2)
set(gca,'xlim',[0 max(frex)*1.05],'ylim',[0 max(empfwhmT)*1000*1.05])
xlabel('Wavelet frequency (Hz)'), ylabel('Empirical FWHM (ms)')
title('Time domain')
% plot FWHM in frequency domain to show how much spectral smoothing occurs
subplot(212)
plot(frex,empfwhmF,'s:','markersize',10,'markerfacecolor','w','linew',2)
set(gca,'xlim',[0 max(frex)*1.05],'ylim',[0 max(empfwhmF)*1.05])
xlabel('Wavelet frequency (Hz)'), ylabel('Empirical FWHM (Hz)')
title('Frequency domain')


% make wavelet
for fi=1:num_frex
    wavelet(fi,:) = exp(2*1i*pi*frex(fi).*wavet) .* exp(-wavet.^2./(2*(s(fi)^2)));
end

% now plot wavelets in groups of 10 at a time..
figure
for fi=1:5
    subplot(5,1,fi)
    plot(wavet,real(wavelet(fi,:)))
    title(sprintf('Wavelet at %fHz, s=%f number of cycles',frex(fi),s(fi)*(2*pi*frex(fi))));
end


%% check the frequency representation of the wavelets

figure(5), clf

% hz vector
hz = linspace(0,srate/2,floor(length(wavtime)/2)-1);

for fi=1:nFrex
    
    % create wavelet
    cmw = exp(  2*1i*pi*frex(fi)*wavtime - (wavtime.^2)/(2*s(fi)^2) );
    cmwX = fft(cmw,nConv);
    
    % plot wavelet and its power spectrum
    subplot(211), plot(real(cmw))
    subplot(212), plot(hz,2*abs(cmwX(1:length(hz))));
    
    title([ 'Frequency = ' num2str(frex(fi)) ])
    set(gca,'xlim',[0 80])
    pause % only the user can continue
end



%% Time-frequency decomposition via chronux
%load monkey data
load('mGoodStableRule1PingRej-split_by_Day_BehResp_and_Chan.mat') 
mwindow = 0.5; %in seconds
%allow for a step size with 75% overlap
movingwin=[mwindow mwindow*.05]; % set the moving window dimensions, 500ms window 50ms step
params.Fs=1000; % sampling frequency
params.fpass=[3 60]; % frequency of interest
params.tapers=[5 9]; % tapers
params.trialave=1; % average over trials
params.err=0; % no error computation

%initialize all variables
i=1; %day 1
j=1; %resp 1
k=1; %chan 1
l=1; %trial 1
resp = fieldnames(monkey(1).day(1));
monkeyTime = -500:1:1310; %total trial time (ms)
times2save = -500+(movingwin(1)*1000)/2-1:movingwin(2)*1000:1310-(movingwin(1)*1000)/2; 
%allows room for half of time_window on front & back

% define baseline period
baselinetime = [ -400 -100 ]; % in ms

% convert baseline window time to indices
[~,baselineidx(1)]=min(abs(times2save-baselinetime(1)));
[~,baselineidx(2)]=min(abs(times2save-baselinetime(2)));

chan = fieldnames(monkey(1).day(i).(resp{j}));
chan_combo = combnk(chan,2); %all chan combinations
size(chan_combo,1); %total number of combinations
%pull out channel combination
data1 = monkey(1).day(i).(resp{j}).(chan_combo{k,1});
data2 = monkey(1).day(i).(resp{j}).(chan_combo{k,2});

%compute coherence (S1 & S2 only hold the power, not analytic signal)
[C,phi,S12,S1,S2,t,f] = cohgramc(data1',data2',movingwin,params);

% dB-baseline corrected
baseline_power = mean(S1(baselineidx(1):baselineidx(2),:,:),1);
dbconverted = 10*log10( bsxfun(@rdivide,S1,baseline_power));

zmin = mean(dbconverted(:)) - 2*std(dbconverted(:));
zmax = mean(dbconverted(:)) + 2*std(dbconverted(:));
zrange = round((zmax - zmin)/2);

%now plot, need baseline corrected dB code
% x = 1 x samples
% y = 1 x frequencies
% z = frequencies x samples
% contourf(x,y,z,...)
W_ = params.tapers(1)/movingwin(1);
figure
contourf(t,f,dbconverted','linecolor','none')
set(gca,'ytick',round(logspace(log10(f(1)),log10(f(end)),10)*100)/100,'yscale','log','clim',[-zrange zrange])
% set(gca,'xlim',[monkeyTime(1) monkeyTime(end)])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title(sprintf('Power via multitaper from Monkey %d, Day %d, Resp %s, Channel %s\n W = %dHz, moving window = %dms, step = %dms',i,day(i),resp{j},chan_combo{k,1},W_,movingwin(1)*1000,movingwin(2)*1000));
cbar = colorbar; set(get(cbar,'label'),'string','dB change from baseline');   
%title(sprintf('%s\n', labels{:}))

[P,T,F]=mtspecgramc(Correct.area8B,movingwin,params);
subplot(2,6,2)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Frontal 8B'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

%other plot
plot_matrix(P,T,F); xlabel([]); % plot spectrogram
%caxis([8 28]); 
colorbar;
set(gca,'FontName','Times New Roman','Fontsize', 14);
title({['LFP 1,  W=' num2str(params.tapers(1)/movingwin(1)) 'Hz']; ['moving window = ' num2str(movingwin(1)) 's, step = ' num2str(movingwin(2)) 's']});
ylabel('frequency Hz');

%contour plot
contourf(times2save,frex,squeeze(mi(i,:,:))-repmat(mean(mi(i,:,baseidx(1):baseidx(2)),3)',1,length(times2save)),40,'linecolor','none')
set(gca,'clim',[-.075 .075],'yscale','log','ytick',round(logspace(log10(frex(1)),log10(frex(end)),6)))
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), title('MI Based on Power')
c = colorbar; set(get(c,'label'),'string','MI (baseline subtracted)');    

%% Time-frequency decomposition via wavelets
%load monkey data