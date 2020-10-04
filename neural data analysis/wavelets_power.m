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

%package into data variable 1: correct & 2: incorrect, then remove
data.mAgoodR1(1) = dataM1goodCorR1; data.mAgoodR1(2) = dataM1goodIncR1; 
data.mBgoodR1(1) = dataM2goodCorR1; data.mBgoodR1(2) = dataM2goodIncR1; 
clear dataM1goodCorR1 dataM2goodCorR1 dataM1goodIncR1 dataM2goodIncR1

%504 samples in baseline
%505 samples in cue
%811 samples in delay
%274 samples in match (check this)
%2094 total samples across all chans

% begin definining convolution parameters
n_wavelet = length(wavet);
half_of_wavelet_size = floor(n_wavelet/2)+1;
% define trial timeline
signalt = -.504:1/srate:1.589; %504 (nonzero sample) + 811 (delay) + 274 (match)=1589ms
% signalt = -.5:1/srate:1.31; % in seconds
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1466; % in ms, 1466 = 505 (sample) + 811 (delay) + 150 (match)
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');
% define baseline time
baset = [-.4 -.1]; % in seconds
baseidx = dsearchn(signalt',baset');

monkeys=fieldnames(data');
monkey=monkeys(end) %testing
dday=alldays(end) %testing
chan=allchans(end) %testing

%one loop for each monkey
monkey=monkeys(1) %which monkey: 1 = mA / Clark, 2 = mB / Betty

tic
%correct + incorrect
for i=1:2 %1 is correct, 2 is incorrect
    alldays = fieldnames(data.(monkey{:})(i))'; %extract all days
    for dday=alldays %loop thru days
        allchans = size(data.(monkey{:})(i).(dday{:}).lfp,1);
        %remove erp to free up space
%         data.(monkey{:})(i).(dday{:}) = rmfield(data.(monkey{:})(i).(dday{:}),'erp');
        for chan=1:allchans %loop thru chans
            %rem single chan dim, convert to µV (1V = 10^6µV =1,000,000µV)
            %and confirm time x trials leftover
            signal = squeeze(data.(monkey{:})(i).(dday{:}).lfp(chan,:,:)).* 1e6;
            reflectsig_all = zeros(size(signal,1)+2*n_wavelet,size(signal,2)); %initialize reflected signals mat
            % reflect all trials
            for signalN=1:size(signal,2) %loop through trials
                reflectsig = [ signal(n_wavelet:-1:1,signalN); signal(:,signalN); signal(end:-1:end-n_wavelet+1,signalN); ];        
                reflectsig_all(:,signalN) = reflectsig;
            end
            % concatenate into a super-trial
            reflectsig_supertri = reshape(reflectsig_all,1,[]); % reshape to 1D time-trials
    %         % plot original signal (example from mB chan 23 trial 848
    %         figure(1), clf
    %         subplot(211)
    %         plot(signal(:,848)','LineWidth',2,'color','b')
    %         set(gca,'xlim',[0 numel(reflectsig)]-n_wavelet)
    %         ylabel('Voltage (\muV)'); title('Original Signal')
    %         ax=gca; ax.FontSize = 25; x1=xticklabels; 
    %         set(gca, 'XTickLabel', []); box off
    %         % plot reflected signal
    %         subplot(212)
    %         p21=plot(n_wavelet+1:length(signal(:,848))+n_wavelet,signal(:,848),...
    %             'LineWidth',3,'color','b');
    %         hold on
    %         p22=plot(reflectsig,'-','LineWidth',2,'color','k');
    %         p22.Color(4) = 0.35; %change transparency
    %         set(gca,'xlim',[0 numel(reflectsig)])
    %         title('Reflected Signal'); xlabel('Time step (ms)')
    %         ylabel('Voltage (\muV)')
    %         ax=gca; ax.FontSize = 25; xticklabels(x1); box off
    %         legend({'original';'reflected'},'FontSize',25,'Location','best','box','off')
    %         export_fig('reflected signal','-png','-transparent'); %save transparent pdf in pwd
            % step 1: finish defining convolution parameters
            n_data = length(reflectsig_supertri); % time*trials
            n_convolution = n_wavelet+n_data-1;
            % step 2: take FFTs
            fft_data = fft(reflectsig_supertri,n_convolution); % all trials for chan
            % which area is this chan
            area = char(data.(monkey{:})(i).(dday{:}).areas(chan));
            %init power mat: freqidx x time x trials:
            pow = zeros(length(frex),length(times2save),size(signal,2)); 
            parfor fi=1:length(frex)
                % FFT of wavelet
                fft_wavelet = fft(wavelets(fi,:),n_convolution);
                % step 3: normalize kernel by scaling amplitudes to one in the 
                % frequency domain. prevents amplitude from decreasing with 
                % increasing frequency. diff from 1/f scaling
                fft_wavelet = fft_wavelet ./ max(fft_wavelet);
                % step 4: point-wise multiply and take iFFT
                as = ifft( fft_data.*fft_wavelet ); % analytic signal
                % step 5: trim wings
                as = as(half_of_wavelet_size:end-half_of_wavelet_size+1);
                % step 6: reshape back to reflected time-by-trials
                as = reshape(as,size(reflectsig_all,1),size(reflectsig_all,2));
                % step 7: chop off the reflections
                as = as(n_wavelet+1:end-n_wavelet,:);
                % as is now a time x trial complex matrix
    %             % compute baseline power averaged over trials then timepoints
    %             basePow = mean( mean( abs( as(baseidx(1):baseidx(2),:).^2 ),2 ),1 );
                % store dB-norm'd down-sampled power for each frequency in freq
                % x time x trials
    %             dbpow(fi,:,:) = 10*log10( abs( as(times2saveidx,:) ) .^2 ./basePow); 
                % save raw power for ea. freqidx x down-sampled time x
                % trial
                pow(fi,:,:) = abs( as(times2saveidx,:) ) .^2;
                % mean( abs( as_ ).^2, 2);
%                     clear fft_wavelet as % start anew with these var's ea. loop
            end
            %save downsampled power as chan x freqidx x time x trials
            data.(monkey{:})(i).(dday{:}).power(chan,:,:,:) = pow;                
            clear pow % start anew with this var ea. loop
        end
%             find(squeeze(data.(monkey{:})(i).(dday{:}).power(13,:,:,:))~=0); %testing
        %remove lfp to free up space in data var since we're done w/ it
%         data.(monkey{:})(i).(dday{:}) = rmfield(data.(monkey{:})(i).(dday{:}),'lfp');
    end
end
toc

%759.832210 seconds. (lab pc) saved average power
% save avg power in var ie: data.mAgoodR1(1).d060406.pow (chan x freqidx x time) in :
% OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data
% also has ie: data.mAgoodR1(1).d060406.areas (areas for each chan)

%only the last channel is showing up 
%find(squeeze(data.mBgoodR1(2).d090625.pow(10,:,:,:))~=0); %testing
% find(squeeze(pow(19,:,:))~=0); %testing

squeeze(data.(monkey{:})(i).(dday{:}).power(13,:,:,:)); %testing 
%freqidx x downsampled time points x trials (monkey B correct day d091001)

%% initialize variables

% define trial timeline
signalt = -.504:1/srate:1.589; %504 (nonzero sample) + 811 (delay) + 274 (match)=1589ms
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1466; % in ms, 1466 = 505 (sample) + 811 (delay) + 150 (match)
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');
% define baseline time
baset = [-.4 -.1]; % in seconds
baseidx = dsearchn(signalt',baset');

% which areas depends on monkeyN
allchans = size(data.(monkey{:})(i).(dday{:}).power,1); %total # of chans
areas = string(data.(monkey{:})(i).(dday{:}).areas); %all areas


%% plotting the raw difference data

% contourf plot template:
% x = 1 x samples
% y = 1 x frequencies
% z = frequencies x samples
% contourf(x,y,z,...)

triggers = [0 0.505 .505+.811]; %cueonset, cueoffset, matchonset

figure(1), clf %open fig and maximize to prepare for viewing

powmBvid = VideoWriter('powmBvid'); %open video file
powmBvid.FrameRate = 5;  %can adjust this, 5 - 10 seems to work
open(powmBvid)

for chanN=1:allchans %loop thru chans
    %pull out correct & incorrect power avg over trials (frex x time/samples)
    cor = mean( squeeze( data.(monkey{:})(1).(dday{:}).power(chanN,:,:,:) ),3 ); 
    % inc = mean( squeeze( data.(monkey{:})(2).(dday{:}).power(chan,:,:,:) ),3 ); 
    % compute the difference in power between the two conditions
    % diffmap = cor - inc;
    figure(1), clf
    contourf(signalt(times2saveidx),frex,cor,100,'linecolor','none')
    yL = get(gca,'YLim'); line([0 triggers(2) triggers(3);0 triggers(2) triggers(3)],...
        yL,'Color','k','LineWidth',4,'LineStyle',':'); 
    set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
    xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
    text(signalt(200),yL(2)*.85,'baseline','color','w','fontsize',20); %baseline
    text(triggers(2)*0.35,yL(2)*.85,'cue','color','w','fontsize',20); %cue
    text(triggers(3)*0.65,yL(2)*.85,'delay','color','w','fontsize',20); %delay
    text(triggers(3)*1.02,yL(2)*.85,'match','color','w','fontsize',20); %delay
    lim = get(cbar,'Limits'); cbar.Ticks=lim;
    cbar.Label.String = 'Raw Power (\muV^2)'; pos = cbar.Label.Position; 
    cbar.Label.Position=[pos(1)-2.5 pos(2)];
    title(sprintf('%s Day %d / %d Chan %d / %d Area %s Correct',monkey{:},...
        find(ismember(alldays,dday{:})),length(alldays),chanN,allchans,areas{chanN}));
    ax=gca; ax.FontSize = 25;
    pause(0.01) %pause to grab frame
    frame = getframe(gcf); %get frame
    writeVideo(powmBvid, frame); %add frame to vid
end

close(powmBvid) %finish with vid
%final video is saved here: 
% OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data

%see about averaging over channels to get a day-level power spectrogram and
%adding that frame to vid
%build out the day loop
