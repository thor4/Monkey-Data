%% Step 1: Define Parameters

%initialize variables
srate = 1000; % 1,000Hz
% -500:0 baseline, 1:500 sample, 501:1310 delay (all in ms)
wavet = -1:1/srate:1; % in seconds (length doesn't matter since multiplying by zeroes won't affect convolution)
min_freq = 4; %in Hz (need several cycles in an epoch, these epochs are 500ms min so 4Hz = 2 cycles)
max_freq = 100; %nothing above 100
num_frex = 35; %50 for 200Hz, 35 for 100Hz, better for statistics mult comp corr, less smooth spectrogram
min_fwhm = .400; % in seconds (350ms)
max_fwhm = .100; % in seconds (75ms)
wavpts = length(wavet);
%there are N/2+1 frequencies between 0 and srate/2:
hz = linspace(0,srate/2,floor(wavpts/2)+1); % pos frequencies (not neg) up to Nyquist
frex = logspace(log10(min_freq),log10(max_freq),num_frex); %total num of freq's
% s    = logspace(log10(3),log10(12),num_frex)./(2*pi*frex); %width of gaussian
%fwhm of gaussian windows used to create wavelets logarithmically spaced:
fwhm = logspace(log10(min_fwhm),log10(max_fwhm),length(frex)); % in seconds
empcycles = 2.667 * (fwhm.*frex); %est # of cycles captured for ea frex (Cohen 2018)


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
plot(wavet,real(wavelets),'linew',1)
xlabel('Time (s)'), ylabel('Amplitude (gain)')
text(wavet(1)*1.19,1.2,'A','fontsize',35); box off
title('Time domain'); ax=gca; ax.FontSize = 25;
% plot wavelets in frequency domain- ensure all are symmetric about their
% peak frequency and they taper to 0 on the ends
subplot(212)
plot(hz,abs(wavelets_fft(:,1:length(hz))).^2,'linew',2)
set(gca,'xlim',[0 120]); text(-hz(2)*24,1.08,'B','fontsize',35);
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
text(3,max(empfwhmT)*1000*1.05+35,'A','fontsize',35);
title('Temporal Resolution log-lin'); ax=gca; ax.FontSize = 25;
ax.XTickLabel=[]; box off

% plot FWHM in frequency domain to show how much spectral smoothing occurs
subplot(212)
% plot(frex,empfwhmF,'s:','markersize',10,'markerfacecolor','w','linew',2)
% semilogx(frex,empfwhmF,'s:','markersize',8,'markerfacecolor','w','linew',2)
semilogx(frex,empfwhmF,'.','markersize',25,'markerfacecolor','b')
set(gca,'ylim',[0 max(empfwhmF)*1.05],'xminortick','off')
text(3,max(empfwhmF)*1.05+1,'B','fontsize',35);
xticks(round(frex(1:2:end),1)); box off
% set(gca,'xlim',[0 max(frex)*1.05],'ylim',[0 max(empfwhmF)*1.05],'xtick',frex)
xlabel('Wavelet frequency (Hz)'), ylabel('FWHM (Hz)')
title('Spectral Resolution log-lin'); ax=gca; ax.FontSize = 25;
% export_fig('fwhm smoothing','-png','-transparent'); %save transparent pdf in pwd

% Prepped figs for paper- check after the T/F analysis, may need
% to change parameters of wavelets to extract analytic signal in diff ways

%% Step 2: Create analytic signal using the complex morlet wavelet family for both (step 1 was ERP)
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

monkeys=fieldnames(data)';

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
%             pow = zeros(length(frex),length(times2save),size(signal,2)); 
            basePow = zeros(length(frex),1)'; 
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
                % compute baseline power averaged over trials then timepoints
                % & save avg baseline power per freq component
                basePow(fi) = mean( mean( abs( as(baseidx(1):baseidx(2),:).^2 ),2 ),1 );
                % store dB-norm'd down-sampled power for each frequency in freq
                % x time x trials
%                 dbpow(fi,:,:) = 10*log10( abs( as(times2saveidx,:) ) .^2 ./basePow); 
                % save raw power for ea. freqidx x down-sampled time x
                % trial
%                 pow(fi,:,:) = abs( as(times2saveidx,:) ) .^2;
                % mean( abs( as_ ).^2, 2);
%                     clear fft_wavelet as % start anew with these var's ea. loop
            end
            %save downsampled power as chan x freqidx x time x trials
%             data.(monkey{:})(i).(dday{:}).power(chan,:,:,:) = pow; 
            %save avg baseline power for all frex as chan x freqidx
            data.(monkey{:})(i).(dday{:}).basepow(chan,:) = basePow; 
            clear basePow % start anew with this var ea. loop
        end
%             find(squeeze(data.(monkey{:})(i).(dday{:}).power(13,:,:,:))~=0); %testing
        %remove lfp to free up space in data var since we're done w/ it
        data.(monkey{:})(i).(dday{:}) = rmfield(data.(monkey{:})(i).(dday{:}),'lfp');
    end
end
toc

%afterwards, need to average the basePow for correct with incorrect to get
%a single number across conditions if/when doing the conditon comparison

mAgoodR1 = data.mAgoodR1; %pull out var to save it as file

%759.832210 seconds. (lab pc) saved average power
%~1.25 hours for correct+incorrect for monkey 2 on lab pc
% save all power trials in var ie: mAgoodR1.d060406.power (chan x freqidx x time x trials) in :
% OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data
% also has ie: mAgoodR1.d060406.areas (areas for each chan)

%find(squeeze(data.mBgoodR1(2).d090625.pow(10,:,:,:))~=0); %testing
% find(squeeze(pow(19,:,:))~=0); %testing

% squeeze(data.(monkey{:})(i).(dday{:}).power(13,:,:,:)); %testing 
%freqidx x downsampled time points x trials (monkey B correct day d091001)

%% Baseline normalize power

%init vars
min_freq = 4; %in Hz (need several cycles in an epoch, these epochs are 500ms min so 4Hz = 2 cycles)
max_freq = 100; %nothing above 100
num_frex = 35; %50 for 200Hz, 35 for 100Hz, better for statistics mult comp corr, less smooth spectrogram
frex = logspace(log10(min_freq),log10(max_freq),num_frex); %total num of freq's
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1466; % in ms, 1466 = 505 (sample) + 811 (delay) + 150 (match)

%load non-normalized power files
monkeys = {'mA','mB'}; %setup monkeys array
mi = 2; %choose which monkey: mA=1, mB=2

monkey=monkeys{mi}; %load data for chosen monkey
if monkey=='mA'
    load('mAgoodR1_pow_erp_lfp.mat'); %monkey A data
    mData = mAgoodR1; clear mAgoodR1
    load('mAgoodR1_basepow_erp.mat'); 
    data=mAgoodR1; clear mAgoodR1
else
    load('mBgoodR1_pow_erp_lfp.mat'); %monkey B data
    mData = mBgoodR1; clear mBgoodR1
    load('mBgoodR1_basepow_erp.mat'); 
    data=mBgoodR1; clear mBgoodR1
end

alldays = fieldnames( mData )'; %extract all days
for dayN=alldays
    allchans = size(mData(1).(dayN{:}).power,1); %total # of chans
    areas = string(mData(1).(dayN{:}).areas); %all areas
    %init dB baseline normalized power array [resp x chan x freqidx x
    %times2saveidx]
    dbnPow = zeros(2,length(allchans),length(frex),length(times2save));
    for chanN = 1:allchans
        for fi = 1:length(frex)
            %extract baseline-normalized power for chanN x freq component
            basePow = (data(1).(dayN{:}).basepow(chanN,fi) + ...
                data(2).(dayN{:}).basepow(chanN,fi))/2; %avg of cor & inc
            for resp=1:2 %correct then incorrect
                %average raw power over trials for chanN x freq combo = 187
                %times2saveidx
                avgPow = squeeze( mean( mData(resp).(dayN{:}).power(chanN,fi,:,:),4 ) );
                %dB baseline normalize the power for cor then inc
                dbnPow(resp,chanN,fi,:) = 10*log10( avgPow ./ basePow );
            end
        end
    end
%             find(dbnPow(2,chanN-8,:,:)~=0); %testing
            size(dbnPower) %testing
%             size(dbnPow) %testing
            find(dbnPower(2,chanN-3,:,:)~=0); %testing
    %now package it all up
    dayy=find(ismember(alldays,dayN{:}));
    if dayy==1 
        dbnPower = dbnPow;
        dbn_pow.chan = [1:allchans]; %1 x allchans vector
        dbn_pow.chant = repmat(allchans,allchans,1)';
        dbn_pow.area = areas;
        dbn_pow.day = repmat(dayy,allchans,1)';
    else
        dbnPower = cat(2,dbnPower,dbnPow);
        dbn_pow.chan = cat(2,dbn_pow.chan,[1:allchans]);
        dbn_pow.chant = cat(2,dbn_pow.chant,repmat(allchans,allchans,1)');
        dbn_pow.area = cat(2,dbn_pow.area,areas);
        dbn_pow.day = cat(2,dbn_pow.day,repmat(dayy,allchans,1)');
    end
end
dbn_pow.dbnPower = dbnPower;

%save dbn_pow as m%goodR1_dBbasenormpow.mat in 
%\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\step 2 TF - power

%143 chans in mA
%318 chans in mB


%% initialize variables

srate=1000;
min_freq = 4; %in Hz (need several cycles in an epoch, these epochs are 500ms min so 4Hz = 2 cycles)
max_freq = 100; %nothing above 100
num_frex = 35; %50 for 200Hz, 35 for 100Hz, better for statistics mult comp corr, less smooth spectrogram
frex = logspace(log10(min_freq),log10(max_freq),num_frex); %total num of freq's
% define trial timeline
signalt = -.504:1/srate:1.589; %504 (nonzero sample) + 811 (delay) + 274 (match)=1589ms
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1466; % in ms, 1466 = 505 (sample) + 811 (delay) + 150 (match)
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');
% define baseline time
baset = [-.4 -.1]; % in seconds
baseidx = dsearchn(signalt',baset');


%% plotting the raw data

% contourf plot template:
% x = 1 x samples
% y = 1 x frequencies
% z = frequencies x samples
% contourf(x,y,z,...)

load('mAgoodR1_pow_erp_lfp.mat'); load('mAgoodR1_dBbasenormpow.mat'); %monkey A data
load('mBgoodR1_pow_erp_lfp.mat'); load('mBgoodR1_dBbasenormpow.mat'); %monkey B data

monkeys={'mAgoodR1(1)','mBgoodR1(1)','mAgoodR1(2)','mBgoodR1(2)'};

i=1; vidname=sprintf('powm%svid','A'); %mA correct trials only
i=2; vidname=sprintf('powm%svid','B'); %mB correct trials only
i=3; vidname=sprintf('powm%svid','A'); %mA incorrect trials only
i=4; vidname=sprintf('powm%svid','B'); %mB incorrect trials only

alldays = fieldnames(eval(monkeys{i}))'; %extract all days
triggers = [0 0.505 .505+.811]; %cueonset, cueoffset, matchonset

figure(1), clf %open fig and maximize to prepare for viewing

powvid = VideoWriter(vidname); %open video file
powvid.FrameRate = 5;  %can adjust this, 5 - 10 seems to work
open(powvid)

for dayN=alldays
    % which areas depends on monkeyN
    allchans = size(eval(monkeys{i}).(dayN{:}).power,1); %total # of chans
    areas = string(eval(monkeys{i}).(dayN{:}).areas); %all areas
    for chanN=1:allchans %loop thru chans
        %pull out correct & incorrect power avg over trials (frex x time/samples)
        cor = mean( squeeze( eval(monkeys{i}).(dayN{:}).power(chanN,:,:,:) ),3 ); 
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
        title(sprintf('%8s Day %d / %d Chan %d / %d Area %s Correct',monkeys{i}(1:8),...
            find(ismember(alldays,dayN{:})),length(alldays),chanN,allchans,areas{chanN}));
        ax=gca; ax.FontSize = 25;
        pause(0.01) %pause to grab frame
        frame = getframe(gcf); %get frame
        writeVideo(powvid, frame); %add frame to vid
    end
    %now calc day power spectrogram avg over all chans
    %pull out correct power avg over trials then chans (frex x time/samples)
    cor = squeeze( mean( mean( eval(monkeys{i}).(dayN{:}).power,4 ),1 ) ); 
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
    title(sprintf('%s Day %d / %d All %d Chans Across %d Areas Correct',monkeys{i}(1:8),...
        find(ismember(alldays,dayN{:})),length(alldays),allchans,length(unique(areas))));
    ax=gca; ax.FontSize = 25;
    pause(0.01) %pause to grab frame
    frame = getframe(gcf); %get frame
    writeVideo(powvid, frame); %add frame to vid
end

close(powvid) %finish with vid
%final video is saved here: 
% OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data

%now get full monkey correct ERP (change i val's and plot title + cbar
%limits for diff monkeys and for diffmap)
for dayN=alldays
    dayy = find(ismember(alldays,dayN{:})); %day number
    % which areas depends on monkeyN
    areas = string(eval(monkeys{i}).(dayN{:}).areas); %all areas
    %now calc day power spectrogram avg over all chans
    %pull out correct power avg over trials then chans (frex x time/samples)
    %build out day avg-over-chans matrix [day freqidx times2saveidx]
    cor(dayy,:,:) = squeeze( mean( mean( eval(monkeys{i}).(dayN{:}).power,4 ),1 ) ); 
    inc(dayy,:,:) = squeeze( mean( mean( eval(monkeys{i+2}).(dayN{:}).power,4 ),1 ) );
end
% find(cor(15,:,:)~=0); %testing
mAvgc = squeeze(mean(cor,1)); mAvgi=squeeze(mean(inc,1));
%baseline normalization: avg over chans= [freqidx x times2saveidx]
mAvgdbnc = squeeze( mean( dbn_pow.dbnPower(1,:,:,:),2 ) );
mAvgdbni = squeeze( mean( dbn_pow.dbnPower(2,:,:,:),2 ) );

%pull out frontal areas only
frontal = unique(dbn_pow.area); 
% frontal = [frontal(1:2), frontal(7:end)]; %mA
frontal = [frontal(1:3), frontal(end)]; %mB
fchanBool = ismember(dbn_pow.area,frontal);
sum(ismember(dbn_pow.area,frontal)) %146 total frontal chans for mB & 69 for mA
%[chans x freqidx x times2saveidx]
fchans = squeeze( dbn_pow.dbnPower(:,fchanBool,:,:) ); %frontal chans only
mAvgdbnf = squeeze( mean( fchans,2 ) ); %avg over all frontal chans
%leaves [resp(2) x freqidx(35) x times2saveidx(187)]



%set subtightplot function parameters up:
make_it_tight = true;
%set ([vert horiz](axes gap),[lower uppper](margins),[left right](margins))
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.05], [0.08 0.04], [0.06 0.035]);
if ~make_it_tight,  clear subplot;  end

figure(1), clf
% subplot(2,1,1)
% contourf(signalt(times2saveidx),frex,mAvgc-mAvgi,100,'linecolor','none')
contourf(signalt(times2saveidx),frex,mAvgc,100,'linecolor','none')
yL = get(gca,'YLim'); line([0 triggers(2) triggers(3);0 triggers(2) triggers(3)],...
    yL,'Color','k','LineWidth',4,'LineStyle',':'); 
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
set(gca,'xticklabel',[])
% xlabel('Time (s)')
ylabel('Frequency (Hz)'), cbar = colorbar; 
text(signalt(200),yL(2)*.85,'baseline','color','w','fontsize',20); %baseline
text(triggers(2)*0.35,yL(2)*.85,'cue','color','w','fontsize',20); %cue
text(triggers(3)*0.65,yL(2)*.85,'delay','color','w','fontsize',20); %delay
text(triggers(3)*1.02,yL(2)*.85,'match','color','w','fontsize',20); %match
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Raw Power (\muV^2)'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-2.5 pos(2)];
% cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('%s %d Days %d Frontal & Parietal Channels Correct',monkeys{i}(1:8),...
    length(alldays),size(dbn_pow.dbnPower,2)));
ax=gca; ax.FontSize = 20;
% export_fig('mB_raw_power_cor','-png','-transparent'); %save transparent pdf in pwd

%work on pulling out baseline normalized power from analytic signal
% subplot(2,1,2)
% contourf(signalt(times2saveidx),frex,squeeze(mAvgdbnf(1,:,:)-mAvgdbnf(2,:,:)),100,'linecolor','none') %diffmap frontal chans
% contourf(signalt(times2saveidx),frex,mAvgdbnc-mAvgdbni,100,'linecolor','none') %diffmap all chans
contourf(signalt(times2saveidx),frex,mAvgdbnc,100,'linecolor','none')
% contourf(signalt(times2saveidx),frex,squeeze(mAvgdbnf(1,:,:)),100,'linecolor','none') %cor frontal chans
yL = get(gca,'YLim'); line([0 triggers(2) triggers(3);0 triggers(2) triggers(3)],...
    yL,'Color','k','LineWidth',4,'LineStyle',':'); 
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
text(signalt(200),yL(2)*.85,'baseline','color','w','fontsize',20); %baseline
text(triggers(2)*0.35,yL(2)*.85,'cue','color','w','fontsize',20); %cue
text(triggers(3)*0.65,yL(2)*.85,'delay','color','w','fontsize',20); %delay
text(triggers(3)*1.02,yL(2)*.85,'match','color','w','fontsize',20); %delay
caxis([-1.5 1.5]) %set colorbar limits
% caxis([-0.25 0.25]) %set colorbar limits
% caxis([-0.4 0.4]) %set colorbar limits
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'dB change from baseline'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.65 pos(2)];
% cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('%s %d Days %d Frontal & Parietal Channels Correct',monkeys{i}(1:8),...
    length(alldays),size(dbn_pow.dbnPower,2)));
% title(sprintf('%s %d Days %d Frontal Channels Correct',monkeys{i}(1:8),...
%     length(alldays),size(fchans,2)));
ax=gca; ax.FontSize = 24;

export_fig('mB_dbn_power_cor-marks','-png','-transparent'); %save transparent png in pwd

% dim = ds2nfu([xo yo wdth hght]); %translate to normalized fig units
% annotation('rectangle',dim,'Color','black')

%% pull out and compute the regional day averages
% doing this to address the issue of channel correlation since that will
% effect generalizability according to Mike Cohen:
% "Your assumption here is that you are randomly sampling from circuits in 
% one animal, so the appropriate generalization is to neural populations in
% that monkey (and then for each monkey)."
% https://discuss.sincxpress.com/t/group-level-strategy-2a-questions/469/2

%init vars
srate = 1000; %1000Hz sampling rate
min_freq = 4; %in Hz (need several cycles in an epoch, these epochs are 500ms min so 4Hz = 2 cycles)
max_freq = 100; %nothing above 100
num_frex = 35; %50 for 200Hz, 35 for 100Hz, better for statistics mult comp corr, less smooth spectrogram
frex = logspace(log10(min_freq),log10(max_freq),num_frex); %total num of freq's
num_frex = 35; %50 for 200Hz, 35 for 100Hz, better for statistics mult comp corr, less smooth spectrogram
% define trial timeline
signalt = -.504:1/srate:1.589; %504 (nonzero sample) + 811 (delay) + 274 (match)=1589ms
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1466; % in ms, 1466 = 505 (sample) + 811 (delay) + 150 (match)
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');
triggers = [0 0.505 .505+.811]; %cueonset, cueoffset, matchonset


monkeys = {'mA','mB'}; %setup monkeys array
mi = 1; %choose which monkey: mA=1, mB=2

monkey=monkeys{mi}; %load data for chosen monkey
if monkey=='mA'
    load('mAgoodR1_dBbasenormpow.mat'); %monkey A data
else %mB
    load('mBgoodR1_dBbasenormpow.mat'); %monkey B data
end

fAreas = {'9L','8B','dPFC','vPFC','6DR','8AD'}; %frontal areas from both monkeys
pAreas = {'PEC','MIP','LIP','PG','PE'}; %parietal areas from both monkeys
mAreas = unique(dbn_pow.area); %all areas in current monkey
frontal = fAreas(ismember(fAreas,mAreas)); %frontal areas for current monkey
parietal = pAreas(ismember(pAreas,mAreas)); %parietal areas for current monkey

for day=1:max(dbn_pow.day) %pull out & avg over regions
    dayReg = zeros(2,2,length(frex),length(times2save)); %init day regional avg mat
    fchanBool = ismember(dbn_pow.area,frontal) & dbn_pow.day==day;
    fchans = dbn_pow.dbnPower(:,fchanBool,:,:); %frontal chans only
    mAvgdbnf = squeeze( mean( fchans,2 ) ); %avg over all frontal chans
    pchanBool = ismember(dbn_pow.area,parietal) & dbn_pow.day==day;
    pchans = dbn_pow.dbnPower(:,pchanBool,:,:); %parietal chans only
    mAvgdbnp = squeeze( mean( pchans,2 ) ); %avg over all parietal chans
    %leaves [resp(2) x freqidx(35) x times2saveidx(187)]
    dayReg(1,:,:,:)=mAvgdbnf; dayReg(2,:,:,:)=mAvgdbnp;
%     find(dayReg(2,:,:,:)~=0); %testing
    if day==1 %save in dbn_pow struct
        dbn_pow.dayRegAvg = dayReg;
        dbn_pow.dayReg = repmat(day,2,1)';
    else
        dbn_pow.dayRegAvg = cat(1,dbn_pow.dayRegAvg,dayReg);
        dbn_pow.dayReg = cat(2,dbn_pow.dayReg,repmat(day,2,1)');
    end
end

size(dbn_pow.dayRegAvg) %testing mA size 46 (23 days * 2 regional avg's)
% find(dbn_pow.dayRegAvg(end,:,:,:)~=0); %testing

%save as new files m%goodR1_dBbasenormpow_dayReg to indicate this struct
%includes the regional averages by day in: 
% OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\step 2 TF - power

%visualize some of the regional day averages
if monkey=='mA'
    load('mAgoodR1_dBbasenormpow_dayReg.mat'); %monkey A dayReg data
else %mB
    load('mBgoodR1_dBbasenormpow_dayReg.mat'); %monkey B dayReg data
end

day = randperm(max(dbn_pow.dayReg),1); %choose a random day
day = ismember(dbn_pow.dayReg,day); %setup bool of where day is
dayRegAvg = dbn_pow.dayRegAvg(day,:,:,:); %extract region-level power for day
dayRegAvgcf = squeeze(dayRegAvg(1,1,:,:)); %extract avg correct frontal signal
dayRegAvgif = squeeze(dayRegAvg(1,2,:,:)); %extract avg incorrect frontal signal
dayRegAvgcp = squeeze(dayRegAvg(2,1,:,:)); %extract avg correct parietal signal
dayRegAvgip = squeeze(dayRegAvg(2,2,:,:)); %extract avg incorrect parietal signal
size(dayRegAvg) %[2(regions) 2(resp) 35(frex) 187(time)]


%set subtightplot function parameters up:
make_it_tight = true;
%set ([vert horiz](axes gap),[lower uppper](margins),[left right](margins))
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.08 0.04], [0.055 0.03]);
if ~make_it_tight,  clear subplot;  end

figure(1), clf
subplot(2,2,1)
contourf(signalt(times2saveidx),frex,dayRegAvgcf,100,'linecolor','none')
yL = get(gca,'YLim'); line([0 triggers(2) triggers(3);0 triggers(2) triggers(3)],...
    yL,'Color','k','LineWidth',4,'LineStyle',':'); 
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
set(gca,'xticklabel',[])
% xlabel('Time (s)')
ylabel('Frequency (Hz)'), cbar = colorbar; 
text(signalt(200),yL(2)*.85,'baseline','color','w','fontsize',14); %baseline
text(triggers(2)*0.35,yL(2)*.85,'cue','color','w','fontsize',14); %cue
text(triggers(3)*0.65,yL(2)*.85,'delay','color','w','fontsize',14); %delay
text(triggers(3)*1.01,yL(2)*.85,'match','color','w','fontsize',14); %match
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'dB change from baseline'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.65 pos(2)];
% cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('%s Day %d Correct Frontal Avg of %d Channels',monkey,...
    unique(dbn_pow.dayReg(day)),...
    sum( dbn_pow.day==unique(dbn_pow.dayReg(day)) & ismember(dbn_pow.area,frontal) )));
ax=gca; ax.FontSize = 18;

subplot(2,2,2)
contourf(signalt(times2saveidx),frex,dayRegAvgif,100,'linecolor','none')
yL = get(gca,'YLim'); line([0 triggers(2) triggers(3);0 triggers(2) triggers(3)],...
    yL,'Color','k','LineWidth',4,'LineStyle',':'); 
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
set(gca,'yticklabel',[]), set(gca,'xticklabel',[])
% xlabel('Time (s)') ylabel('Frequency (Hz)'), 
cbar = colorbar; 
text(signalt(200),yL(2)*.85,'baseline','color','w','fontsize',14); %baseline
text(triggers(2)*0.35,yL(2)*.85,'cue','color','w','fontsize',14); %cue
text(triggers(3)*0.65,yL(2)*.85,'delay','color','w','fontsize',14); %delay
text(triggers(3)*1.01,yL(2)*.85,'match','color','w','fontsize',14); %match
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'dB change from baseline'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.65 pos(2)];
% cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('%s Day %d Incorrect Frontal Avg of %d Channels',monkey,...
    unique(dbn_pow.dayReg(day)),...
    sum( dbn_pow.day==unique(dbn_pow.dayReg(day)) & ismember(dbn_pow.area,frontal) )));
ax=gca; ax.FontSize = 18;

subplot(2,2,3)
contourf(signalt(times2saveidx),frex,dayRegAvgcp,100,'linecolor','none')
yL = get(gca,'YLim'); line([0 triggers(2) triggers(3);0 triggers(2) triggers(3)],...
    yL,'Color','k','LineWidth',4,'LineStyle',':'); 
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
% set(gca,'xticklabel',[])
xlabel('Time (s)') 
ylabel('Frequency (Hz)'), cbar = colorbar; 
text(signalt(200),yL(2)*.85,'baseline','color','w','fontsize',14); %baseline
text(triggers(2)*0.35,yL(2)*.85,'cue','color','w','fontsize',14); %cue
text(triggers(3)*0.65,yL(2)*.85,'delay','color','w','fontsize',14); %delay
text(triggers(3)*1.01,yL(2)*.85,'match','color','w','fontsize',14); %match
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'dB change from baseline'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.65 pos(2)];
% cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('%s Day %d Correct Parietal Avg of %d Channels',monkey,...
    unique(dbn_pow.dayReg(day)),...
    sum( dbn_pow.day==unique(dbn_pow.dayReg(day)) & ismember(dbn_pow.area,parietal) )));
ax=gca; ax.FontSize = 18;

subplot(2,2,4)
contourf(signalt(times2saveidx),frex,dayRegAvgip,100,'linecolor','none')
yL = get(gca,'YLim'); line([0 triggers(2) triggers(3);0 triggers(2) triggers(3)],...
    yL,'Color','k','LineWidth',4,'LineStyle',':'); 
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
set(gca,'yticklabel',[])
xlabel('Time (s)') 
%ylabel('Frequency (Hz)'), 
cbar = colorbar; 
text(signalt(200),yL(2)*.85,'baseline','color','w','fontsize',14); %baseline
text(triggers(2)*0.35,yL(2)*.85,'cue','color','w','fontsize',14); %cue
text(triggers(3)*0.65,yL(2)*.85,'delay','color','w','fontsize',14); %delay
text(triggers(3)*1.01,yL(2)*.85,'match','color','w','fontsize',14); %match
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'dB change from baseline'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.65 pos(2)];
% cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('%s Day %d Incorrect Parietal Avg of %d Channels',monkey,...
    unique(dbn_pow.dayReg(day)),...
    sum( dbn_pow.day==unique(dbn_pow.dayReg(day)) & ismember(dbn_pow.area,frontal) )));
ax=gca; ax.FontSize = 18;

export_fig('mB_dbn_power_RegAvg_day4','-png','-transparent'); %save transparent png in pwd

%first, verify this calculation is correct. next, make a movie of all days 
%for each monkey doing this. 

%% test hypotheses using one-sample directional t-tests

%init vars
srate = 1000; %1000Hz sampling rate
min_freq = 4; %in Hz (need several cycles in an epoch, these epochs are 500ms min so 4Hz = 2 cycles)
max_freq = 100; %nothing above 100
num_frex = 35; %50 for 200Hz, 35 for 100Hz, better for statistics mult comp corr, less smooth spectrogram
frex = logspace(log10(min_freq),log10(max_freq),num_frex); %total num of freq's
num_frex = 35; %50 for 200Hz, 35 for 100Hz, better for statistics mult comp corr, less smooth spectrogram
% define trial timeline
signalt = -.504:1/srate:1.589; %504 (nonzero sample) + 811 (delay) + 274 (match)=1589ms
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1466; % in ms, 1466 = 505 (sample) + 811 (delay) + 150 (match)
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');
% define baseline time (for db-normalized to baseline to make it normally
% distributed so it's appropriate for parametric statistics)
baset = [-.4 -.1]; % in seconds
baseidx = dsearchn(signalt',baset');

%set the alpha level
ntests = 3; %how many t-tests?
pval = 0.05/ntests; %.05 bonferonni-corrected for multiple comparisons

monkeys = {'mA','mB'}; %setup monkeys array
mi = 1; %choose which monkey: mA=1, mB=2

monkey=monkeys{mi}; %load data for chosen monkey
if monkey=='mA'
    load('mAgoodR1_dBbasenormpow.mat'); %monkey A data
    %hyp1: f & p regions show increased low-band LFP power delay
    x1=580; x2=1230; y1=4.3972; y2=18.1934; %coords of hyp1 rect ROI
    hyp1t = [find(times2save==x1):find(times2save==x2)]; %hyp1 timepoints
    hyp1f = [dsearchn(frex',y1):dsearchn(frex',y2)]; %hyp1 freq components
    x1=580; x2=1230; y1=21.986; y2=38.8008; %coords of hyp2 rect ROI
    hyp2t = [find(times2save==x1):find(times2save==x2)]; %hyp2 timepoints
    hyp2f = [dsearchn(frex',y1):dsearchn(frex',y2)]; %hyp2 freq components
    x1=1020; x2=1230; y1=46.8892; y2=90.9671; %coords of hyp3 rect ROIs
    hyp3t = [find(times2save==x1):find(times2save==x2)]; %hyp3 timepoints
    hyp3f = [dsearchn(frex',y1):dsearchn(frex',y2)]; %hyp3 freq components
else %mB
    load('mBgoodR1_dBbasenormpow.mat'); %monkey B data
    x1=580; x2=1230; y1=4.3972; y2=12.4581; %coords of hyp1 rect ROI
    hyp1t = [find(times2save==x1):find(times2save==x2)]; %hyp1 timepoints
    hyp1f = [dsearchn(frex',y1):dsearchn(frex',y2)]; %hyp1 freq components
    x1=580; x2=1230; y1=15.0551; y2=32.1077; %coords of hyp2 rect ROI
    hyp2t = [find(times2save==x1):find(times2save==x2)]; %hyp2 timepoints
    hyp2f = [dsearchn(frex',y1):dsearchn(frex',y2)]; %hyp2 freq components
    x1=1020; x2=1230; y1=38.8008; y2=90.9671; %coords of hyp3 rect ROIs
    hyp3t = [find(times2save==x1):find(times2save==x2)]; %hyp3 timepoints
    hyp3f = [dsearchn(frex',y1):dsearchn(frex',y2)]; %hyp3 freq components
end

%hyp1: f & p regions show increased low-band LFP power delay
%hyp2: f & p regions show decreased mid-band LFP power delay 
%hyp3: f region high-band LFP power bursts near end of delay

size(dbn_pow.dbnPower) % [resp(2) chant(143A/318B) frexidx(35) times2saveidx(187)]

allchans = size(dbn_pow.dbnPower,2); %total # of chans





mhyp1 = zeros(allchans,1)'; mhyp2 = zeros(allchans,1)'; mhyp3 = zeros(allchans,1)';
for chanN=1:allchans
    %avg over time then frequencies
    mhyp1(chanN)=mean( mean( squeeze( dbn_pow.dbnPower(1,chanN,hyp1f,hyp1t ) ),2 ),1 ); %correct
    mhyp2(chanN)=mean( mean( squeeze( dbn_pow.dbnPower(1,chanN,hyp2f,hyp2t ) ),2 ),1 ); %correct
end

%hypothesis 3 needs to be tested with only frontal channels
mhyp3(chanN)=mean( mean( squeeze( dbn_pow.dbnPower(1,chanN,hyp3f,hyp3t ) ),2 ),1 ); %correct;

squeeze(dbn_pow.dbnPower(1,chanN,hyp3f,hyp3t)); %testing
test=squeeze(dbn_pow.dbnPower(1,chanN,:,:)); %testing
mean(mhyp1)

[h1,p1,ci1,stats1] = ttest(mhyp1',0,'Tail','right','Alpha',pval); %hyp1 is sig for both monkeys
[h2,p2,ci2,stats2] = ttest(mhyp2',0,'Tail','left','Alpha',pval); %hyp2 is sig for both monkeys
[h3,p3,ci3,stats3] = ttest(mhyp3',0,'Tail','right','Alpha',pval); %hyp3 is not sig for both monkeys

%save stats for all 3 hypotheses as m#ttest_hyp1-3.mat in:
% OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\step 3 - stats


%finish writing up methods for stats and re-check code sections to ensure
%everything is included (see about adding fig showing correct vs db-norm'd)
%write up results section and include the figure
%write up discussion and conclusion

%% statistics via permutation testing

%init vars
srate = 1000; %1000Hz sampling rate
num_frex = 35; %50 for 200Hz, 35 for 100Hz, better for statistics mult comp corr, less smooth spectrogram
% define trial timeline
signalt = -.504:1/srate:1.589; %504 (nonzero sample) + 811 (delay) + 274 (match)=1589ms
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1466; % in ms, 1466 = 505 (sample) + 811 (delay) + 150 (match)
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');


% p-value, p<0.05 two-tailed is 0.025
pval = 0.025;

% convert p-value to Z value
zval = abs(norminv(pval));

% number of permutations
n_permutes = 1000;
n_mpermutes = 4; %meta-permutations
monkeys = {'mA','mB'};

mi = 1; %choose which monkey: mA=1, mB=2
monkey=monkeys{mi}; 
if monkey=='mA'
    load('mAgoodR1_pow_erp_lfp.mat'); %monkey A data
    mData = mAgoodR1; clear mAgoodR1
else
    load('mBgoodR1_pow_erp_lfp.mat'); %monkey B data
    mData = mBgoodR1; clear mBgoodR1
end
    
ppc = parallel.pool.Constant(mData);

tic
% meta-permutation test
parfor permN = 7:n_mpermutes
    mperm = sprintf('%sp%d.mat',monkey,permN);
    metaperm = permmapper(mData,n_permutes,num_frex,times2save);
    parsave(mperm,metaperm)    
end
toc

find(mAchan_diffpermmaps(5,:,:,:)~=0); %testing

permN=1; mperm = sprintf('p%d',permN); %choose permutation
%visualize TF point from diffmap over permutations histogram to see if it's a gaussian
figure(1), clf
% subplot(211)
total_chans = size(mA_metaperm.(mperm).diffmap,1); rand_chan = randperm(total_chans,1);
rand_frex = randperm(length(frex),1); rand_time = randperm(length(times2save),1);
histogram(mAchan_diffpermmaps(rand_chan,:,rand_frex,rand_time));
sprintf('Monkey A Day %d Chan %d Freq %fHz %dms Timepoint Cor-Inc Diffmap %d Permutations',...
    find(ismember(mA_metaperm.(mperm).day,mA_metaperm.(mperm).day(rand_chan))),...

%can't easily get the channel # from this structure, think about how to add
%in the chan# and total chans into permutation structure
allchans = size(mAgoodR1(1).(dayN{:}).power,1); %total # of chans

figure(1), clf
% subplot(211)
total_chans = size(mAchan_diffpermmaps,1); rand_chan = randperm(total_chans,1);
histogram(mAchan_diffpermmaps(rand_chan,:,randperm(length(frex),1),...

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
title(sprintf('%s Day %d / %d All %d Chans Across %d Areas Correct',monkeys{i}(1:8),...
    find(ismember(alldays,dayN{:})),length(alldays),allchans,length(unique(areas))));
ax=gca; ax.FontSize = 25;

%size(mAchan_diffpermmaps) %[143 1000 35 187] [allchans perm freqidx time]

%fix variables to reflect monkey-resp-chan-perm_diffmaps
%clean up data folder, tag files with section headers in this code
%plot test perm diffmap to see if its a gaussian

% mike x cohen convo about meta perm testing: 
% https://groups.google.com/g/analyzingneuraltimeseriesdata/c/5j-ej09p1DI/m/boVUBFaHAwAJ

