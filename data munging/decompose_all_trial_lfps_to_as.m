%% Step 1: Define Parameters for Wavelets

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
% export_fig('fwhm smoothing','-png','-transparent'); %save transparent pdf in pwd

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

%% Step 3: Extract all LFPs from both monkeys, all trials and save as separate file in each session
% First, ensure extract_and_save_LFPs.m is followed to get uniform epoch-based LFPs
% epochs: 
%   504 samples in baseline
%   505 samples in cue
%   811 samples in delay
%   259 samples in match (check this)
%   2079 total samples across all chans

% init vars, run this twice (once for each monkey):
%   homepc:
path = 'G:\\monkey_data\\';
%   labpc:
%   path = 'C:\\Users\\bryan\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\'
monk = 2; %1 = clark, 2 = betty
days_betty = { '090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001' };
days_clark = { '060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221' };
if monk==1
    monkey='clark';
    alldays = days_clark;
    days = { days_clark, "session02", "session03" };
else
    monkey='betty';
    alldays = days_betty;
    days = { days_betty, "session01" };
end
trial_info_path = strcat(path,'%s\\%s\\%s\\trial_info.mat'); %build trial_info path
recording_info_path = strcat(path,'%s\\%s\\%s\\recording_info.mat'); %build recording_info path
lfp_path = strcat(path,'%s\\%s\\%s\\%s%s%s.all.mat'); %build lfp raw data path

% begin definining convolution parameters
n_wavelet = length(wavet);
half_of_wavelet_size = floor(n_wavelet/2)+1;
% define trial timeline
signalt = -.504:1/srate:1.574; %504 (nonzero sample) + 811 (delay) + 259 (match)=1589ms (this is in seconds)
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1466; % in ms, 1466 = 505 (sample) + 811 (delay) + 150 (match)
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');
% % define baseline time
% baset = [-.4 -.1]; % in seconds
% baseidx = dsearchn(signalt',baset');
tic
for day=alldays %cycle through all days
    for j=2:3
        if (j==3) && (monkey=="betty") %only one session for betty
            continue %skip rest of loop
        end
        trial_infoN = sprintf(trial_info_path,monkey,day{:},days{j}); %create full path to trial_info.mat
        load(trial_infoN,'trial_info'); %load trial_info for day's trials
        recording_infoN = sprintf(recording_info_path,monkey,day{:},days{j}); %create full path to recording_info.mat
        load(recording_infoN,'recording_info'); %load recording_info for day's trials
        areas = recording_info.area;
        trial_lfp = sprintf(lfp_path,monkey,day{:},days{j},monkey,day{:},days{j}{1}(8:9));
        %load all trials across all chans in a single variable 'lfp':
        load(trial_lfp,'lfp'); %chan (recording_info.numChannels) x time (2049) x trial (trial_info.numTrials)
        for k=1:recording_info.numChannels %parse all channels                       
            signal = squeeze(lfp(k,:,:)); %rem single chan dim and confirm time x trials leftover
            reflectsig_all = zeros(size(signal,1)+2*n_wavelet,size(signal,2)); %initialize reflected signals mat
            % reflect all trials
            for signalN=1:size(signal,2) %loop through trials
                reflectsig = [ signal(n_wavelet:-1:1,signalN); signal(:,signalN); signal(end:-1:end-n_wavelet+1,signalN); ];        
                reflectsig_all(:,signalN) = reflectsig;
            end
            % concatenate into a super-trial
            reflectsig_supertri = reshape(reflectsig_all,1,[]); % reshape to 1D time-trials
%             % confirm this reflection actually works by visualization:
%             % plot original signal (example from mB chan 23 trial 500
%             figure(1), clf
%             subplot(211)
%             plot(signal(:,500)','LineWidth',2,'color','b')
%             set(gca,'xlim',[0 numel(reflectsig)]-n_wavelet)
%             ylabel('Voltage (\muV)'); title('Original Signal')
%             ax=gca; ax.FontSize = 25; x1=xticklabels; 
%             set(gca, 'XTickLabel', []); box off
%             % plot reflected signal
%             subplot(212)
%             p21=plot(n_wavelet+1:length(signal(:,500))+n_wavelet,signal(:,500),...
%                 'LineWidth',3,'color','b');
%             hold on
%             p22=plot(reflectsig,'-','LineWidth',2,'color','k');
%             p22.Color(4) = 0.35; %change transparency
%             set(gca,'xlim',[0 numel(reflectsig)])
%             title('Reflected Signal'); xlabel('Time step (ms)')
%             ylabel('Voltage (\muV)')
%             ax=gca; ax.FontSize = 25; xticklabels(x1); box off
%             legend({'original';'reflected'},'FontSize',25,'Location','best','box','off')
%             export_fig('reflected signal','-png','-transparent'); %save transparent pdf in pwd
            % step 1: finish defining convolution parameters
            n_data = length(reflectsig_supertri); % time*trials
            n_convolution = n_wavelet+n_data-1;
            % step 2: take FFTs
            fft_data = fft(reflectsig_supertri,n_convolution); % all trials for chan        
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
                ansig(k,fi,:,:) = as;
            end
        end
        lfp_as_path = strcat(path,'%s\\%s\\%s\\%s%s%s.ansig.mat'); %build lfp all analytic signals path
        lfp_as = sprintf(lfp_as_path,monkey,day{:},days{j},monkey,day{:},days{j}{1}(8:9));
        save(lfp_as,'ansig','-v7.3') %save chan x freq (35) x time (2079) x trials ansig matrix
        clear ansig
%         % testing:
%         find(ansig(:,:,:,:)==0);
%         size(ansig)
    end
end %all trial lfp's are now saved in a single file per session
toc