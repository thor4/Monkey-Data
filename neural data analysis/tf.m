%% Step 1: Get Signal

%initialize variables
srate = 1000; % 1,000Hz
% -500:0 baseline, 1:500 sample, 501:1310 delay (all in ms)
wavet = -.5:1/srate:.5; % in seconds 
min_freq = 3.5;
max_freq = 200;
num_frex = 50;
min_fwhm = .350; % in seconds
max_fwhm = .050; % in seconds
wavpts = length(wavet);
hz = linspace(0,srate/2,floor(wavpts/2)+1); % pos frequencies (not neg) up to Nyquist
% hzf = linspace(0,srate,N); % total frequencies from fft
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
% s    = logspace(log10(3),log10(12),num_frex)./(2*pi*frex); %width of gaussian
fwhm = logspace(log10(min_fwhm),log10(max_fwhm),length(frex)); % in seconds **TRY SPECTROGRAMS WITH 25ms and with 50ms**

% %% Make wavelets, ensure they taper to 0 at either end in time domain and
% % ensure Gaussian in frequency domain
% 
% % Using n-cycles (Cohen's code from Google group)
% srate = 1000;
% frex  = linspace(.5,200,50);
% s2 = 2*(linspace(5,25,length(frex))./(2*pi*frex)).^2;
% 
% wavet = -.5:1/srate:.5;
% 
% for fi=1:length(frex)
%     waves(fi,:) = exp(2*1i*pi*frex(fi)*wavet).*exp(-wavet.^2/s2(fi));
% end
% 
% subplot(211), plot(wavet,real(waves))
% 
% % fft all the wavelets
% for ffti=1:size(waves,1)
%     waves_fft(ffti,:) = fft(waves(ffti,:));
% end
% 
% subplot(212)
% plot(linspace(0,srate,length(wavet)),abs(waves_fft).^2)
% set(gca,'xlim',[30 62])
% % now plot wavelets in groups of 10 at a time..
% figure
% for fi=1:5
%     subplot(5,1,fi)
%     plot(wavet,real(wavelet(fi,:)))
%     title(sprintf('Wavelet at %fHz, s=%f number of cycles',frex(fi),s(fi)*(2*pi*frex(fi))));
% end

%% Make wavelets, ensure they taper to 0 at either end in time domain and
% ensure Gaussian in frequency domain

% now using fwhm-specified in time domain
midp = dsearchn(wavet',0);
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
set(gca,'xlim',[0 225])
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

%% Create analytic signal using the complex morlet wavelet family for both 
% monkeys, both conditions, all days, all channels and all trials

% load monkey data
load('mGoodStableRule1PingRej-split_by_Day_BehResp_and_Chan.mat')

% initialize variables
signalt = -.5:1/srate:1.31; % in seconds
% vector of time points to save in post-analysis downsampling
times2save = -500:10:1310; % in ms
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');

% setup response structs
% monkey(1).correct=[];
% monkey(1).incorrect=[];
% begin definining convolution parameters
n_wavelet = length(wavet);
half_of_wavelet_size = floor(n_wavelet/2)+1;
monkeyN = 1; % which monkey (1 or 2)
m1chans = {'8B', '9L', 'dPFC', 'vPFC', 'LIP', 'MIP', 'PEC', 'PG'};
m1chansa = {'a8B', 'a9L', 'adPFC', 'avPFC', 'aLIP', 'aMIP', 'aPEC', 'aPG'};

% get signal
% d = 'd%i';
% numel(monkey(monkeyN).day)
tic
%correct
for i=1:numel(monkey(monkeyN).day)
    chan = fieldnames(monkey(monkeyN).day(i).correct);
%     day = sprintf(d, i);
    for j=1:numel(chan)
        signal = monkey(monkeyN).day(i).correct.(chan{j})'; % change to time-by-trials
        reflectsig_all = zeros(size(signal,1)+2*n_wavelet,size(signal,2)); %initialize reflected signals mat
        % reflect all signals
        for signalN=1:size(signal,2)
            reflectsig = [ signal(n_wavelet:-1:1,signalN); signal(:,signalN); signal(end:-1:end-n_wavelet+1,signalN); ];        
            reflectsig_all(:,signalN) = reflectsig;
        end
        % concatenate into a super-trial
        reflectsig_supertri = reshape(reflectsig_all,1,[]); % reshape to 1D time-trials
%         % plot original signal
%         figure(1), clf
%         subplot(211)
%         plot(signal(:,1)','LineWidth',2)
%         set(gca,'xlim',[0 numel(reflectsig)]-n_wavelet)
%         ylabel('Voltage (\muV)','FontSize',14)
%         title('Original Signal','FontSize',16)
%         % plot reflected signal
%         subplot(212)
%         plot(n_wavelet+1:length(signal(:,1))+n_wavelet,signal(:,1),'LineWidth',5,'color','g')
%         hold on
%         plot(reflectsig,'LineWidth',2)
%         set(gca,'xlim',[0 numel(reflectsig)])
%         title('Reflected Signal')
%         xlabel('Time step (ms)','FontSize',14), ylabel('Voltage (\muV)','FontSize',14)
%         title('Reflected Signal','FontSize',16)
%         legend({'original';'reflected'},'FontSize',14,'Location','best')
        % step 1: finish defining convolution parameters
        n_data = length(reflectsig_supertri); % time*trials
        n_convolution = n_wavelet+n_data-1;
        % step 2: take FFTs
        fft_data = fft(reflectsig_supertri,n_convolution); % all trials for chan
        for k=1:numel(m1chans)
            if endsWith(chan{j},m1chans{k}) %which chan
                for fi=1:length(frex)
                    % FFT of wavelet
                    fft_wavelet = fft(wavelets(fi,:),n_convolution);
                    % step 3: normalize kernel by scaling amplitudes to one in the 
                    % frequency domain. prevents amplitude from decreasing with 
                    % increasing frequency. diff from 1/f scaling
                    fft_wavelet = fft_wavelet ./ max(fft_wavelet);
                    % step 4: point-wise multiply and take iFFT
                    as_ = ifft( fft_data.*fft_wavelet ); % analytic signal
                    % step 5: trim wings
                    as_ = as_(half_of_wavelet_size:end-half_of_wavelet_size+1);
                    % step 6: reshape back to reflected time-by-trials
                    as_ = reshape(as_,size(reflectsig_all,1),size(reflectsig_all,2));
                    % step 7: chop off the reflections
                    as_ = as_(n_wavelet+1:end-n_wavelet,:);
                    as(fi,:,:) = abs( as_ ).^2; % store mean power for each frequency
                    clear fft_wavelet as_ % start anew with these var's ea. loop
                end
                if isfield(monkey(monkeyN).correct,(m1chansa{k})) % field exists?
                    aTrials = size(monkey(monkeyN).correct.(m1chansa{k}),3);
                    % append trials for area
                    monkey(monkeyN).correct.(m1chansa{k})(:,:,aTrials+1:aTrials+size(signal,2)) = as;
                else
                    monkey(monkeyN).correct.(m1chansa{k}) = as; % freq x time x trials
                end
                clear as % start anew with this var ea. loop
            end
        end
            % store analytic signal for all frequencies for each channel
%         monkey(monkeyN).day(i).correct.(chan{j}) = as;
%         clear as % start anew with this var for ea. loop
    end
end
toc

tic
%incorrect
for i=1:numel(monkey(monkeyN).day)
    chan = fieldnames(monkey(monkeyN).day(i).incorrect);
%     day = sprintf(d, i);
    for j=1:numel(chan)
        signal = monkey(monkeyN).day(i).incorrect.(chan{j})'; % change to time-by-trials
        signal_alltrials = reshape(signal,1,[]); % reshape to 1D time-trials
        % step 1: finish defining convolution parameters
        n_data = length(signal_alltrials); % time*trials
        n_convolution = n_wavelet+n_data-1;
        % step 2: take FFTs
        fft_data = fft(signal_alltrials,n_convolution); % all trials for chan
        for k=1:numel(m1chans)
            [x]=endsWith(chan{j},m1chans);
            if endsWith(chan{j},m1chans{k}) %which chan
                for fi=1:length(frex)
                    % FFT of wavelet
                    fft_wavelet = fft(wavelets(fi,:),n_convolution);
                    % step 3: normalize kernel by scaling amplitudes to one in the 
                    % frequency domain. prevents amplitude from decreasing with 
                    % increasing frequency. diff from 1/f scaling
                    fft_wavelet = fft_wavelet ./ max(fft_wavelet);
                    % step 4: point-wise multiply and take iFFT
                    as_ = ifft( fft_data.*fft_wavelet ); % analytic signal
                    % step 5: trim wings
                    as_ = as_(half_of_wavelet_size:end-half_of_wavelet_size+1);
                    % step 6: reshape back to time-by-trials
                    as_ = reshape(as_,size(signal,1),size(signal,2));
                    as(k,fi,:) = mean( abs( as_ ).^2, 2 ); % store trial-averaged power for each frequency
                    clear fft_wavelet as_ % start anew with these var's ea. loop
                end
                % how to store as sessions? endsWith skips earlier values
                % of k if parietal area comes first in chan1 ie
                if isfield(monkey(monkeyN).incorrect,(m1chansa{k})) % field exists?
                    aTrials = size(monkey(monkeyN).incorrect.(m1chansa{k}),3);
                    monkey(monkeyN).incorrect.(m1chansa{k})(:,:,aTrials+1:aTrials+size(signal,2)) = as;
                else
                    monkey(monkeyN).incorrect.(m1chansa{k}) = as; % freq x time x trials
                end
                clear as % start anew with this var ea. loop
            end
        end
            % store analytic signal for all frequencies for each channel
%         monkey(monkeyN).day(i).correct.(chan{j}) = as;
%         clear as % start anew with this var for ea. loop
    end
end
toc

%% dB-normalize the session data


%% Aggregate power for each area

% load data
% monkey 1 all days correct + incorrect
tic
load('m1GoodStableRule1PingRejAnSig-AllDays_Cor_Inc_Allchans.mat')
toc

monkeyN = 1; %which monkey
m1chans = {'8B', '9L', 'dPFC', 'vPFC', 'LIP', 'MIP', 'PEC', 'PG'};
m1chansa = {'a8B', 'a9L', 'adPFC', 'avPFC', 'aLIP', 'aMIP', 'aPEC', 'aPG'};
% setup matrix as area x condition x frequency x time x trial using fake
% data so that cat function below will work
monkey18B = repmat(1i,[8 2 50 1811 2]);
% % remove extra rows created during initialization
% monkey18B(2,:,:,:,:)=[]; monkey18B(:,2,:,:,:)=[]; monkey18B(:,:,2,:,:)=[]; 
% monkey18B(:,:,:,2,:)=[]; monkey18B(:,:,:,:,2)=[]; monkey18B(1,:,:,:,:)=[]; 
% monkey18B(:,1,:,:,:)=[]; monkey18B(:,:,1,:,:)=[]; monkey18B(:,:,:,1,:)=[]; 
% monkey18B(:,:,:,:,1)=[]; 

% initialize signal variable
signal = zeros(8,2,50,1811,2);

% for k=1:numel(m1chans) % initialize variables
%     m1.correct.(m1chansa{k})=ones(1811,1,50);
%     m1.incorrect.(m1chansa{k})=ones(1811,1,50);
% end

tic
% correct
for i=1:numel(monkey1.day)
    chan = fieldnames(monkey1.day(i).correct);
    for j=1:numel(chan)
        k=1; % which chan
        if endsWith(chan{j},m1chans{k}) %which chan
            % change to area x cond x freq x time x trial
            data = permute(monkey1.day(i).correct.(chan{j}),[3 1 2]);
            signal = zeros(8,2,50,1811,size(data,3));
            signal(7,1,:,:,:) = permute(monkey1.day(i).correct.(chan{j}),[3 1 2]);
            % save in aggregate data structure
            monkey18B = cat(5,monkey18B,signal);
            m1.correct.(m1chansa{k}) = cat(2,m1.correct.(m1chansa{k}),);
        end
        end
    end
    monkey1.day(i).correct = []; %remove day to free-up memory
end
signal = monkey1.day(1).correct.chan1PEC;
signal_ = permute(signal,[3 2 1]);
signal(1,1,45)
signal_(45,1,1)
% incorrect
for i=1:numel(monkey1.day)
    chan = fieldnames(monkey1.day(i).incorrect);
    for j=1:numel(chan)
        for k=1:numel(m1chans) % initialize variables
            if endsWith(chan{j},m1chans{k}) %which chan
                % save in aggregate data structure
                m1.incorrect.(m1chansa{k}) = cat(2,m1.incorrect.(m1chansa{k}),monkey1.day(i).incorrect.(chan{j}));
            end
        end
    end
    monkey1.day(i).incorrect = []; %remove day to free-up memory
end
toc

% initialize output time-frequency data
% notice that here we save all trials
tf = zeros(2,num_frex,length(times2save),ntrials);

for fi=1:num_frex
    
    % run convolution
    as1 = ifft(cmwX(fi,:).*dataX1);
    as1 = as1(half_wav+1:end-half_wav);
    as1 = reshape(as1,npnts,ntrials);
    
    % power on all trials from channel "1"
    % only from times2saveidx!
    tf(1,fi,:,:) = abs( as1(times2saveidx,:) ).^2;
    
    
    % run convolution
    as2 = ifft(cmwX(fi,:).*dataX2);
    as2 = as2(half_wav+1:end-half_wav);
    as2 = reshape(as2,npnts,ntrials);
    
    % power on all trials from channel "2"
    % only from times2saveidx!
    tf(2,fi,:,:) = 
end

% for convenience, compute the difference in power between the two channels
diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));


%% plotting the raw data

clim = [0 20000];

figure(1), clf
subplot(221)
imagesc(times2save,frex,squeeze(mean( tf(1,:,:,:),4 )))
set(gca,'clim',clim,'ydir','n','xlim',xlim)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Channel ' num2str(chan1idx) ])

subplot(222)
imagesc(times2save,frex,squeeze(mean( tf(2,:,:,:),4 )))
set(gca,'clim',clim,'ydir','n','xlim',xlim)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Channel ' num2str(chan2idx) ])

subplot(223)
imagesc(times2save,frex,diffmap)
set(gca,'clim',[-mean(clim) mean(clim)],'ydir','n','xlim',xlim)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Difference: channels ' num2str(chan2idx) ' - ' num2str(chan1idx) ])

%% statistics via permutation testing

% p-value
pval = 0.05;

% convert p-value to Z value
zval = abs(norminv(pval));

% number of permutations
n_permutes = 1000;

% initialize null hypothesis maps
permmaps = zeros(n_permutes,num_frex,length(times2save));

% for convenience, tf power maps are concatenated
%   in this matrix, trials 1:ntrials are from channel "1" 
%   and trials ntrials+1:end are from channel "2"
tf3d = cat(3,squeeze(tf(1,:,:,:)),squeeze(tf(2,:,:,:)));


% generate maps under the null hypothesis
for permi = 1:n_permutes
    
    % randomize trials, which also randomly assigns trials to channels
    randorder = randperm(size(tf3d,3));
    temp_tf3d = tf3d(:,:,randorder);
    
    % compute the "difference" map
    % what is the difference under the null hypothesis?
    permmaps(permi,:,:) = squeeze( mean(temp_tf3d(:,:,1:ntrials),3) - mean(temp_tf3d(:,:,ntrials+1:end),3) );
end

%% show non-corrected thresholded maps

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps));
std_h0  = squeeze(std(permmaps));

% now threshold real data...
% first Z-score
zmap = (diffmap-mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;


% now some plotting...

figure(2), clf

subplot(221)
imagesc(times2save,frex,diffmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','nor')
title('TF map of real power values')

subplot(222)
imagesc(times2save,frex,diffmap);
hold on
contour(times2save,frex,logical(zmap),1,'linecolor','k');
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')
title('Power values and outlined significance regions')

subplot(223)
imagesc(times2save,frex,zmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-10 10],'xlim',xlim,'ydir','no')
title('Thresholded TF map of Z-values')

%% 






%% corrections for multiple comparisons

% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);
% ... and for maximum-pixel based correction
max_val = zeros(n_permutes,2); % "2" for min/max

% loop through permutations
for permi = 1:n_permutes
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(permi,:,:));
    threshimg = (threshimg-mean_h0)./std_h0;
    
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    
    
    % find clusters (need image processing toolbox for this!)
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
    end
    
    
    % get extreme values (smallest and largest)
    temp = sort( reshape(permmaps(permi,:,:),1,[] ));
    max_val(permi,:) = [ min(temp) max(temp) ];
    
end

%% show histograph of maximum cluster sizes

figure(3), clf
hist(max_cluster_sizes,20);
xlabel('Maximum cluster sizes'), ylabel('Number of observations')
title('Expected cluster sizes under the null hypothesis')


% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));

%% plots with multiple comparisons corrections

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<cluster_thresh
        zmap(islands.PixelIdxList{i})=0;
    end
end

% plot tresholded results
figure(4), clf
subplot(221)
imagesc(times2save,frex,diffmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power, no thresholding') 
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')

subplot(222)
imagesc(times2save,frex,diffmap)
hold on
contour(times2save,frex,logical(zmap),1,'linecolor','k')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power with contour')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')

subplot(223)
imagesc(times2save,frex,zmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('z-map, thresholded')
set(gca,'clim',[-13 13],'xlim',xlim,'ydir','normal')

%% now with max-pixel-based thresholding

% find the threshold for lower and upper values
thresh_lo = prctile(max_val(:,1),100-100*pval); % what is the
thresh_hi = prctile(max_val(:,2),100-100*pval); % true p-value?

% threshold real data
zmap = diffmap;
zmap(zmap>thresh_lo & zmap<thresh_hi) = 0;

figure(5), clf
subplot(221)
imagesc(times2save,frex,diffmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('tf power map, no thresholding') 
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','n')

subplot(222)
imagesc(times2save,frex,diffmap)
hold on
contour(times2save,frex,logical(zmap),1,'linecolor','k')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('tf power map with contour')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','normal')

subplot(223)
imagesc(times2save,frex,zmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('tf power map, thresholded')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','no')

%% end.

%% Visualizations

figure(8), clf

subplot(211)
plot(signalt,mean(signal,2))
% set(gca,'ylim',[-.05 1.05])
ylabel('Amplitude'), title('ERP')

subplot(212)
contourf(sigtime,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 .1])
xlabel('Time (s)'), ylabel('Frequency (Hz)'), 
% ylabel('Raw power (\muV^2)') colorbar label


% plot for comparison
figure
subplot(311)
plot(EEG.times,real(convolution_result_fft))
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
title([ 'Projection onto real axis is filtered signal at ' num2str(frequency) ' Hz.' ])

subplot(312)
plot(EEG.times,abs(convolution_result_fft).^2)
xlabel('Time (ms)'), ylabel('Power (\muV^2)')
title([ 'Magnitude of projection vector squared is power at ' num2str(frequency) ' Hz.' ])

subplot(313)
plot(EEG.times,angle(convolution_result_fft))
xlabel('Time (ms)'), ylabel('Phase angle (rad.)')
title([ 'Angle of vector is phase angle time series at ' num2str(frequency) ' Hz.' ])


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

