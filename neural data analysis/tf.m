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
plot(wavet,real(wavelets))
axis off
export_fig test2.png -transparent % no background
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

% begin definining convolution parameters
n_wavelet = length(wavet);
half_of_wavelet_size = floor(n_wavelet/2)+1;
monkeyN = 2; % which monkey (1 or 2)
m1chans = {'8B', '9L', 'dPFC', 'vPFC', 'LIP', 'MIP', 'PEC', 'PG'};
m1areas = {'a8B', 'a9L', 'adPFC', 'avPFC', 'aLIP', 'aMIP', 'aPEC', 'aPG'};
m2chans = {'6DR', '8AD', '8B', 'dPFC', 'LIP', 'PE', 'PEC', 'PG'};
m2areas = {'a6DR', 'a8AD', 'a8B', 'adPFC', 'aLIP', 'aPE', 'aPEC', 'aPG'};

% define trial timeline
signalt = -.5:1/srate:1.31; % in seconds
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1200; % in ms
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');
% define baseline time
baset = [-.4 -.1]; % in seconds
baseidx = dsearchn(signalt',baset');

% setup response structs
monkey(monkeyN).correct=[];
monkey(monkeyN).incorrect=[];

tic
%correct
for i=1:numel(monkey(monkeyN).day)
    chan = fieldnames(monkey(monkeyN).day(i).correct);
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
        % which area is this chan
        for k=1:numel(m2areas)
            if endsWith(chan{j},m2chans{k}) %which chan
                area=k;
            end
        end
        for fi=1:length(frex)
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
            % save raw power for ea. freq x down-sampled time x trials 
            pow(fi,:,:) = abs( as(times2saveidx,:) ) .^2;
            % mean( abs( as_ ).^2, 2);
            clear fft_wavelet as % start anew with these var's ea. loop
        end
        if isfield(monkey(monkeyN).correct,(m2areas{area})) % field exists?
            % yes, field exists..
            aTrials = size(monkey(monkeyN).correct.(m2areas{area}),3);
            % append session for area in freq x time x trial struct
            monkey(monkeyN).correct.(m2areas{area})(:,:,aTrials+1:aTrials+size(pow,3)) = pow;
        else
            monkey(monkeyN).correct.(m2areas{area}) = pow; % freq x time x trial
        end
        clear pow % start anew with this var ea. loop
    end
end

%incorrect
for i=1:numel(monkey(monkeyN).day)
    chan = fieldnames(monkey(monkeyN).day(i).incorrect);
    for j=1:numel(chan)
        signal = monkey(monkeyN).day(i).incorrect.(chan{j})'; % change to time-by-trials
        reflectsig_all = zeros(size(signal,1)+2*n_wavelet,size(signal,2)); %initialize reflected signals mat
        % reflect all signals
        for signalN=1:size(signal,2)
            reflectsig = [ signal(n_wavelet:-1:1,signalN); signal(:,signalN); signal(end:-1:end-n_wavelet+1,signalN); ];        
            reflectsig_all(:,signalN) = reflectsig;
        end
        % concatenate into a super-trial
        reflectsig_supertri = reshape(reflectsig_all,1,[]); % reshape to 1D time-trials
        % step 1: finish defining convolution parameters
        n_data = length(reflectsig_supertri); % time*trials
        n_convolution = n_wavelet+n_data-1;
        % step 2: take FFTs
        fft_data = fft(reflectsig_supertri,n_convolution); % all trials for chan
        % which area is this chan
        for k=1:numel(m2areas)
            if endsWith(chan{j},m2chans{k}) %which chan
                area=k;
            end
        end
        for fi=1:length(frex)
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
            % as_ is now a time x trial complex matrix
%             % compute baseline power averaged over trials then timepoints
%             basePow = mean( mean( abs( as(baseidx(1):baseidx(2),:).^2 ),2 ),1 );
            % store dB-norm'd down-sampled power for each frequency in freq
            % x time x trials
%             dbpow(fi,:,:) = 10*log10( abs( as(times2saveidx,:) ) .^2 ./basePow); 
            % save raw power for ea. freq x down-sampled time x trials 
            pow(fi,:,:) = abs( as(times2saveidx,:) ) .^2;
            % mean( abs( as_ ).^2, 2);
            clear fft_wavelet as % start anew with these var's ea. loop
        end
        if isfield(monkey(monkeyN).incorrect,(m2areas{area})) % field exists?
            % yes, field exists..
            aTrials = size(monkey(monkeyN).incorrect.(m2areas{area}),3);
            % append session for area in freq x time x trial struct
            monkey(monkeyN).incorrect.(m2areas{area})(:,:,aTrials+1:aTrials+size(pow,3)) = pow;
        else
            monkey(monkeyN).incorrect.(m2areas{area}) = pow; % freq x time x trial
        end
        clear pow % start anew with this var ea. loop
    end
end

toc

% %% dB-normalize the session data
% 
% % load raw mean power monkey data
% load('mGoodStableRule1PingRejRawMeanPow-AllDays_Cor_Inc_Allchans.mat')
% 
% % re-define areas if skipped analytic signal area
% m1areas = {'a8B', 'a9L', 'adPFC', 'avPFC', 'aLIP', 'aMIP', 'aPEC', 'aPG'};
% m2areas = {'a6DR', 'a8AD', 'a8B', 'adPFC', 'aLIP', 'aPE', 'aPEC', 'aPG'};
% 
% % define trial timeline
% signalt = -.5:1/srate:1.31; % in seconds
% % vector of time points to save in post-analysis downsampling
% times2save = -500:10:1310; % in ms
% % time vector converted to indices
% times2saveidx = dsearchn((signalt.*1000)',times2save');
% % define baseline time
% baset = [-.4 -.1]; % in seconds
% baseidx = dsearchn(signalt',baset');
% 
% % compute condition-averaged baseline for each area in each monkey
% % mean over sessions then baseline time period
% % stored as freq x time x session
% monkeyN=1; % monkey 1
% for areaN=1:numel(m1areas)
%     cor = mean( mean( monkey(monkeyN).correct.(m1areas{areaN})(:,baseidx(1):baseidx(2),:),3 ),2 );
%     inc = mean( mean( monkey(monkeyN).incorrect.(m1areas{areaN})(:,baseidx(1):baseidx(2),:),3 ),2 );
%     base = mean([cor inc],2); % condition-averaged (c-a) baseline
%     % dB-normalize the session data with condition-averaged baseline
%     raw_pow_cor = monkey(monkeyN).correct.(m1areas{areaN});
%     raw_pow_inc = monkey(monkeyN).correct.(m1areas{areaN});
%     % divide ea frequency's time x sesh by its respective c-a baseline,
%     % then db transform
%     db_base_pow_cor = 10*log10( raw_pow_cor./base ); 
%     db_base_pow_inc = 10*log10( raw_pow_inc./base ); 
%     monkey_(monkeyN).correct.(m1areas{areaN}) = db_base_pow_cor;
%     monkey_(monkeyN).incorrect.(m1areas{areaN}) = db_base_pow_inc ;
% end
% 
% monkeyN=2; % monkey 2
% for areaN=1:numel(m2areas)
%     cor = mean( mean( monkey(monkeyN).correct.(m2areas{areaN})(:,baseidx(1):baseidx(2),:),3 ),2 );
%     inc = mean( mean( monkey(monkeyN).incorrect.(m2areas{areaN})(:,baseidx(1):baseidx(2),:),3 ),2 );
%     base = mean([cor inc],2);
%     % dB-normalize the session data with condition-averaged baseline
%     raw_pow_cor = monkey(monkeyN).correct.(m2areas{areaN});
%     raw_pow_inc = monkey(monkeyN).correct.(m2areas{areaN});
%     % divide ea frequency's time x sesh by its respective c-a baseline,
%     % then db transform
%     db_base_pow_cor = 10*log10( raw_pow_cor./base ); 
%     db_base_pow_inc = 10*log10( raw_pow_inc./base ); 
%     monkey_(monkeyN).correct.(m2areas{areaN}) = db_base_pow_cor;
%     monkey_(monkeyN).incorrect.(m2areas{areaN}) = db_base_pow_inc ;
% end
% 
% % average db normalized power over sessions yielding freq x time
% rawpowC = mean( monkey_pow(monkeyN).correct.(m2areas{areaN}),3 );
% rawpowI = mean( monkey_pow(monkeyN).incorrect.(m2areas{areaN}),3 );
% 
% %% Visualize power
% 
% % load data
% % both monkeys dB-normalized power by area
% tic
% load('mGoodStableRule1PingRejdBPow-AllDays_Cor_Inc_Allchans.mat')
% toc
% % both monkeys raw power by area
% tic
% load('mGoodStableRule1PingRejRawPow-AllDays_Cor_Inc_Allchans.mat')
% toc
% 
% % contourf plot template:
% % x = 1 x samples
% % y = 1 x frequencies
% % z = frequencies x samples
% % contourf(x,y,z,...)
% monkeyN = 1;
% areaN = 1;
% 
% figure(6), clf
% subplot(221)
% contourf(signalt(times2saveidx),frex,mean( monkey(monkeyN).correct.(m1areas{areaN}),3 ),'linecolor','none')
% % ,'clim',[-3 3]
% set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
% ylabel('Frequency (Hz)')
% title(sprintf('Raw Power via Morlet Wavelet from Monkey %d, Area %s, Resp Correct',monkeyN,m1areas{areaN}));
% cbar = colorbar; set(get(cbar,'label'),'string','Raw power (\muV^2)');   
% subplot(223)
% contourf(signalt(times2saveidx),frex,mean( monkey(monkeyN).incorrect.(m1areas{areaN}),3 ),'linecolor','none')
% % contourf(signalt(times2saveidx),frex,pow(:,times2saveidx),'linecolor','none')
% set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
% % set(gca,'xlim',[monkeyTime(1) monkeyTime(end)])
% xlabel('Time (ms)'), ylabel('Frequency (Hz)')
% title(sprintf('Raw Power via Morlet Wavelet from Monkey %d, Area %s, Resp Incorrect',monkeyN,m1areas{areaN}));
% cbar = colorbar; set(get(cbar,'label'),'string','Raw power (\muV^2)');   
% subplot(222)
% contourf(signalt(times2saveidx),frex,mean( monkey(monkeyN).correct.(m1areas{areaN+5}),3 ),'linecolor','none')
% set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
% title(sprintf('Raw Power via Morlet Wavelet from Monkey %d, Area %s, Resp Correct',monkeyN,m1areas{areaN+5}));
% % cbar = colorbar; set(get(cbar,'label'),'string','dB change from baseline');   
% cbar = colorbar; set(get(cbar,'label'),'string','Raw power (\muV^2)');   
% subplot(224)
% contourf(signalt(times2saveidx),frex,mean( monkey(monkeyN).incorrect.(m1areas{areaN+5}),3 ),'linecolor','none')
% set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
% xlabel('Time (ms)')
% title(sprintf('Raw Power via Morlet Wavelet from Monkey %d, Area %s, Resp Incorrect',monkeyN,m1areas{areaN+5}));
% % title(sprintf('Raw Power via Morlet Wavelet from Monkey %d, Area %s, Resp Incorrect \n Freq Step = %d, temporal FWHM %dms to %.0dms, spectral FWHM %dHz - %dHz',monkeyN,m2areas{areaN},length(frex),empfwhmT(1)*1000,empfwhmT(end)*1000,empfwhmF(1),empfwhmF(end)));
% cbar = colorbar; set(get(cbar,'label'),'string','Raw power (\muV^2)');   
% suptitle(sprintf('Freq Step = %d, temporal FWHM %dms to %1.0fms, spectral FWHM %dHz to %dHz',length(frex),empfwhmT(1)*1000,empfwhmT(end)*1000,empfwhmF(1),empfwhmF(end)));

%% initialize variables

load('mGoodStableRule1PingRejRawPow-AllDays_Cor_Inc_Allchans.mat')

monkeyN = 2; %which monkey? 1 or 2
responses = [{'correct'},{'incorrect'}]; 
responseN = 1; %which response? 1 or 2

% define trial timeline
signalt = -.5:1/srate:1.31; % in seconds
times2save = -400:10:1200; % in ms
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');
% define baseline time
baset = [-400 -100]; % in ms
baseidx = dsearchn(times2save',baset'); % search down-sampled time

% which areas depends on monkeyN
areas = fieldnames(monkey(monkeyN).(responses{responseN}));
n_areas = length(areas);


%% plotting the raw difference data

areaN = 2; %plot which area?
clim = [0 90];

%pull out correct & incorrect power avg over trials
cor = mean( monkey(monkeyN).(responses{1}).(areas{areaN}),3 ); 
inc = mean( monkey(monkeyN).(responses{2}).(areas{areaN}),3 );
% compute the difference in power between the two conditions
% diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));
diffmap = cor - inc;

figure(9), clf
subplot(221)
imagesc(times2save,[],cor)
hold on %vertical line
yL = get(gca,'YLim'); line([0 500;0 500],yL,'Color','k','LineStyle',':'); 
set(gca,'clim',clim,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Raw Power'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1 pos(2)];
title(sprintf('Monkey %d, Area %s, Resp Correct',monkeyN,areas{areaN}(2:end)));

subplot(222)
imagesc(times2save,[],inc)
hold on %vertical line
yL = get(gca,'YLim'); line([0 500;0 500],yL,'Color','k','LineStyle',':'); 
set(gca,'clim',clim,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Raw Power'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1 pos(2)];
title(sprintf('Monkey %d, Area %s, Resp Incorrect',monkeyN,areas{areaN}(2:end)));

subplot(223)
imagesc(times2save,[], diffmap)
hold on %vertical line
yL = get(gca,'YLim'); line([0 500;0 500],yL,'Color','k','LineStyle',':'); 
set(gca,'clim',[-mean(clim)/5 mean(clim)/5],'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Raw Power'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1 pos(2)];
title(sprintf('Monkey %d, Area %s, Correct - Incorrect',monkeyN,areas{areaN}(2:end)));

%% statistics via permutation testing

% number of tests
n_tests = 8;

% p-value, p<0.05 two-tailed is 0.025
pval = 0.025/n_tests;

% convert p-value to Z value
zval = abs(norminv(pval));

% number of permutations
n_permutes = 1000;

tic
% meta-permutation test
for permN = 1:20
	monkeyN = 1; % first monkey 1
    % generate maps under the null hypothesis
    [m1_permmaps, ~] = permmapper(monkey,monkeyN,n_permutes,num_frex,times2save);
    % append permutations for monkey in chan x permutation x freq x time diffmap
    m1_meta_permmaps(:,(permN-1)*n_permutes+1:permN*n_permutes,:,:) = m1_permmaps;
    monkeyN = 2; % now monkey 2
    % generate maps under the null hypothesis
    [m2_permmaps, ~] = permmapper(monkey,monkeyN,n_permutes,num_frex,times2save);
    m2_meta_permmaps(:,(permN-1)*n_permutes+1:permN*n_permutes,:,:) = m2_permmaps;
end
toc

%% show non-corrected thresholded maps

% now some plotting...

figure(8), clf

% contourf plot template:
% x = 1 x samples
% y = 1 x frequencies
% z = frequencies x samples
% contourf(x,y,z,...)

figure(8), clf
subplot(221)
contourf(signalt(times2saveidx),frex,diffmap,'linecolor','none')
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
pos = get(cbar,'Position'); lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Power (\muV^2)'; cbar.Label.Position=[pos(1)+1 pos(2)];
cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('TF Map Monkey %d, Area %s, Correct > Incorrect',monkeyN,m2areas{areaN+3}(2:end)));

subplot(222)
contourf(signalt(times2saveidx),frex,diffmap,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap),1,'linecolor','k');
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
pos = get(cbar,'Position'); lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Power (\muV^2)'; cbar.Label.Position=[pos(1)+1 pos(2)];
cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('Significant Regions Outlined Monkey %d, Area %s, Correct > Incorrect',monkeyN,m2areas{areaN+3}(2:end)));

subplot(223)
contourf(signalt(times2saveidx),frex,zmap,'linecolor','none')
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
pos = get(cbar,'Position'); lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Power (\muV^2)'; cbar.Label.Position=[pos(1)+1 pos(2)];
cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('Thresholded Z-values Monkey %d, Area %s, Correct > Incorrect',monkeyN,m2areas{areaN+3}(2:end)));


%% corrections for multiple comparisons

% load permutation maps for all areas. 
% 4-d: area x permutation x freq x down-sampled time points
load('permmaps.mat')

monkeyN = 2; % which monkey: 1 or 2

if monkeyN==1
    areas = {'a8B', 'a9L', 'adPFC', 'avPFC', 'aLIP', 'aMIP', 'aPEC', 'aPG'}; %monkey1
    permmaps = m1_meta_permmaps;
else
    areas = {'a6DR', 'a8AD', 'a8B', 'adPFC', 'aLIP', 'aPE', 'aPEC', 'aPG'}; %monkey2
    permmaps = m2_meta_permmaps;
end

for areaN = 1:numel(areas)
    cor = mean( monkey(monkeyN).correct.(areas{areaN}),3 );
    inc = mean( monkey(monkeyN).incorrect.(areas{areaN}),3 );

    % for convenience, compute the difference in power between the two
    % conditions, correct-incorrect
    % diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));
    diffmap = cor - inc;
    % compute mean and standard deviation maps under the null hypothesis
    mean_h0 = squeeze(mean(permmaps(areaN,:,:,:)));
    std_h0  = squeeze(std(permmaps(areaN,:,:,:)));
    % now change data to zmap
    zmap = (diffmap-mean_h0) ./ std_h0;
    % threshold image at p-value, by setting subthreshold values to 0, this
    % is uncorrected zmap, abs() for two-tailed
    zmap(abs(zmap)<zval) = 0;
    
    % now for cluster and pixel correction
    % initialize matrices for cluster-based correction
    max_cluster_sizes = zeros(1,n_permutes*20);
    % ... and for maximum-pixel based correction
    max_val = zeros(n_permutes*20,2); % "2" for min/max
    zmap_cluster = zmap; % initialize for cluster-correction

    % loop through permutations
    for permi = 1:n_permutes*20

        % take each permutation map, and transform to Z
        threshimg = squeeze(permmaps(areaN,permi,:,:));
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
        temp = sort( reshape(permmaps(areaN,permi,:,:),1,[] ));
        max_val(permi,:) = [ min(temp) max(temp) ];

    end
    % find cluster threshold (need image processing toolbox for this!)
    % based on p-value and null hypothesis distribution
    cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));
    % now find clusters in the real thresholded zmap
    % if they are "too small" set them to zero
    islands = bwconncomp(zmap);
    for i=1:islands.NumObjects
        % if real clusters are too small, remove them by setting to zero!
        if numel(islands.PixelIdxList{i}==i)<cluster_thresh
            zmap_cluster(islands.PixelIdxList{i})=0;
        end
    end
    % find the pixel threshold for lower and upper values, two-tailed
    thresh_lo = prctile(max_val(:,1),100*(pval/2)); % pval/2 percentile of smallest values
    thresh_hi = prctile(max_val(:,2),100-100*(pval/2)); % pval/2 percentile of largest values
    % threshold real data
    pixel_threshmap = diffmap; % initialize for pixel correction
    pixel_threshmap(pixel_threshmap>thresh_lo & pixel_threshmap<thresh_hi) = 0;
    if find(pixel_threshmap~=0) % significance left after pixel correction
        signif.(areas{areaN}).zmap_pixel = pixel_threshmap; 
        if find(zmap_cluster~=0) % now check cluster correction
            signif.(areas{areaN}).zmap_cluster = zmap_cluster;
        end
    elseif find(zmap_cluster~=0) % sig left after cluster correction
        signif.(areas{areaN}).zmap_cluster = zmap_cluster;
    end
end

%% show histograph of maximum cluster sizes

figure(1), clf
hist(max_cluster_sizes,20);
hold on % vert line
yL = get(gca,'YLim'); line([cluster_thresh cluster_thresh],yL,'Color','k','LineStyle',':','LineWidth',2); 
axis off
xlabel('Maximum cluster sizes'), ylabel('Number of observations')
title('Expected cluster sizes under the null hypothesis')
histogram(max_val,500,'FaceColor','b'); % bimodal
hold on % vertical lines
yL = get(gca,'YLim'); line([thresh_lo thresh_hi;thresh_lo thresh_hi],yL,'Color','k','LineStyle',':','LineWidth',2); 

%% plots with multiple comparisons corrections

% load significant areas that survived multiple comparison correction
load('m2_significant_areas.mat')

monkeyN = 2; % only monkey 2 was significant
% significant areas: 8B, dPFC

if monkeyN==1
    areas = {'a8B', 'a9L', 'adPFC', 'avPFC', 'aLIP', 'aMIP', 'aPEC', 'aPG'}; %monkey1
    permmaps = m1_permmaps;
else
    areas = {'a6DR', 'a8AD', 'a8B', 'adPFC', 'aLIP', 'aPE', 'aPEC', 'aPG'}; %monkey2
    permmaps = m2_permmaps;
end

areaN = 3; % area 8B
areaN = 4; % area dPFC

cor = mean( monkey(monkeyN).correct.(areas{areaN}),3 );
inc = mean( monkey(monkeyN).incorrect.(areas{areaN}),3 );

% compute the difference in power between the two conditions
diffmap = cor - inc;
% extract mult comp correction maps if they exist
zmap_cluster = signif.(areas{areaN}).zmap_cluster; % 8B, dPFC 
zmap_pixel = signif.(areas{areaN}).zmap_pixel; %  dPFC

figure(13), clf
subplot(221)
contourf(signalt(times2saveidx),frex,diffmap,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_cluster),1,'linecolor','k');
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
pos = get(cbar,'Position'); lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Power (\muV^2)'; cbar.Label.Position=[pos(1)+1 pos(2)];
cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('Cluster-corrected Sig Regions Outlined Monkey %d, Area %s, p=%f',monkeyN,areas{areaN}(2:end),pval));

subplot(222)
contourf(signalt(times2saveidx),frex,diffmap,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_cluster),1,'linecolor','k');
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
pos = get(cbar,'Position'); lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Power (\muV^2)'; cbar.Label.Position=[pos(1)+1 pos(2)];
cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('Cluster-corrected Sig Regions Outlined Monkey %d, Area %s, p=%f',monkeyN,areas{areaN}(2:end),pval));

subplot(223)
contourf(signalt(times2saveidx),frex,diffmap,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_pixel),1,'linecolor','k');
set(gca,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
% Create arrow
xa = [0.409895833333333 0.427395833333333]; ya = [0.158878504672897 0.114730010384216];
annotation('arrow',xa,ya)
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
pos = get(cbar,'Position'); lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Power (\muV^2)'; cbar.Label.Position=[pos(1)+1 pos(2)];
cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('Pixel-corrected Significant Regions Outlined Monkey %d, Area %s, p=%f',monkeyN,areas{areaN}(2:end),pval));

%% Poster Images
% plot correct and incorrect then diffmap with sig regions outlined after mcc

% load relevant workspace var's
load('m2_pow_significant_areas.mat')

figure(9), clf
% subplot(221)
contourf(signalt(times2saveidx),frex,cor,100,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_pixel),1,'linecolor','w','linewidth',3);
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineWidth',2,'LineStyle',':'); 
set(gca,'clim',[0 25],'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off','FontSize',30)
xlabel('Time (s)','FontSize',34), ylabel('Frequency (Hz)','FontSize',34), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Raw Power (\muV^2)'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-2.5 pos(2)];
% title('Correct, Significant Regions Outlined','FontSize',32);
export_fig test2.png -transparent % no background

% zoom in on clustered area (delay only, 2.5-8.6Hz) area 8B, [-250 130,
% 3.5-11Hz area dPFC), pixel dfpc: [1150 1200ms], [3.5 4Hz] too small
clustert = [1150 1200]; % in ms
clusteridx = dsearchn(times2save',clustert');
freqz = [3.5 4]; % in Hz
freqzidx = dsearchn(frex',freqz');

figure(10), clf
% subplot(221)
contourf(times2save(clusteridx(1):clusteridx(2)),frex(freqzidx(1):freqzidx(2)),diffmap(freqzidx(1):freqzidx(2),clusteridx(1):clusteridx(2)),100,'linecolor','none')
hold on
contour(times2save(clusteridx(1):clusteridx(2)),frex(freqzidx(1):freqzidx(2)),logical(zmap_pixel(freqzidx(1):freqzidx(2),clusteridx(1):clusteridx(2))),1,'linecolor','k','linewidth',3);
% yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineWidth',2,'LineStyle',':'); 
set(gca,'clim',[0 2],'ytick',round(logspace(log10(frex(freqzidx(1))),log10(frex(freqzidx(2))),10)*100)/100,'yscale','log','YMinorTick','off','FontSize',30)
xlabel('Time (ms)','FontSize',34), ylabel('Frequency (Hz)','FontSize',34), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Power Difference (\muV^2)'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
export_fig test2.png -transparent % no background

figure(9), clf
% subplot(222)
contourf(signalt(times2saveidx),frex,inc,100,'linecolor','none');
hold on
contour(signalt(times2saveidx),frex,logical(zmap_pixel),1,'linecolor','w','linewidth',3);
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineWidth',2,'LineStyle',':'); 
set(gca,'clim',[0 25],'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off','FontSize',30)
xlabel('Time (s)','FontSize',34), ylabel('Frequency (Hz)','FontSize',34), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Raw Power (\muV^2)'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-2.5 pos(2)];
% title('Incorrect, Significant Regions Outlined','FontSize',32);
export_fig test2.png -transparent % no background

figure(9), clf
% subplot(223)
contourf(signalt(times2saveidx),frex,diffmap,100,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_pixel),1,'linecolor','k','linewidth',3);
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineWidth',2,'LineStyle',':'); 
set(gca,'clim',[-2 2],'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off','FontSize',30)
xlabel('Time (s)','FontSize',34), ylabel('Frequency (Hz)','FontSize',34), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Power Difference (\muV^2)'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-2.3 pos(2)];
% title(sprintf('M%d Avg FP Coh Pixel-corrected Sig Regions Outlined Cor > Inc p<%1.2f',monkeyN,pval));
export_fig test2.png -transparent % no background

%% now with max-pixel-based thresholding

% find the threshold for lower and upper values, two-tailed
thresh_lo = prctile(max_val(:,1),100*(pval/2)); % pval/2 percentile of smallest values
thresh_hi = prctile(max_val(:,2),100-100*(pval/2)); % pval/2 percentile of largest values

% threshold real data
zmap = diffmap;
zmap(zmap>thresh_lo & zmap<thresh_hi) = 0;

figure(11), clf
subplot(221)
imagesc(times2save,[],diffmap)
set(gca,'clim',[-mean(clim)/5 mean(clim)/5],'ydir','norm')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title(sprintf('TF power map, no thresholding Monkey %d, Area %s, Correct - Incorrect',monkeyN,m2areas{areaN+3}(2:end)));

subplot(222)
imagesc(times2save,[],diffmap)
hold on
contour(times2save,frex,logical(zmap),1,'linecolor','k')
set(gca,'clim',[-mean(clim)/5 mean(clim)/5],'ydir','norm')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title(sprintf('Max-pixel-corrected TF power with contour Monkey %d, Area %s',monkeyN,m2areas{areaN+3}(2:end)));

subplot(223)
imagesc(times2save,[],zmap)
set(gca,'clim',[-mean(clim)/5 mean(clim)/5],'ydir','norm')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title(sprintf('Max-pixel-corrected & thresholded TF z-map Monkey %d, Area %s',monkeyN,m2areas{areaN+3}(2:end)));

