% time frequency coherence
%% Step 1: Make Wavelets

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
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
fwhm = logspace(log10(min_fwhm),log10(max_fwhm),length(frex)); % in seconds 

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

%% Create analytic signal using the complex morlet wavelet family for both 
% monkeys, both conditions, all days, all channels and all trials. pull out
% and compute spectral coherence for all channel combinations

% load monkey data
load('mGoodStableRule1PingRej-split_by_Day_BehResp_and_Chan.mat')

% choose which monkey and response you want to pull out
monkeyN = 2; % which monkey (1 or 2)
% resp = 2; % 1 = correct, 0 = incorrect

% begin definining convolution parameters
n_wavelet = length(wavet);
half_of_wavelet_size = floor(n_wavelet/2)+1;

if monkeyN==1
    areas = {'a8B', 'a9L', 'adPFC', 'avPFC', 'aLIP', 'aMIP', 'aPEC', 'aPG'}; %monkey1
    chans = {'8B', '9L', 'dPFC', 'vPFC', 'LIP', 'MIP', 'PEC', 'PG'};
else
    areas = {'a6DR', 'a8AD', 'a8B', 'adPFC', 'aLIP', 'aPE', 'aPEC', 'aPG'}; %monkey2
    chans = {'6DR', '8AD', '8B', 'dPFC', 'LIP', 'PE', 'PEC', 'PG'};
end
% if resp==1, response = 'correct';
% else, response = 'incorrect';
% end

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
% setup response vector
responses=[{'correct'},{'incorrect'}];

tic
for responseN=1:length(responses) % loop through correct + incorrect
    response = responses{responseN}; % which response
    % pull out spectral coherence across all days and append per session,
    % per condition, averaged across trials
    for dayN=1:numel(monkey(monkeyN).day)
        chan = fieldnames(monkey(monkeyN).day(dayN).(response));
        % initialize mat to save ea analytic signal for all channels
        % chan x freq x time x trials
        chan_as = zeros(numel(chan),num_frex,length(signalt),size(monkey(monkeyN).day(dayN).(response).(chan{1})',2));
        for j=1:numel(chan)
            signal = monkey(monkeyN).day(dayN).(response).(chan{j})'; % change to time-by-trials
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
            % initialize mat to save ea analytic signal for all frex
            % freq x time x trials
            as_ = zeros(num_frex,length(signalt),size(signal,2));
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
                % save as for ea freq x time x trials 
                as_(fi,:,:) = as;
            end
            chan_as(j,:,:,:) = as_; % save as for ea chan x freq x time x trial
        end
        chancombo = combnk(chan,2); %all chan combinations
        chancomboidx = combnk(1:numel(chan),2); %all chan index combinations
        % count number of incorrect trials (any chan doesn't matter)
        ntrials = size( monkey(monkeyN).day(dayN).(responses{2}).(chan{j}),1 );
        npermutes=20; %how many times to randomly sample correct trials?
        for comboN=1:size(chancombo,1) %total number of combinations
            areaF=[]; areaP=[]; %initialize areas
            for areaN=1:numel(areas) %which area pair
                if endsWith(chancombo{comboN,1},chans{areaN}) %which chan
                    if areaN<=4
                        areaF=areaN;
                    else, areaP=areaN;
                    end
                end
                if endsWith(chancombo{comboN,2},chans{areaN}) %which chan
                    if areaN<=4
                        areaF=areaN;
                    else, areaP=areaN;
                    end
                end
            end
            if ~isempty(areaF) && ~isempty(areaP) % ensure not finding coherence within same area
                % initialize coherence matrix for all frex for area
                spectcoh = zeros(length(frex),length(times2save));
                % pull out analytic signal from each channel
                as1 = squeeze(chan_as(chancomboidx(comboN,1),:,:,:)); % squeeze to frex x time x trial
                as2 = squeeze(chan_as(chancomboidx(comboN,2),:,:,:)); % squeeze to frex x time x trial
                if responseN==1 % correct trials?
                    spectcohpermmat = zeros(npermutes,length(frex),length(times2save));
                    for permN=1:npermutes
                        % initialize coherence matrix for all frex for area
                        spectcohperm = zeros(length(frex),length(times2save));
                        % randomly sample from correct trials (decimation) 
                        % to match incorrect
                        randcoridx = randperm( size( as1,3 ),ntrials );
                        % create permutation for each analytic signal using
                        % same random indices, still frex x time x trial
                        as1perm = as1(:,:,randcoridx);
                        as2perm = as2(:,:,randcoridx);
                        for fi=1:length(frex) %compute spectral coherence for ea freq
                            sig1 = squeeze(as1perm(fi,:,:));
                            sig2 = squeeze(as2perm(fi,:,:));
                            % sig1 & sig2 are now a time x trial complex matrix
                            % compute avg power and avg cross-spectral power over trials
                            spec1 = mean(sig1.*conj(sig1),2); % faster to compute power than abs.^2
                            spec2 = mean(sig2.*conj(sig2),2);
                            specX = abs(mean(sig1.*conj(sig2),2)).^2;
                            % yields three time x 1 matrices, setting up for coherence over trials
                            % to show connectivity evolution over time
                            % compute spectral coherence, using only requested time points
                            spectcohperm(fi,:) = specX(times2saveidx)./(spec1(times2saveidx).*spec2(times2saveidx));
                        end
                        % save as permutation x frex x down-sampled time
                        spectcohpermmat(permN,:,:) = spectcohperm;
                    end
                    %take the mean over permutations to yield avg decimated
                    %correct frex x trial matrix
                    spectcoh = squeeze(mean(spectcohpermmat,1)); 
                else
                    for fi=1:length(frex) %compute spectral coherence for ea freq
                        sig1 = squeeze(as1(fi,:,:));
                        sig2 = squeeze(as2(fi,:,:));
                        % sig1 & sig2 are now a time x trial complex matrix
                        % compute avg power and avg cross-spectral power over trials
                        spec1 = mean(sig1.*conj(sig1),2); % faster to compute power than abs.^2
                        spec2 = mean(sig2.*conj(sig2),2);
                        specX = abs(mean(sig1.*conj(sig2),2)).^2;
                        % yields three time x 1 matrices, setting up for coherence over trials
                        % to show connectivity evolution over time
                        % compute spectral coherence, using only requested time points
                        spectcoh(fi,:) = specX(times2saveidx)./(spec1(times2saveidx).*spec2(times2saveidx));
                    end
                end
                combo = sprintf('%s_%s',areas{areaF},areas{areaP});
                if isfield(monkey(monkeyN).(response),(combo)) % field exists?
                    % yes, field exists..
                    aTrials = size(monkey(monkeyN).(response).(combo),3);
                    % append session for area in freq x time struct
                    monkey(monkeyN).(response).(combo)(:,:,aTrials+1) = spectcoh;
                else
                    monkey(monkeyN).(response).(combo) = spectcoh; % freq x time x trial
                end
            end
        end
    end
end
toc

%% setting up for plots & permutation testing

% load coherence data for both monkeys
% load('mGoodStableRule1PingRejCoh-AllDays_Cor_Inc_Allchans.mat') % legacy
load('mcohfp_corinc.mat') %has only FP combinations saved

monkeyN = 2; %which monkey? 1 or 2
responses = [{'correct'},{'incorrect'}]; 
responseN = 1; %which response? 1 or 2

times2save = -400:10:1200; % in ms
% time vector converted to indices
times2saveidx = dsearchn((signalt.*1000)',times2save');
% define baseline time
baset = [-400 -100]; % in ms
baseidx = dsearchn(times2save',baset'); % search down-sampled time

if monkeyN==1
    areas = {'a8B', 'a9L', 'adPFC', 'avPFC', 'aLIP', 'aMIP', 'aPEC', 'aPG'}; %monkey1
    combos = fieldnames(monkey(monkeyN).(responses{responseN}));
else
    areas = {'a6DR', 'a8AD', 'a8B', 'adPFC', 'aLIP', 'aPE', 'aPEC', 'aPG'}; %monkey2
    combos = fieldnames(monkey(monkeyN).(responses{responseN}));
end
% legacy for when all combos are saved, not sure FP
% fpidx = []; %initialize fp_pairs matrix for each monkey
% % pull out only frontal-parietal or parietal-frontal pairs for testing
% for comboN=1:numel(combos) %parse all combinations
%     for areaN=1:numel(areas) %which area pair
%         if startsWith(combos{comboN},areas{areaN})
%             area1=areaN;
%         elseif endsWith(combos{comboN},areas{areaN})
%             area2=areaN;
%         end
%     end
%     if area1<=4 && area2>=5 %frontal-parietal
%         fpidx = [fpidx comboN]; %save idx
%     elseif area1>=5 && area2<=4 %parietal-frontal
%         fpidx = [fpidx comboN];
%     end
% end

%% plotting the raw difference data from one FP combo

%pull out correct & incorrect coherence avg over trials
cor = mean( monkey(monkeyN).(responses{1}).(combos{1}),3 ); 
inc = mean( monkey(monkeyN).(responses{2}).(combos{1}),3 );
% compute the difference in power between the two conditions
% diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));
diffmap = cor - inc;

% plot examples showing correct, incorrect and their difference in raw
% coherence
figure(24), clf
subplot(221)
imagesc(times2save,[],cor)
set(gca,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), colorbar
title(sprintf('Coherence from Monkey %d, Area %s, Resp Correct',monkeyN,combos{1}(2:end)));

subplot(222)
imagesc(times2save,[],inc)
set(gca,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), colorbar
title(sprintf('Coherence from Monkey %d, Area %s, Resp Incorrect',monkeyN,combos{1}(2:end)));

subplot(223)
imagesc(times2save,[],diffmap)
set(gca,'clim',[-.04,.04],'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
pos = get(cbar,'Position'); lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; cbar.Label.Position=[pos(1)+1 pos(3)-.02];
cbar.TickLabels = ({'-0.04','0.04'});
title(sprintf('Coherence Difference Monkey %d, Area %s, Correct > Incorrect',monkeyN,combos{1}(2:end)));

% plot examples showing correct, incorrect and their difference for
% baseline-corrected
% compute baseline for correct & incorrect responses averaged over trials 
% then timepoints, leaving freq x 1 matrices of baseline
cor_base = mean( cor(:,baseidx(1):baseidx(2)),2 );
inc_base = mean( inc(:,baseidx(1):baseidx(2)),2 );
% create appropriately-sized matrices to subtract from avg coherence, where
% cor is freq x time matrix
cor_base = repmat(cor_base,1,size(cor,2));
inc_base = repmat(inc_base,1,size(inc,2));
% calculate diffmap of baseline-subtracted conditions
diffmap_base = (cor-cor_base) - (inc-inc_base);
clim = [-.1 .1]; %set color limits

figure(1), clf
subplot(221)
imagesc(times2save,[],cor-cor_base)
set(gca,'clim',clim,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), colorbar
title(sprintf('Baseline-subtracted Coherence from Monkey %d, Area %s, Resp Correct',monkeyN,combos{1}(2:end)));

subplot(222)
imagesc(times2save,[],inc-inc_base)
set(gca,'clim',clim,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), colorbar
title(sprintf('Baseline-subtracted Coherence from Monkey %d, Area %s, Resp Incorrect',monkeyN,combos{1}(2:end)));

subplot(223)
imagesc(times2save,[],diffmap_base)
set(gca,'clim',clim./2,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
pos = get(cbar,'Position'); lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Coherence vs baseline'; cbar.Label.Position=[pos(1)+1 pos(3)-.02];
% cbar.TickLabels = ({'-0.04','0.04'});
title(sprintf('Baseline-subtracted Coherence Difference Monkey %d, Area %s, Correct > Incorrect',monkeyN,combos{1}(2:end)));

%% pull out avg coherence for all frontoparietal pairs under correct, 
% incorrect & difference conditions + baseline-corrected of same

% number of frontoparietal (FP) pairs 
% n_fppairs = length(fpidx); %legacy
n_fppairs = length(combos);

% initialize FP pair matrix for monkey 1, [frontal parietal]
m2fppairs = zeros(n_fppairs,2);

% initialize FP pair for correct, incorrect & differences for monkey
m2fpcor = zeros(num_frex,length(times2save),n_fppairs);
m2fpinc = zeros(num_frex,length(times2save),n_fppairs);
m2fpdiffmap = zeros(num_frex,length(times2save),n_fppairs);

% initialize baseline-corrected FP pair for correct, incorrect & 
% differences for monkey
m2fpcor_base = zeros(num_frex,length(times2save),n_fppairs);
m2fpinc_base = zeros(num_frex,length(times2save),n_fppairs);
m2fpdiffmap_base = zeros(num_frex,length(times2save),n_fppairs);

for fpN=1:n_fppairs
    %pull out correct & incorrect coherence avg over sessions & channels,
    %previously avg'd over trials, leaves freq x time
    m2fpcor(:,:,fpN) = mean( monkey(monkeyN).(responses{1}).(combos{fpN}),3 ); 
    m2fpinc(:,:,fpN) = mean( monkey(monkeyN).(responses{2}).(combos{fpN}),3 );
    % compute the difference in power between the two conditions
    % diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));
    m2fpdiffmap(:,:,fpN) = m2fpcor(:,:,fpN) - m2fpinc(:,:,fpN);
    % compute baseline for correct & incorrect responses averaged over trials 
    % then timepoints, leaving freq x 1 matrices of baseline
    cor_base = mean(m2fpcor(:,baseidx(1):baseidx(2),fpN),2);
    inc_base = mean(m2fpinc(:,baseidx(1):baseidx(2),fpN),2);
    % create appropriately-sized matrices to subtract from avg coherence, 
    % where cor is freq x time matrix
    cor_base = repmat(cor_base,1,size(m2fpcor,2));
    inc_base = repmat(inc_base,1,size(m2fpinc,2));
    % save baseline-corrected avg coh for all fp pairs across all frex for
    % both conditions: correct & incorrect
    m2fpcor_base(:,:,fpN) = m2fpcor(:,:,fpN)-cor_base;
    m2fpinc_base(:,:,fpN) = m2fpinc(:,:,fpN)-inc_base;
    % calculate diffmap of baseline-subtracted conditions
    m2fpdiffmap_base(:,:,fpN) = m2fpcor_base(:,:,fpN) - m2fpinc_base(:,:,fpN);
    % pull out FP pair
    for areaN=1:numel(areas) %which area pair
        if startsWith(combos{fpN},areas{areaN}) %1st area
            if areaN<=4
                areaF=areaN;
            else, areaP=areaN;
            end
        elseif endsWith(combos{fpN},areas{areaN}) %2nd area
            if areaN<=4
                areaF=areaN;
            else, areaP=areaN;
            end
        end
    end
    m2fppairs(fpN,:) = [areaF areaP]; % save FP pair    
end

%% plot raw correct, incorrect & difference coherence spectrograms

clim = [0 0.2];

figure(5), clf
subplot(221)
contourf(signalt(times2saveidx),frex,mean(m2fpcor,3),15,'linecolor','none')
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
title(sprintf('Monkey %d TF Map, Avg Frontoparietal Coherence, Correct',monkeyN));

subplot(222)
contourf(signalt(times2saveidx),frex,mean(m2fpinc,3),15,'linecolor','none')
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1 pos(2)];
title(sprintf('Monkey %d TF Map, Avg Frontoparietal Coherence, Incorrect',monkeyN));

subplot(223)
contourf(signalt(times2saveidx),frex,mean(m2fpdiffmap,3),15,'linecolor','none')
set(gca,'clim',[-0.04 0.04],'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
title(sprintf('Monkey %d TF Map, Avg Frontoparietal Coherence, Correct - Incorrect',monkeyN));

% plot examples showing correct, incorrect and their difference for
% baseline-corrected

clim = [-0.1 0.1];

figure(6), clf
subplot(221)
contourf(signalt(times2saveidx),frex,mean(m2fpcor_base,3),15,'linecolor','none')
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Baseline-subtracted Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
title(sprintf('Monkey %d TF Map, Avg Frontoparietal Baseline-subtracted Coherence, Correct',monkeyN));

subplot(222)
contourf(signalt(times2saveidx),frex,mean(m2fpinc_base,3),15,'linecolor','none')
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Baseline-subtracted Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1 pos(2)];
title(sprintf('Monkey %d TF Map, Avg Frontoparietal Baseline-subtracted Coherence, Incorrect',monkeyN));

subplot(223)
contourf(signalt(times2saveidx),frex,mean(m2fpdiffmap_base,3),15,'linecolor','none')
hold on %vertical line
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineStyle',':'); 
set(gca,'clim',[-0.04 0.04],'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Baseline-subtracted Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
title(sprintf('Monkey %d TF Map, Avg Frontoparietal Baseline-subtracted Coherence, Correct > Incorrect',monkeyN));

%% statistics via permutation testing

% p-value
pval = 0.05;

% convert p-value to Z value
zval = abs(norminv(pval));

% meta-permutation test 
nmeta = 20; % 20 permutation tests
% number of permutations
n_permutes = 1000;

% initialize meta-permutation matrix
zmap_groupmeta = zeros(nmeta,num_frex,length(times2save));
% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);
% initialize cluster meta permutation matrix to store all vals 
max_cluster_meta = [];
% initialize max val meta permutation matrix
max_val_meta = [];


for metaN=1:nmeta
    % initialize permutation matrix for monkey 1
    m2permN = zeros(num_frex,length(times2save),n_permutes);
    % initialize max val permutation matrix
    max_val_perm = zeros(n_permutes,2);
    % now permutation test, on each permutation, create group-avg coherence
    % across all FP pairs. this results in maps under the null hypothesis
    for permN=1:n_permutes
        % create random vector of 1's and -1's to apply to each FP pair map
        shuffinv = sign(randn(n_fppairs,1));
        shuffinv = reshape(shuffinv,[1 1 n_fppairs]);
        % apply shuffling through multiplication because multiplying by -1 is 
        % algebraically the same as changing from cor - inc to inc - cor, then 
        % calculate mean across all FP pairs
        m2permN(:,:,permN) = mean( m2fpdiffmap.*shuffinv,3 );
        % this represents freq x time points under the null hypothesis
        % get extreme values (smallest and largest) (not z score)
        temp = sort( reshape(m2permN(:,:,permN),1,[] ));
        max_val_perm(permN,:) = [ min(temp) max(temp) ];
    end
    % compute mean and standard deviation maps under the null hypothesis
    mean_h0 = mean(m2permN,3);
    std_h0  = std(m2permN,0,3);
    % calculate my z-score based on the null hypothesis distribution created
    % through group-level permutations
    zmap_groupmeta(metaN,:,:) = ( mean(m2fpdiffmap,3) - mean_h0 ) ./ std_h0;
    % take each permutation map, and transform to Z
    threshimg = (m2permN-mean_h0)./std_h0;
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    for threshN=1:n_permutes
        threshimgperm = threshimg(:,:,threshN);
        % find clusters (need image processing toolbox for this!)
        islands = bwconncomp(threshimgperm);
        if numel(islands.PixelIdxList)>0
            % count sizes of clusters
            tempclustsizes = cellfun(@length,islands.PixelIdxList);
            % store size of biggest cluster
            max_cluster_sizes(threshN) = max(tempclustsizes);
        end
    end
    % save max values from each permutation test for pixel & cluster corr
    max_val_meta = [max_val_meta; max_val_perm]; % pixel correction
    max_cluster_meta = [max_cluster_meta max_cluster_sizes]; % cluster corr
end

% complete meta permutation test by avg over all permutation tests, leaves 
% frex x time of my observed z-scores
zmap_group = squeeze(mean(zmap_groupmeta,1));

% threshold image at p-value, by setting subthreshold values to 0
zmap_group(abs(zmap_group)<zval) = 0;

%% show non-corrected thresholded maps

% contourf plot template:
% x = 1 x samples
% y = 1 x frequencies
% z = frequencies x samples
% contourf(x,y,z,...)
clim=[-.04,.04];

figure(3), clf
subplot(221)
contourf(signalt(times2saveidx),frex,mean(m2fpdiffmap,3),15,'linecolor','none')
hold on %vertical line
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineStyle',':'); 
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
title(sprintf('M%d TF Map, Avg Frontoparietal Coherence, Correct - Incorrect',monkeyN));

subplot(222)
contourf(signalt(times2saveidx),frex,mean(m2fpdiffmap,3),15,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_group),1,'linecolor','k');
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineStyle',':'); 
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
title(sprintf('M%d Sig Regions Outlined, Avg FP Coherence, Cor-Inc p<%1.2f',monkeyN,pval));

subplot(223)
contourf(signalt(times2saveidx),frex,zmap_group,15,'linecolor','none')
hold on % vertical line
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineStyle',':'); 
set(gca,'clim',clim.*100,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Z-score'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-.8 pos(2)];
title(sprintf('M%d Thresholded Uncorrected Z-values Avg FP Coh, Cor-Inc p<%1.2f',monkeyN,pval));

%% corrections for multiple comparisons

% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_meta,100-(100*pval));
% initialize cluster zmap to real thresholded map
zmap_cluster = zmap_group;
% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap_cluster);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<cluster_thresh
        % cluster-corrected threshold
        zmap_cluster(islands.PixelIdxList{i})=0;
    end
end
% find the pixel threshold for lower and upper values, two-tailed
thresh_lo = prctile(max_val_meta(:,1),100*(pval/2)); % pval/2 percentile of smallest values
thresh_hi = prctile(max_val_meta(:,2),100-100*(pval/2)); % pval/2 percentile of largest values
% threshold real data
pixel_threshmap = mean(m2fpdiffmap,3);
% pixel-corrected threshold
pixel_threshmap(pixel_threshmap>thresh_lo & pixel_threshmap<thresh_hi) = 0;

%% show histograph of maximum cluster sizes

figure(1), clf
histogram(max_cluster_meta,200,'FaceColor','b');
hold on % vert line
yL = get(gca,'YLim'); line([cluster_thresh cluster_thresh],[yL(1) yL(2)/4],'Color','k','LineStyle',':','LineWidth',2); 
xlabel('Maximum cluster sizes'), ylabel('Number of observations')
title('Expected cluster sizes under the null hypothesis')
hist(max_val_meta,'FaceColor','b',1000); %bimodal
histogram(max_val_meta,500,'FaceColor','b');
hold on % vertical lines
yL = get(gca,'YLim'); line([thresh_lo thresh_hi;thresh_lo thresh_hi],yL,'Color','k','LineStyle',':','LineWidth',2); 
axis off
export_fig test2.png -transparent % no background

%% plots with multiple comparisons corrections


figure(8), clf
subplot(221)
contourf(signalt(times2saveidx),frex,mean(m2fpdiffmap_base,3),15,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_cluster),1,'linecolor','k');
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineStyle',':'); 
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
title(sprintf('M%d Avg Base-sub FP Coh Cluster-corrected Sig Regions Outlined Cor-Inc p<%1.2f',monkeyN,pval));

subplot(222)
contourf(signalt(times2saveidx),frex,zmap_cluster,15,'linecolor','none');
hold on %vertical line
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineStyle',':'); 
set(gca,'clim',clim.*100,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Z-Score'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-.7 pos(2)];
title(sprintf('M%d Cluster-corrected & thresholded Base-sub FP Coherence z-map, Cor-Inc p<%1.2f',monkeyN,pval));

subplot(223)
contourf(signalt(times2saveidx),frex,mean(m2fpdiffmap_base,3),15,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(pixel_threshmap),1,'linecolor','k');
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineStyle',':'); 
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
title(sprintf('M%d Avg Base-sub FP Coh Pixel-corrected Sig Regions Outlined Cor-Inc p<%1.2f',monkeyN,pval));

%% Poster Images
% plot correct and incorrect then diffmap with sig regions outlined after mcc

% zmap_cluster for cluster correction var
% pixel_threshmap for pixel correction var

load('sig_coherence.mat') % contains all coherence var's from analysis, correction var's are for baseline-sub
% to go back, be sure to remove the _base from Lines 510, 521 & 614

figure(11), clf
% subplot(221)
contourf(signalt(times2saveidx),frex,mean(m2fpcor,3),100,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_cluster),1,'linecolor','w','linewidth',3);
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineWidth',2,'LineStyle',':'); 
set(gca,'clim',[0 .13],'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off','FontSize',30)
xlabel('Time (s)','FontSize',34), ylabel('Frequency (Hz)','FontSize',34), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-3.5 pos(2)];
% title('Correct, Significant Regions Outlined','FontSize',32);
export_fig test2.png -transparent % no background

% zoom in on clustered area ([50 1200], 3.5-11Hz baseline-subtracted coh), 
% ([-400 1200], 3.5-18Hz raw coh)
clustert = [-400 1200]; % in ms
clusteridx = dsearchn(times2save',clustert');
freqz = [3.5 18]; % in Hz
freqzidx = dsearchn(frex',freqz');

figure(12), clf
% subplot(221)
contourf(times2save(clusteridx(1):clusteridx(2)),frex(freqzidx(1):freqzidx(2)),mean(m2fpinc(freqzidx(1):freqzidx(2),clusteridx(1):clusteridx(2),:),3),100,'linecolor','none')
hold on
contour(times2save(clusteridx(1):clusteridx(2)),frex(freqzidx(1):freqzidx(2)),logical(zmap_cluster(freqzidx(1):freqzidx(2),clusteridx(1):clusteridx(2))),1,'linecolor','w','linewidth',3);
yL = get(gca,'YLim'); line([0 500;0 500],yL,'Color','k','LineWidth',2,'LineStyle',':'); 
set(gca,'clim',[0 0.13],'ytick',round(logspace(log10(frex(freqzidx(1))),log10(frex(freqzidx(2))),10)*100)/100,'yscale','log','YMinorTick','off','FontSize',30)
xlabel('Time (ms)','FontSize',34), ylabel('Frequency (Hz)','FontSize',34), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-3.5 pos(2)];
export_fig test2.png -transparent % no background

figure(11), clf
% subplot(221)
contourf(signalt(times2saveidx),frex,mean(m2fpinc,3),100,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_cluster),1,'linecolor','w','linewidth',3);
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineWidth',2,'LineStyle',':'); 
set(gca,'clim',[0 0.13],'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off','FontSize',30)
xlabel('Time (s)','FontSize',34), ylabel('Frequency (Hz)','FontSize',34), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-3.5 pos(2)];
% title('Correct, Significant Regions Outlined','FontSize',32);
export_fig test2.png -transparent % no background

figure(11), clf
% subplot(221)
contourf(signalt(times2saveidx),frex,mean(m2fpdiffmap_base,3),100,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_cluster),1,'linecolor','k','linewidth',3);
yL = get(gca,'YLim'); line([0 .5;0 .5],yL,'Color','k','LineWidth',2,'LineStyle',':'); 
set(gca,'clim',[-0.02 0.02],'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off','FontSize',30)
xlabel('Time (s)','FontSize',34), ylabel('Frequency (Hz)','FontSize',34), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Coherence Difference'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-4 pos(2)];
% title('Correct, Significant Regions Outlined','FontSize',32);
export_fig test2.png -transparent % no background