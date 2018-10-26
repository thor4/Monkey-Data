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
    for i=1:numel(monkey(monkeyN).day)
        chan = fieldnames(monkey(monkeyN).day(i).(response));
        % initialize mat to save ea analytic signal for all channels
        % chan x freq x time x trials
        chan_as = zeros(numel(chan),num_frex,length(signalt),size(monkey(monkeyN).day(i).(response).(chan{1})',2));
        for j=1:numel(chan)
            signal = monkey(monkeyN).day(i).(response).(chan{j})'; % change to time-by-trials
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
        for j=1:size(chancombo,1) %total number of combinations
            for k=1:numel(areas)
                if endsWith(chancombo{j,1},chans{k}) %which chan
                    area1=k;
                end
                if endsWith(chancombo{j,2},chans{k}) %which chan
                    area2=k;
                end
            end
            if area1~=area2 % ensure not finding coherence within same area
                % initialize coherence matrix for all frex for area
                spectcoh = zeros(length(frex),length(times2save));
                % pull out analytic signal from each channel
                as1 = squeeze(chan_as(chancomboidx(2,1),:,:,:)); % squeeze to frex x time x trial
                as2 = squeeze(chan_as(chancomboidx(2,2),:,:,:)); % squeeze to frex x time x trial
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
                combo = sprintf('%s_%s',areas{area1},areas{area2});
                combo_rev = sprintf('%s_%s',areas{area2},areas{area1});  
                if isfield(monkey(monkeyN).(response),(combo)) % field exists?
                    % yes, field exists..
                    aTrials = size(monkey(monkeyN).(response).(combo),3);
                    % append session for area in freq x time struct
                    monkey(monkeyN).(response).(combo)(:,:,aTrials+1) = spectcoh;
                elseif isfield(monkey(monkeyN).(response),(combo_rev))
                    aTrials = size(monkey(monkeyN).(response).(combo_rev),3);
                    % append session for area in freq x time struct
                    monkey(monkeyN).(response).(combo_rev)(:,:,aTrials+1) = spectcoh;
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
load('mGoodStableRule1PingRejCoh-AllDays_Cor_Inc_Allchans.mat')

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
fpidx = []; %initialize fp_pairs matrix for each monkey
% pull out only frontal-parietal or parietal-frontal pairs for testing
for comboN=1:numel(combos) %parse all combinations
    for areaN=1:numel(areas) %which area pair
        if startsWith(combos{comboN},areas{areaN})
            area1=areaN;
        elseif endsWith(combos{comboN},areas{areaN})
            area2=areaN;
        end
    end
    if area1<=4 && area2>=5 %frontal-parietal
        fpidx = [fpidx comboN]; %save idx
    elseif area1>=5 && area2<=4 %parietal-frontal
        fpidx = [fpidx comboN];
    end
end

%% plotting the raw difference data

%pull out correct & incorrect coherence avg over trials
cor = mean( monkey(monkeyN).(responses{1}).(combos{fpidx(1)}),3 ); 
inc = mean( monkey(monkeyN).(responses{2}).(combos{fpidx(1)}),3 );
% compute the difference in power between the two conditions
% diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));
diffmap = cor - inc;

% plot examples showing correct, incorrect and their difference in raw
% coherence
figure(15), clf
subplot(221)
imagesc(times2save,[],cor)
set(gca,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), colorbar
title(sprintf('Coherence from Monkey %d, Area %s, Resp Correct',monkeyN,combos{fpidx(1)}(2:end)));

subplot(222)
imagesc(times2save,[],inc)
set(gca,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), colorbar
title(sprintf('Coherence from Monkey %d, Area %s, Resp Incorrect',monkeyN,combos{fpidx(1)}(2:end)));

subplot(223)
imagesc(times2save,[],diffmap)
set(gca,'clim',[-.04,.04],'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
pos = get(cbar,'Position'); lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; cbar.Label.Position=[pos(1)+1 pos(3)-.02];
cbar.TickLabels = ({'-0.04','0.04'});
title(sprintf('Coherence Difference Monkey %d, Area %s, Correct > Incorrect',monkeyN,combos{fpidx(1)}(2:end)));

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

figure(16), clf
subplot(221)
imagesc(times2save,[],cor-cor_base)
set(gca,'clim',clim,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), colorbar
title(sprintf('Baseline-subtracted Coherence from Monkey %d, Area %s, Resp Correct',monkeyN,combos{subplotN}(2:end)));

subplot(222)
imagesc(times2save,[],inc-inc_base)
set(gca,'clim',clim,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), colorbar
title(sprintf('Baseline-subtracted Coherence from Monkey %d, Area %s, Resp Incorrect',monkeyN,combos{subplotN}(2:end)));

subplot(223)
imagesc(times2save,[],diffmap_base)
set(gca,'clim',clim./2,'ydir','n')
set(gca,'ytick',1:4:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),13)*10)/10)
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
pos = get(cbar,'Position'); lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Coherence vs baseline'; cbar.Label.Position=[pos(1)+1 pos(3)-.02];
% cbar.TickLabels = ({'-0.04','0.04'});
title(sprintf('Baseline-subtracted Coherence Difference Monkey %d, Area %s, Correct > Incorrect',monkeyN,combos{subplotN}(2:end)));

%% statistics via permutation testing

% p-value
pval = 0.01;

% convert p-value to Z value
zval = abs(norminv(pval));

% number of permutations
n_permutes = 1000;

% number of frontoparietal (FP) pairs 
n_fppairs = length(fpidx);

% initialize FP pair matrix for monkey 1, [frontal parietal]
m1fppairs = zeros(n_fppairs,2);

% innitialize FP pair diffmap for monkey 1
m1fpdiffmap = zeros(num_frex,length(times2save),n_fppairs);

% get cohmap of all avg coherence over sessions for all frontoparietal
% pairs for each monkey
for fpN=1:n_fppairs
    %pull out correct & incorrect coherence avg over trials
    cor = mean( monkey(monkeyN).(responses{1}).(combos{fpidx(fpN)}),3 ); 
    inc = mean( monkey(monkeyN).(responses{2}).(combos{fpidx(fpN)}),3 );
    % compute the difference in power between the two conditions
    % diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));
    m1fpdiffmap(:,:,fpN) = cor - inc;
    % pull out FP pair
    for areaN=1:numel(areas) %which area pair
        if startsWith(combos{fpidx(fpN)},areas{areaN}) %1st area
            if areaN<=4
                areaF=areaN;
            else, areaP=areaN;
            end
        elseif endsWith(combos{fpidx(fpN)},areas{areaN}) %2nd area
            if areaN<=4
                areaF=areaN;
            else, areaP=areaN;
            end
        end
    end
    m1fppairs(fpN,:) = [areaF areaP]; % save FP pair    
end

% initialize permutation matrix for monkey 1
m1permN = zeros(num_frex,length(times2save),n_permutes);

% now permutation test, on each permutation, create group-avg coherence
% across all FP pairs. this results in maps under the null hypothesis
for permN=1:n_permutes
    % create random vector of 1's and -1's to apply to each FP pair map
    shuffinv = sign(randn(n_fppairs,1));
    shuffinv = reshape(shuffinv,[1 1 n_fppairs]);
    % apply shuffling through multiplication because multiplying by -1 is 
    % algebraically the same as changing from cor - inc to inc - cor, then 
    % calculate mean across all FP pairs
    m1permN(:,:,permN) = mean( m1fpdiffmap.*shuffinv,3 );
    % this represents freq x time points under the null hypothesis
end

% calculate my z-score based on the null hypothesis distribution created
% through group-level permutations
zmap_group = ( mean(m1fpdiffmap,3) - mean(m1permN,3) ) ./ std(m1permN,0,3);

% threshold image at p-value, by setting subthreshold values to 0
zmap_group(abs(zmap_group)<zval) = 0;

%% show non-corrected thresholded maps

% contourf plot template:
% x = 1 x samples
% y = 1 x frequencies
% z = frequencies x samples
% contourf(x,y,z,...)
clim=[-.04,.04];

figure(19), clf
subplot(221)
contourf(signalt(times2saveidx),frex,mean(m1fpdiffmap,3),15,'linecolor','none')
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
% cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('TF Map Monkey %d, Avg Frontoparietal Coherence, Correct > Incorrect',monkeyN));

%areas{m1fppairs(fppair,1)}(2:end),areas{m1fppairs(fppair,2)}(2:end)

subplot(222)
contourf(signalt(times2saveidx),frex,mean(m1fpdiffmap,3),15,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_group),1,'linecolor','k');
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
% cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('Significant Regions Outlined, Avg FP Coherence, Monkey %d, Cor > Inc p=%1.2f',monkeyN,pval));

subplot(223)
contourf(signalt(times2saveidx),frex,zmap_group,15,'linecolor','none')
set(gca,'clim',clim.*100,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Z-score'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-.8 pos(2)];
% cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('Thresholded Uncorrected Z-values Avg FP Coh, Monkey %d, Correct > Incorrect p=%1.2f',monkeyN,pval));

%% corrections for multiple comparisons

if monkeyN==1, permmaps = m1permN;
else, permmaps = m2permN;
end

% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);
% ... and for maximum-pixel based correction
max_val = zeros(n_permutes,2); % "2" for min/max

% compute mean and standard deviation maps under the null hypothesis
mean_h0 = mean(permmaps,3);
std_h0  = std(permmaps,0,3);

% initialize cluster-correction threshold map
zmap_cluster = zmap_group; 

% loop through permutations
for permi = 1:n_permutes

    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(:,:,permi));
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


    % get extreme values (smallest and largest) (not z score)
    temp = sort( reshape(permmaps(:,:,permi),1,[] ));
    max_val(permi,:) = [ min(temp) max(temp) ];

end
% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));
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
thresh_lo = prctile(max_val(:,1),100*(pval/2)); % pval/2 percentile of smallest values
thresh_hi = prctile(max_val(:,2),100-100*(pval/2)); % pval/2 percentile of largest values
% threshold real data
pixel_threshmap = mean(m1fpdiffmap,3);
% pixel-corrected threshold
pixel_threshmap(pixel_threshmap>thresh_lo & pixel_threshmap<thresh_hi) = 0;

%% show histograph of maximum cluster sizes

figure(20), clf
hist(max_cluster_sizes,20);
xlabel('Maximum cluster sizes'), ylabel('Number of observations')
title('Expected cluster sizes under the null hypothesis')
% hist(max_val,50); %bimodal

%% plots with multiple comparisons corrections

figure(20), clf
subplot(221)
contourf(signalt(times2saveidx),frex,mean(m1fpdiffmap,3),15,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(zmap_cluster),1,'linecolor','k');
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
title(sprintf('Monkey %d Avg FP Coh Cluster-corrected Sig Regions Outlined Cor > Inc p=%1.2f',monkeyN,pval));

subplot(222)
contourf(signalt(times2saveidx),frex,zmap_cluster,15,'linecolor','none');
set(gca,'clim',clim.*100,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Z-Score'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-.7 pos(2)];
% cbar.TickLabels = ({'Incorrect','Correct'});
title(sprintf('Monkey %d Cluster-corrected & thresholded FP Coherence z-map, Cor > Inc p=%1.2f',monkeyN,pval));

subplot(223)
contourf(signalt(times2saveidx),frex,mean(m1fpdiffmap,3),15,'linecolor','none')
hold on
contour(signalt(times2saveidx),frex,logical(pixel_threshmap),1,'linecolor','k');
set(gca,'clim',clim,'ytick',round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100,'yscale','log','YMinorTick','off')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), cbar = colorbar; 
lim = get(cbar,'Limits'); cbar.Ticks=lim;
cbar.Label.String = 'Spectral Coherence'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-1.5 pos(2)];
title(sprintf('Monkey %d Avg FP Coh Pixel-corrected Sig Regions Outlined Cor > Inc p=%1.2f',monkeyN,pval));
