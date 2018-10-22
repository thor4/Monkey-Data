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
% monkeys, both conditions, all days, all channels and all trials

% load monkey data
load('mGoodStableRule1PingRej-split_by_Day_BehResp_and_Chan.mat')

% choose which monkey and response you want to pull out
monkeyN = 1; % which monkey (1 or 2)
resp = 1; % 1 = correct, 0 = incorrect

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
if resp==1, response = 'correct';
else, response = 'incorrect';
end

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

% pull out spectral coherence across all days and append per session,
% per condition, averaged across trials
for i=1:numel(monkey(monkeyN).day)
    chan = fieldnames(monkey(monkeyN).day(i).correct);
    chancombo = combnk(chan,2); %all chan combinations
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
            % pull out signal from each channel
            sig1 = monkey(monkeyN).day(i).(response).(chancombo{j,1})'; % change to time-by-trials
            sig2 = monkey(monkeyN).day(i).(response).(chancombo{j,2})';
            reflectsig1_all = zeros(size(sig1,1)+2*n_wavelet,size(sig1,2)); %initialize reflected signals mat
            reflectsig2_all = zeros(size(sig2,1)+2*n_wavelet,size(sig2,2)); 
            % reflect all signals
            for signalN=1:size(sig1,2)
                reflectsig1 = [ sig1(n_wavelet:-1:1,signalN); sig1(:,signalN); sig1(end:-1:end-n_wavelet+1,signalN); ];        
                reflectsig1_all(:,signalN) = reflectsig1;
                reflectsig2 = [ sig2(n_wavelet:-1:1,signalN); sig2(:,signalN); sig2(end:-1:end-n_wavelet+1,signalN); ];        
                reflectsig2_all(:,signalN) = reflectsig2;
            end
            % concatenate into a super-trial
            reflectsig1_supertri = reshape(reflectsig1_all,1,[]); % reshape to 1D time-trials
            reflectsig2_supertri = reshape(reflectsig2_all,1,[]); 
            % step 1: finish defining convolution parameters
            n_data = length(reflectsig1_supertri); % time*trials
            n_convolution = n_wavelet+n_data-1;
            % step 2: take FFTs
            fft_data1 = fft(reflectsig1_supertri,n_convolution); % all trials for chan
            fft_data2 = fft(reflectsig2_supertri,n_convolution); % all trials for chan
            for fi=1:length(frex)
                % FFT of wavelet
                fft_wavelet = fft(wavelets(fi,:),n_convolution);
                % step 3: normalize kernel by scaling amplitudes to one in the 
                % frequency domain. prevents amplitude from decreasing with 
                % increasing frequency. diff from 1/f scaling
                fft_wavelet = fft_wavelet ./ max(fft_wavelet);
                % step 4: point-wise multiply and take iFFT
                as1 = ifft( fft_data1.*fft_wavelet ); % analytic signal
                as2 = ifft( fft_data2.*fft_wavelet ); % analytic signal
                % step 5: trim wings
                as1 = as1(half_of_wavelet_size:end-half_of_wavelet_size+1);
                as2 = as2(half_of_wavelet_size:end-half_of_wavelet_size+1);
                % step 6: reshape back to reflected time-by-trials
                as1 = reshape(as1,size(reflectsig1_all,1),size(reflectsig1_all,2));
                as2 = reshape(as2,size(reflectsig2_all,1),size(reflectsig2_all,2));
                % step 7: chop off the reflections
                as1 = as1(n_wavelet+1:end-n_wavelet,:);
                as2 = as2(n_wavelet+1:end-n_wavelet,:);
                % as is now a time x trial complex matrix
                % compute avg power and avg cross-spectral power over trials
                spec1 = mean(as1.*conj(as1),2); % faster to compute power than abs.^2
                spec2 = mean(as2.*conj(as2),2);
                specX = abs(mean(as1.*conj(as2),2)).^2;
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
