% New analysis on the Alpha Inhibition Hypothesis 4/8/21
% 
% Question: does the power during different stages of the task depend on the location of the stimulus?
% 
% use analytic signal created from decompose_all_trial_lfps_to_as.m
%
% identify for every trial the rule and where the stimulus was presented, 
% which corner of the triangle. the triangle was flipped per day
% [0,1,2] is one orientation of the triangle and [3,4,5] is the other
%
% output: a single figure per channel with 3 subplots showing the avg power
% for each stimulus location for a specific rule across all correct trials



% init vars, run this twice (once for each monkey):
%   homepc:
path = 'G:\\monkey_data\\';
%   labpc:
%   path = 'C:\\Users\\bryan\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\'
monk = 1; %1 = clark, 2 = betty
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
        %STOPPED HERE
        
        
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
end



