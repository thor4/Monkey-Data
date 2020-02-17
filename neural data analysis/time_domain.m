%%%% Time-domain analysis %%%%
% %% Step 1a: Pre-process data using baseline normalization
% 
% %% Step 1b: Pre-process data using a low-pass butterworth filter
% 
% filtOrder = 3; % filter order, higher values = steeper roll off around the 
% %cutoff but too high and you're gonna get phase distortions
% x = 0;
% y = 20;
% sr = 1000;
% 
% %[b,a] = butter(filtOrder,[xx,yy]/(sr/2),'bandpass'); % construct filter; xx 
% %and yy are the lower and upper components of the bandpass, sr is your sampling rate
% [b,a] = butter(filtOrder,20/(sr/2),'low'); %low pass filter for components up to 20Hz
% 
% % filter your data using the coefficients you computed above and save the 
% % output. data should be an #samples x #trials matrix
% filtCFdpfc = filtfilt(b,a,CFdpfc);  
% filtCPpec = filtfilt(b,a,CPpec);  
% filtIFdpfc = filtfilt(b,a,IFdpfc);  
% filtIPpec = filtfilt(b,a,IPpec);
% 
% %baseline
% filtbCFdpfc = filtfilt(b,a,bCFdpfc);  
% filtbCPpec = filtfilt(b,a,bCPpec);  
% filtbIFdpfc = filtfilt(b,a,bIFdpfc);  
% filtbIPpec = filtfilt(b,a,bIPpec);

%% ERP %%
%initialize variables
% homepc:
path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\';
% labpc:
% path = 'C:\\Users\\bryan\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\';
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
dayy = append('d',lp{:});
% [data.(dayy).lfp,data.(dayy).areas] = extractDay(path,monkey,lp{:},1,2,1,1,'match'); %correct match
day = alldays{21}; %assign day
[lfp,areas] = extractDay(path,monkey,day,1,2,1,1,'match'); %correct match
lfp_avg_mv = mean(lfp,3) .* 1e6; % avg over trials & convert to µV (1V = 10^6µV = 1,000,000µV)