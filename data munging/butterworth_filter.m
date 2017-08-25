filtOrder = 3; % filter order, higher values = steeper roll off around the 
%cutoff but too high and you're gonna get phase distortions
x = 0;
y = 20;
sr = 1000;

[b,a] = butter(filtOrder,[0,20]/(sr/2),'bandpass'); % construct filter; xx 
%and yy are the lower and upper components of the bandpass, sr is your sampling rate
[b,a] = butter(filtOrder,20/(sr/2),'low');

% filter your data using the coefficients you computed above and save the 
% output. data should be an #samples x #trials matrix
filtI = filtfilt(b,a,incorrect);  
filtC = filtfilt(b,a,correct);
filtIb = filtfilt(b,a,incorrectBase);  
filtCb = filtfilt(b,a,correctBase);