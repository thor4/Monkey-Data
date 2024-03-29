function [pp]=arpower(x,Nr,Nl,porder,fs)
% Use Geweke's method to estimate the power spectral density for a data
% matrix
%
%   x is a matrix whose every row is one channel's time series.
%   Nr is the number of realizations
%   Nl is the length of every realization
%   porder is the order of the AR model
%   fs is the sampling frequency
%   [When porder < Nl, the model is repeatedly estimated for each point in the
%   estimation window (having length Nl)]
%
%   pp is the power spectral density (PSD) of this channel 
%   The power spectrum is computed for f=0:fs/2 (the Nyquist range)

% written by SB on 04-25-2017


[L,N] = size(x); %L is the number of channels, and N is the total number of points for the channel       (Nr x Nl)

% fit the ar model   for this channel with estimation window of length Nl
[A2,Z2] = armorf(x,Nr,Nl,porder); 
 % compute the spectral matrix and power for all channels
f_ind=0;
for f = 0:fs/2
    f_ind=f_ind + 1;
    [S2,H2] = spectrum(A2,Z2,porder,f,fs);
 for ichan=1:L
    pp(ichan,f_ind) = abs(S2(ichan,ichan))*2; %1-sided power spectrum; corrected by SB 9/06
 end
end

