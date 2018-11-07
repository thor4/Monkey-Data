% run tf.mat up through Line 149 when signal is defined
signal=signal(:,1);
n_data = length(signal);
n_convolution = n_wavelet+n_data-1;
% step 2: take FFTs
fft_data = fft(signal,n_convolution); % all trials for chan
fi=8; %which frequency frex(8) ~= 6.2Hz
fft_wavelet = fft(wavelets(fi,:),n_convolution)'; % transpose since only 1 trial
% step 3: normalize kernel by scaling amplitudes to one in the 
% frequency domain. prevents amplitude from decreasing with 
% increasing frequency. diff from 1/f scaling
fft_wavelet = fft_wavelet ./ max(fft_wavelet);
% step 4: point-wise multiply and take iFFT
as = ifft( fft_data.*fft_wavelet ); % analytic signal
% step 5: trim wings
as = as(half_of_wavelet_size:end-half_of_wavelet_size+1);



figure(1), clf

plot(signalt,signal,'k','linew',2)
axis off
export_fig m1_day1_8b_cor_t1.png -transparent % no background

plot(wavet,real(wavelets(fi,:)),'m','linew',5)
axis off
export_fig wavelet6.2Hz_win.png -transparent % no background

plot(signalt,real(as),'r','linew',3)
xlabel('Time (s)'), ylabel('Amplitude')
title([ 'Real part of wavelet at ' num2str(frex) ' Hz' ])
axis off
export_fig m1_day1_8b_cor_t1_6.2hz-filter_real.png -transparent % no background

'm',time,imag(cmw)

subplot(212)
plot(time,imag(cmw),'r--','linew',3)
xlabel('Time (s)'), ylabel('Amplitude')
title([ 'Imaginary part of wavelet at ' num2str(frex) ' Hz' ])
export_fig wavelet_6hz-imag.png -transparent % no background