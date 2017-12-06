figure 
plot(deltaIdx, mean(bdPFCpowCorR1deltasampled,1), deltaIdx, mean(bdPFCpowIncR1delta,1), ':', 'LineWidth', 2);
figure 
plot(bdPFCfreqCorR1(highgamma3Idx), mean(bdPFCpowCorR1highgamma3,1), bdPFCfreqIncR1(highgamma3Idx), mean(bdPFCpowIncR1highgamma3,1), ':', 'LineWidth', 2);
% shadedErrorBar(x,y,{@mean,@std},'lineprops','-r','transparent',1);
figure
shadedErrorBar(deltaIdx,mean(10.*log10(bdPFCpowCorR1deltasampled),1),std(10.*log10(bdPFCpowCorR1deltasampled)),'lineprops','-b');
hold on
shadedErrorBar(deltaIdx,mean(10.*log10(bdPFCpowIncR1delta),1),std(10.*log10(bdPFCpowIncR1delta)),'lineprops','-r');
hold off
title('LFP FFT/DFT Spectrum During Delay');
xlabel('Frequency (Hz)');
ylabel('Baseline-normalized Power (dB)');
legend('Correct Trials','Incorrect Trials');

figure
shadedErrorBar(highgamma3Idx,mean(10.*log10(bdPFCpowCorR1highgamma3),1),std(10.*log10(bdPFCpowCorR1highgamma3)),'lineprops','-b');
hold on
shadedErrorBar(highgamma3Idx,mean(10.*log10(bdPFCpowIncR1highgamma3),1),std(10.*log10(bdPFCpowIncR1highgamma3)),'lineprops','-r');
hold off

clf


rng(42);
%randomly sample with replacement and take the mean across trials for all
%bands. 
bdPFCidxCorR1 = randi(size(bdPFCpowCorR1delta,1),size(bdPFCpowIncR1delta,1),1);
bdPFCpowCorR1deltasampled = bdPFCpowCorR1delta(bdPFCidxCorR1,:);