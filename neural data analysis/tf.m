%Time-frequency decomposition via chronux
%load monkey data
load('mGoodStableRule1PingRej-split_by_Day_BehResp_and_Chan.mat') 
movingwin=[0.5 0.05]; % set the moving window dimensions, 500ms window 50ms step
params.Fs=1000; % sampling frequency
params.fpass=[2 60]; % frequency of interest
params.tapers=[3 5]; % tapers
params.trialave=0; % no average over trials
params.err=0; % no error computation

%initialize all variables
i=1; %day 1
j=1; %resp 1
k=1; %chan 1
l=1; %trial 1
resp = fieldnames(monkey(1).day(1));
times2save = -500+(movingwin(1)*1000)/2-1:movingwin(2)*1000:1310-(movingwin(1)*1000)/2; 
%allows room for half of time_window on front & back

d = 'd%i'; %day

chan = fieldnames(monkey(1).day(i).(resp{j}));
C = combnk(chan,2); %all chan combinations
size(C,1); %total number of combinations
%pull out channel combination
data1 = monkey(1).day(i).(resp{j}).(C{k,1});
data2 = monkey(1).day(i).(resp{j}).(C{k,2});

%compute coherence 
[C,phi,S12,S1,S2,t,f] = cohgramc(data1',data2',movingwin,params);

%now plot, need baseline corrected dB code


[P,T,F]=mtspecgramc(Correct.area8B,movingwin,params);
subplot(2,6,2)
imagesc(T,F,10*log10(P)) %Plot power in dB
axis xy; title('Frontal 8B'); xlabel('Time(s)'); ylabel('Freq (Hz)'); colormap jet; colorbar;

%other plot
plot_matrix(P,T,F); xlabel([]); % plot spectrogram
%caxis([8 28]); 
colorbar;
set(gca,'FontName','Times New Roman','Fontsize', 14);
title({['LFP 1,  W=' num2str(params.tapers(1)/movingwin(1)) 'Hz']; ['moving window = ' num2str(movingwin(1)) 's, step = ' num2str(movingwin(2)) 's']});
ylabel('frequency Hz');

%contour plot
contourf(times2save,frex,squeeze(mi(i,:,:))-repmat(mean(mi(i,:,baseidx(1):baseidx(2)),3)',1,length(times2save)),40,'linecolor','none')
set(gca,'clim',[-.075 .075],'yscale','log','ytick',round(logspace(log10(frex(1)),log10(frex(end)),6)))
xlabel('Time (ms)'), ylabel('Frequency (Hz)'), title('MI Based on Power')
c = colorbar; set(get(c,'label'),'string','MI (baseline subtracted)');    