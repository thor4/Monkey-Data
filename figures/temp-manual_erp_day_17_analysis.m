%% Segment the trials

% navigate to monkey B day 17 (suspicious):
% OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\betty\090917\session01

load('trial_info.mat')

idx = 0; %init counter
for k=1:trial_info.numTrials %parse all trials for day
%stable not specified, give stable perf and transition trials
    if (trial_info.good_trials(k) == 1) && ...%artifacts/none 0/1
            (trial_info.BehResp(k) == 0) && ... %correct(1)/incorrect(0)
            (trial_info.rule(k) == 1) %identity(1)/location(2)
        idx = idx + 1;  
        trial_lfp = sprintf('betty09091701.%04d.mat',k);
        load(trial_lfp,'lfp_data');
        base = lfp_data(:,floor(trial_info.CueOnset(k))-504:floor(trial_info.CueOnset(k))-1);
        sample = lfp_data(:,floor(trial_info.CueOnset(k)):floor(trial_info.CueOnset(k))+504);
        delay = lfp_data(:,floor(trial_info.CueOffset(k)):floor(trial_info.CueOffset(k))+810);
        match = lfp_data(:,floor(trial_info.MatchOnset(k)):floor(trial_info.MatchOnset(k))+199);
        lfp(:,:,idx) = cat(2,base,sample,delay,match);
    end
end

%save lfp as lfpC or lfpI dep on behResp

%504 samples in baseline
%505 samples in cue
%811 samples in delay
%200 samples in match
%2020 total samples across all chans
%gives all trials cut up and stitched together, ready for ERPing

%% ERP 

%init variables
time_x  = -504:size(lfpC,2)-505; % time, from -504ms baseline
triggers = [0 505 1316]; %epoch switches base/sample, sample/delay, delay/match

% avg over trials & convert to µV (1V = 10^6µV = 1,000,000µV) for ERP
correct = mean(mean(lfpC,3),1) .* 1e6; %1st mean across trials 2nd mean across chans
incorrect = mean(mean(lfpI,3),1) .* 1e6; %1st mean across trials 2nd mean across chans
erptitle = sprintf('%d %s & %d %s Trial-averaged Monkey B Day 17 / 24 All Areas & Chans',...
        size(lfpC,3),"Correct",size(lfpI,3),"Incorrect");
figure(1), clf
%add in the (1:end-73) to ensure only 200ms of match period shows
erpD = plot(time_x,correct-correct(1),time_x,incorrect-incorrect(1),':', 'LineWidth', 2);
set(gca,'box','off','Xlim',[time_x(1);time_x(end)]);
y1 = get(gca,'ylim'); hold on
epochs = plot([triggers(1) triggers(1)],y1,'--', ...
[triggers(2) triggers(2)],y1,'--',[triggers(3) triggers(3)],y1,'--'); 
epochs(1).Color = [0.5 0.5 0.5]; epochs(2).Color = [0.5 0.5 0.5];
epochs(3).Color = [0.5 0.5 0.5];
title(erptitle); xlabel('Time (ms)'); ylabel('Voltage (µV)'); 
text(time_x(1)+100,y1(2)-1,'baseline','FontSize',14); 
text(triggers(1)+100,y1(2)-1,'sample','FontSize',14);
text(triggers(2)+100,y1(2)-1,'delay','FontSize',14); 
text(triggers(3)+75,y1(2)-1,'match','FontSize',14);
ax=gca; ax.FontSize = 18;

%it's the same


