%% Method Figure - LFP traces

%init variables
path = 'H:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\';
monkey='betty'; %only betty
days_betty = { '090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001' };
day = 17; epochs = {'Baseline','Sample','Delay','Match'}; %identify day + epochs
trial_info_path = strcat(path,'%s\\%s\\%s\\trial_info.mat'); %build trial_info path
trial_infoN = sprintf(trial_info_path,monkey,days_betty{day},"session01"); %create full path to trial_info.mat
load(trial_infoN,'trial_info'); %load trial_info for day's trials

%extract entire trial lfp for good/correct/rule1
[dayN,areasN] = extractDay(path,monkey,days_betty{day},1,2,1,1,'entire'); 
chans = length(areasN)+3; %total number of channels + key + 2 eye chans
% length(fieldnames(dayN)) %total number of good/correct/rule1 trials
trialNames = fieldnames(dayN);
tcurrTrial = char(randsample(trialNames,1));
currTrial = str2double(strip(tcurrTrial,'left','t')); %identify trials used

% define time, from -ms baseline prior to stimulus onset
time  = 0-trial_info.CueOnset(currTrial):size(dayN.(tcurrTrial)(1,:),2)-trial_info.CueOnset(currTrial)-1;
triggers = [0 ... %epoch switch base/sample
    trial_info.CueOffset(currTrial)-trial_info.CueOnset(currTrial) ... %sample/delay
    trial_info.MatchOnset(currTrial)-trial_info.CueOnset(currTrial)]; %delay/match

%load eye channels
trial_path = strcat(path,'%s\\%s\\%s\\%s%s%s.%04d.mat'); %build trial data path
trial_data = sprintf(trial_path,monkey,days_betty{day},"session01",monkey,days_betty{day},'01',currTrial);
load(trial_data,'vertical_eye'); load(trial_data,'horizontal_eye');
%subsample eye channels to get them to 10kHz
downH = downsample(horizontal_eye,30); downV = downsample(vertical_eye,30);
%repeat last eye measurement until vector length equals total time
fillerH = repmat(downH(end),1,length(time)-length(downH));
fillerV = repmat(downV(end),1,length(time)-length(downV));
h_eye = [downH fillerH]; v_eye = [downV fillerV]; %repeat last measurement to = lfp

 %subplot version
clf
cmap = colormap; %get current colormap
allColors = cmap(chans:chans:chans*chans,:); %split colormap into diff colors per chan
for subplotN=1:chans
    if subplotN<chans-2 %plot raw lfp traces
        hSub(subplotN) = subplot(chans,1,subplotN);
        x = [find(time==-500),find(time==triggers(3)+200)]; %set cutoff 500ms before sample & 200ms after match
        y = dayN.(tcurrTrial)(subplotN,:).* 1e6; %convert to µV (1V = 10^6µV = 1,000,000µV)
        y = y(trial_info.CueOnset(currTrial)-500:trial_info.MatchOnset(currTrial)+200); %match to time
        plot(time(x(1):x(2)),y,'LineWidth',2,'Color',allColors(subplotN,:)) 
        line([time(x(2)+5) time(x(2)+5)],[min(y) max(y)],'LineWidth',2,'Color','k'); %y line
        scale = sprintf('%0.2f mV',(max(y)-min(y))/1000); %scale of y-axis
        text(time(x(2)+10),(max(y)-min(y))/2*0.5,scale,'FontSize',12,'HorizontalAlignment','left')
        set(gca,'xlim',[time(x(1));time(x(2)+5)],'ylim',[min(y) max(y)],'Visible','off');
        text(time(x(1))-10,(max(y)-min(y))/2*0.5,areasN{subplotN},'FontSize',14,'HorizontalAlignment','right')
    elseif subplotN==chans-2 %plot vertical eye movements
        hSub(subplotN) = subplot(chans,1,subplotN);
        vem_y = v_eye(trial_info.CueOnset(currTrial)-500:trial_info.MatchOnset(currTrial)+200); %match to time
        % find exponent automatically and divide by it when plotting
        vem_y_scaled = floor(log10(max(vem_y)));
        plot(time(x(1):x(2)),vem_y/10^vem_y_scaled,'LineWidth',2,'color','b') %vertical eye movements
        temp_ylim = ylim;
        text(time(x(1))-10,(temp_ylim(1)+temp_ylim(2))/2,'V','FontSize',14,'HorizontalAlignment','right')
        set(gca,'xlim',[time(x(1));time(x(2))],'Visible','off');      
    elseif subplotN==chans-1 %plot horizontal eye movements
        hSub(subplotN) = subplot(chans,1,subplotN);
        hem_y = h_eye(trial_info.CueOnset(currTrial)-500:trial_info.MatchOnset(currTrial)+200); %match to time
        % find exponent automatically and divide by it when plotting
        hem_y_scaled = floor(log10(max(hem_y)));
        plot(time(x(1):x(2)),hem_y/10^hem_y_scaled,'LineWidth',2,'color','r') %horizontal eye movements
        temp_ylim = ylim;
        text(time(x(1))-10,(temp_ylim(1)+temp_ylim(2))/2,'H','FontSize',14,'HorizontalAlignment','right')
        set(gca,'xlim',[time(x(1));time(x(2))],'Visible','off');      
    else %draw key
        hSub(subplotN) = subplot(chans,1,subplotN);
        line([time(x(2)-201) time(x(2)-201)+200],[1 1],'LineWidth',2,'Color','k'); %x line
        text(time(x(2)-201)+50,0.5,'200 ms','FontSize',12,'HorizontalAlignment','left')
        set(gca,'xlim',[time(x(1));time(x(2))],'Visible','off');      
    end
end

%figure position: [left bottom width height]
allsubPos=reshape([hSub.Position],4,[])'; %positions for each subplot allsubPos(1:end)=top:bottom
allsubPos(:,4)=repmat(1/(numel(hSub)+1),numel(hSub),1); %*height* assign based on total plots
for i=1:numel(hSub)
    allsubPos(i,2)=allsubPos(i,4)*(numel(hSub)-i);
    if (allsubPos(i,2) < 0)
        disp('subplot ',i,' is messed up');
    end
end
for subPos=1:numel(hSub) %set each subplot to new position... 
    hSub(subPos).Position=allsubPos(subPos,:); %top grows down, bottom grows up
end

 %draw epoch lines 
pos(1:3) = hSub(1).Position(1:3); %get left/bottom/width of bottom axis
pos(4) = sum(hSub(end).Position([2 4])); %top of top axis
xlims = hSub(2).XLim; %get limits of x-axis
%express cumsum of x-coords from pos as a function of x-axis limits, then
%identify where the lines should be drawn on this scale
line_locations_norm = interp1(xlims, cumsum(pos([1 3])), triggers);
for i=1:numel(line_locations_norm) %draw the lines
    annotation('line', ...
        [line_locations_norm(i) line_locations_norm(i)], ...
        [pos(2)+hSub(1).Position(4) pos(4)], ... %start from top of bottom plot
        'Color', 'k', ...
        'LineWidth', 2);
end

%add epoch labels
text(hSub(1),triggers(1)-200,hSub(1).YLim(2)*1.2,epochs{1},...
    'FontSize',16,'HorizontalAlignment','right','color',[0.5 0.5 0.5])
text(hSub(1),triggers(1)+200,hSub(1).YLim(2)*1.2,epochs{2},...
    'FontSize',16,'HorizontalAlignment','right','color',[0.5 0.5 0.5])
text(hSub(1),triggers(3)-500,hSub(1).YLim(2)*1.2,epochs{3},...
    'FontSize',16,'HorizontalAlignment','right','color',[0.5 0.5 0.5])
text(hSub(1),triggers(3)+75,hSub(1).YLim(2)*1.2,epochs{4},...
    'FontSize',16,'HorizontalAlignment','left','color',[0.5 0.5 0.5])

trial_title = sprintf('Monk2 D%d Good Cor Rule1 Trial %d / %d',...
                day,trialN,length(trials));
sgtitle(trial_title,'FontSize',18);  

fileName = sprintf('mB_day_%d_trial_%d_lfps.pdf',day,currTrial);
fullFileName = fullfile(pwd, fileName); %create png filename in pwd
export_fig(fullFileName); %save fig
