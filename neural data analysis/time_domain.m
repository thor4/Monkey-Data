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

%% Step 1: Extract good, cor+inc rule 1 trials from both monkeys %%
%include both stable & unstable (monkey maybe not sure which rule in play)
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
dayN = 13; %day number
day = alldays{dayN}; %assign day
dday = append('d',day); %setup day for test indexing
i=0; %init counter

%change var name per monkey/resp combo (4x)
for lp=alldays %cycle through all days
    i = i+1;
    dayy = append('d',lp{:});
    [dataM2goodCorR1.(dayy).lfp,dataM2goodCorR1.(dayy).areas] = extractDay(path,monkey,lp{:},1,2,1,1,'all');
    % avg over trials & convert to µV (1V = 10^6µV = 1,000,000µV) for ERP
    dataM2goodCorR1.(dayy).erp = mean(dataM2goodCorR1.(dayy).lfp,3) .* 1e6;
end

%test
dataM2goodCorR1.(dayy).erp == mean(dataM2goodCorR1.(dayy).lfp,3) .* 1e6

save('time_domain-m1','dataM1goodCorR1','dataM1goodIncR1')
save('time_domain-m2','dataM2goodCorR1','dataM2goodIncR1','-v7.3') %v7.3 since filesize so big
%home pc: D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data

%% Step 2: Visualize ERPs for each chan from every day in same fig

load('time_domain-m1.mat')
load('time_domain-m2.mat')

%package into data variable 1: correct & 2: incorrect
data.M1goodR1(1) = dataM1goodCorR1; data.M1goodR1(2) = dataM1goodIncR1; 
data.M2goodR1(1) = dataM2goodCorR1; data.M2goodR1(2) = dataM2goodIncR1; 

figure(1), clf %open fig and maximize to prepare for viewing

%init variables
time  = -504:size(dataM2goodCorR1.d090709.erp(2,:),2)-505; % time, from -504ms baseline
triggers = [0 505 1316]; %epoch switches base/sample, sample/delay, delay/match
monkeys = fieldnames(data)'; %extract both monkeys (2)
%init video
erpVid = VideoWriter('erpVid'); %open video file
erpVid.FrameRate = 5;  %can adjust this, 5 - 10 seems to work
open(erpVid)

% fileID = fopen('total_erps.txt','a'); %record all erp titles
for monkey=monkeys
%     behresponses = size(data.(monkey{:}),2); %extract # of responses (2)
    alldays = fieldnames(data.(monkey{:})(1))'; %extract all days
    for dday=alldays
        allchans = size(data.(monkey{:})(1).(dday{:}).erp,1);
        for chan=1:allchans
            erptitle = sprintf('%d %s & %d %s Trial-averaged Monkey %d Day %d / %d Area %s Chan %d / %d',...
                size(data.(monkey{:})(1).(dday{:}).lfp,3),"Correct",...
                size(data.(monkey{:})(2).(dday{:}).lfp,3),"Incorrect",...
                find(strcmp(monkeys,monkey{:})),...
                find(strcmp(alldays,dday{:})),...
                size(fieldnames(data.(monkey{:})(1)),1),...
                data.(monkey{:})(1).(dday{:}).areas{chan},chan,...
                size(data.(monkey{:})(1).(dday{:}).erp,1));
%                 fprintf(fileID,'%s\n',erptitle); %record erp title
            correct = data.(monkey{:})(1).(dday{:}).erp(chan,:) + (-data.(monkey{:})(1).(dday{:}).erp(chan,1));
            incorrect = data.(monkey{:})(2).(dday{:}).erp(chan,:) + (-data.(monkey{:})(2).(dday{:}).erp(chan,1));
            figure(1), clf
            erp = plot(time,correct,time,incorrect,':', 'LineWidth', 2);
            set(gca,'box','off','Xlim',[time(1);time(end)]);
            y1 = get(gca,'ylim'); hold on
            epochs = plot([triggers(1) triggers(1)],y1,'--', ...
            [triggers(2) triggers(2)],y1,'--',[triggers(3) triggers(3)],y1,'--'); 
            epochs(1).Color = [0.5 0.5 0.5]; epochs(2).Color = [0.5 0.5 0.5];
            epochs(3).Color = [0.5 0.5 0.5];
            title(erptitle); xlabel('Time (ms)'); ylabel('Voltage (µV)'); 
            text(time(1)+100,y1(2)-1,'baseline'); text(triggers(1)+100,y1(2)-1,'sample');
            text(triggers(2)+100,y1(2)-1,'delay'); text(triggers(3)+100,y1(2)-1,'match');
%             pause(1) %pause 1 second
            pause(0.01) %pause to grab frame
            frame = getframe(gcf); %get frame
            writeVideo(erpVid, frame); %add frame to vid
        end
        %make day-level erp
        correct = mean(data.(monkey{:})(1).(dday{:}).erp,1) + (-mean(data.(monkey{:})(1).(dday{:}).erp(:,1)));
        incorrect = mean(data.(monkey{:})(2).(dday{:}).erp,1) + (-mean(data.(monkey{:})(2).(dday{:}).erp(:,1)));
        erptitle = sprintf('%d %s & %d %s Trial-averaged Monkey %d Day %d / %d All Areas & Chans',...
                size(data.(monkey{:})(1).(dday{:}).lfp,3),"Correct",...
                size(data.(monkey{:})(2).(dday{:}).lfp,3),"Incorrect",...
                find(strcmp(monkeys,monkey{:})),...
                find(strcmp(alldays,dday{:})),...
                size(fieldnames(data.(monkey{:})(1)),1));
        figure(1), clf
        erpD = plot(time,correct,time,incorrect,':', 'LineWidth', 2);
        set(gca,'box','off','Xlim',[time(1);time(end)]);
        y1 = get(gca,'ylim'); hold on
        epochs = plot([triggers(1) triggers(1)],y1,'--', ...
        [triggers(2) triggers(2)],y1,'--',[triggers(3) triggers(3)],y1,'--'); 
        epochs(1).Color = [0.5 0.5 0.5]; epochs(2).Color = [0.5 0.5 0.5];
        epochs(3).Color = [0.5 0.5 0.5];
        title(erptitle); xlabel('Time (ms)'); ylabel('Voltage (µV)'); 
        text(time(1)+100,y1(2)-1,'baseline'); text(triggers(1)+100,y1(2)-1,'sample');
        text(triggers(2)+100,y1(2)-1,'delay'); text(triggers(3)+100,y1(2)-1,'match');
%         pause(1) %pause 1 second
        pause(0.01) %pause to grab frame
        frame = getframe(gcf); %get frame
        writeVideo(erpVid, frame); %add frame to vid
    end
end
% fclose(fileID);
close(erpVid)

%saved total_erps.txt in D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\results
%saved erpVid.avi in same

%% Step 3 Explore Days 17-24 from Monkey 2 where there is clear separation
% between correct and incorrect trials, esp during delay that doesn't
% appear physiological in origin. Focus on correct+good trials only

% 30Khz = 1/30,000 sec %eye tracking signal
% 1Khz = 1/1,000 sec %lfp signal
% 
% 30,000 / 30 = 1000
% [1:30:30001] % subsample the 30Khz signal to get it down to 1Khz

%init variables
path = 'D:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\';
monkey='betty'; %only betty
days_betty = { '090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001' };
day = 17; %look at 17-24
epochs = {'Baseline','Sample','Delay','Match'};
trial = 450; %look at second trial for test
trial_info_path = strcat(path,'%s\\%s\\%s\\trial_info.mat'); %build trial_info path
trial_infoN = sprintf(trial_info_path,monkey,days_betty{day},"session01"); %create full path to trial_info.mat
load(trial_infoN,'trial_info'); %load trial_info for day's trials

%extract entire trial lfp for good/correct/rule1
[dayN,areasN] = extractDay(path,monkey,days_betty{day},1,2,1,1,'entire'); 
chans = length(areasN)+3; %total number of channels + key + 2 eye chans
% length(fieldnames(dayN)) %total number of good/correct/rule1 trials
trialNames = fieldnames(dayN);
trials = str2double(strip(trialNames,'left','t')); %identify trials used
currTrial = trials(trial); %identify current trial

% define time, from -ms baseline prior to stimulus onset
time  = 0-trial_info.CueOnset(currTrial):size(dayN.(trialNames{trial})(1,:),2)-trial_info.CueOnset(currTrial)-1;
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

%temp figure for Charlie to see eye movements and comment on scale
figure, clf
tiledlayout(2,1); %setup tile for all subplots
nexttile
plot(time,h_eye,'LineWidth',2)
title('Horizontal eye movements')
nexttile
plot(time,v_eye,'LineWidth',2)
xlabel('Time (ms)')
title('Vertical eye movements')
% export_fig eye_movements_trial_302.png -transparent % no background

%next: see if Charlie wants to keep the y-axis key on the right

% % legacy stackedplot attempt
% % hf=figure; clf
% % cmap = colormap; %get current colormap
% % allColors = cmap(chans:chans:chans*chans,:); %split colormap into diff colors per chan
% % ylimit = [-100,100]; %set y-limit for scaling
% y = zeros(length([find(time==-500):find(time==triggers(3)+200)]),chans-1); %init
% xmarks = [find(time==-500),find(time==triggers(3)+200)]; %set cutoff 500ms before sample & 200ms after match
% x = time(xmarks(1):xmarks(2));
% for subplotN=1:chans-3 %plot first set of time-series
%     temp = dayN.(trialNames{trial})(subplotN,:).* 1e6; %convert to µV (1V = 10^6µV = 1,000,000µV)
%     y(:,subplotN) = temp(trial_info.CueOnset(currTrial)-500:trial_info.MatchOnset(currTrial)+200); %match to time
% end
% % add eye movement plots: vert then horiz
% y(:,subplotN+1) = v_eye(trial_info.CueOnset(currTrial)-500:trial_info.MatchOnset(currTrial)+200); %match to time
% y(:,subplotN+2) = h_eye(trial_info.CueOnset(currTrial)-500:trial_info.MatchOnset(currTrial)+200); %match to time
% % add key
% y(time(xmarks(2)-201):time(xmarks(2)-201)+200,subplotN+3) = ones(1,length(time(xmarks(2)-201):time(xmarks(2)-201)+200)); 
% key=y(:,subplotN+3); key(key==0)=nan; y(:,subplotN+3)=key; %replace 0's with NaN
% sp = stackedplot(x,y,'LineWidth',2,'FontSize',12); %plot them all
% spPos = sp.Position; %get positions to create subplot version

% %not sure about removing the white background, tickmarks and alter the tick
% %labels. also need to add epoch text
% sp.LineProperties(1:16) %for changing colors
% 
% 
% sp.DisplayLabels = [areasN{:} {'V'} {'H'} {''}]; %add labels

% export_fig temp_disp.png

%legacy tiledplot version
% hf=figure; clf
% t = tiledlayout(chans,1); %setup tile for all subplots
% cmap = colormap; %get current colormap
% allColors = cmap(chans:chans:chans*chans,:); %split colormap into diff colors per chan
% ylimit = [-100,100]; %set y-limit for scaling
% for subplotN=1:chans-3
%     nexttile
%     x = [find(time==-500),find(time==triggers(3)+200)]; %set cutoff 500ms before sample & 200ms after match
%     y = dayN.(trialNames{trial})(subplotN,:).* 1e6; %convert to µV (1V = 10^6µV = 1,000,000µV)
%     y = y(trial_info.CueOnset(currTrial)-500:trial_info.MatchOnset(currTrial)+200); %match to time
%     plot(time(x(1):x(2)),y,'LineWidth',2,'Color',allColors(subplotN,:)) 
% %     ylim(ylimit)
% end

%subplot version
hf=figure; clf
cmap = colormap; %get current colormap
allColors = cmap(chans:chans:chans*chans,:); %split colormap into diff colors per chan
% ylimit = [-100,100]; %set y-limit for scaling
%legacy subplot version
for subplotN=1:chans
%     nexttile
    if subplotN<chans-2 %plot raw lfp traces
        hSub(subplotN) = subplot(chans,1,subplotN);
        x = [find(time==-500),find(time==triggers(3)+200)]; %set cutoff 500ms before sample & 200ms after match
        y = dayN.(trialNames{trial})(subplotN,:).* 1e6; %convert to µV (1V = 10^6µV = 1,000,000µV)
        y = y(trial_info.CueOnset(currTrial)-500:trial_info.MatchOnset(currTrial)+200); %match to time
        plot(time(x(1):x(2)),y,'LineWidth',2,'Color',allColors(subplotN,:)) 
        line([time(x(2)+5) time(x(2)+5)],[min(y) max(y)],'LineWidth',2,'Color','k'); %y line
        scale = sprintf('%0.2f mV',(max(y)-min(y))/1000); %scale of y-axis
        text(time(x(2)+10),(max(y)-min(y))/2*0.5,scale,'FontSize',12,'HorizontalAlignment','left')
%         set(gca,'ylim',ylimit,'xlim',[time(x(1));time(x(2))],'Visible','off');      
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
%         line([time(x(2)-201) time(x(2)-201)],[ylimit(1) ylimit(1)+200],'LineWidth',2,'Color','k'); %y line
%         line([time(x(2)-201) time(x(2)-201)+200],[ylimit(1) ylimit(1)],'LineWidth',2,'Color','k'); %x line
        line([time(x(2)-201) time(x(2)-201)+200],[1 1],'LineWidth',2,'Color','k'); %x line
        text(time(x(2)-201)+50,0.5,'200 ms','FontSize',12,'HorizontalAlignment','left')
%         text(time(x(2)-201)-10,ylimit(1)+25,'0 mV','FontSize',12,'HorizontalAlignment','right')
%         text(time(x(2)-201)-10,ylimit(1)+175,'0.2 mV','FontSize',12,'HorizontalAlignment','right')
        set(gca,'xlim',[time(x(1));time(x(2))],'Visible','off');      
    end
end

allsubPos=reshape([hSub.Position],4,[])'; %positions for each subplot allsubPos(1:end)=top:bottom
% dHt=spPos(4)-sum(allsubPos(:,4)); %height of stackedplot - height of subplots
allsubPos(:,4)=repmat(1/(numel(hSub)+1),numel(hSub),1); %*height* assign based on total plots
for i=1:numel(hSub)
    allsubPos(i,2)=allsubPos(i,4)*(numel(hSub)-i);
%     allsubPos(end-i+1,2)=allsubPos(i,4)*(i-1);
    if (allsubPos(i,2) < 0)
        disp('subplot ',i,' is messed up');
    end
end
% AXposn(:,4)=AXposn(:,4)+dHt/numel(hSub); %allocate difference to grow each subplot
% AXposn(:,2)=AXposn(:,2)-dHt/numel(hSub); %lower the bottom of each subplot to make up for added height
for subPos=1:numel(hSub) %set each subplot to new position... 
    hSub(subPos).Position=allsubPos(subPos,:); %top grows down, bottom grows up
end

%group the raw plots 
% allaxes = findobj(gcf,'type','axes'); %aggregate all axes from all subplots
% allaxes(end) is top plot. allaxes(1) is bottom plot
% linkaxes(allaxes,'xy') %link all tiles so axes are on same scale

% nexttile %add eye movements
% vem_y = v_eye(trial_info.CueOnset(currTrial)-500:trial_info.MatchOnset(currTrial)+200); %match to time
% plot(time(x(1):x(2)),vem_y,'LineWidth',2,'color','b') %vertical eye movements
% text(time(x(1))-10,(ylimit(1)+ylimit(2))/2,'V','FontSize',14,'HorizontalAlignment','right')
% ylim(gca,[min(vem_y) max(vem_y)])
% nexttile 
% hem_y = h_eye(trial_info.CueOnset(currTrial)-500:trial_info.MatchOnset(currTrial)+200); %match to time
% plot(time(x(1):x(2)),hem_y,'LineWidth',2,'color','r') %horizontal eye movements
% text(time(x(1))-10,(ylimit(1)+ylimit(2))/2,'H','FontSize',14,'HorizontalAlignment','right')

%set x-axis scale + turn off box + xtick labels and rename y labels
% set(allaxes,'Xlim',[time(x(1)) time(x(2))],'Visible','off'); 
%minimize the spacing around the perimeter of the layout & around each tile
% t.Padding = 'compact'; t.TileSpacing = 'none';

% nexttile %add key for scale
% line([time(x(2)-201) time(x(2)-201)],[ylimit(1) ylimit(1)+200],'LineWidth',2,'Color','k'); %y line
% line([time(x(2)-201) time(x(2)-201)+200],[ylimit(1) ylimit(1)],'LineWidth',2,'Color','k'); %x line
% set(gca,'ylim',ylimit,'xlim',[time(x(1));time(x(2))],'Visible','off');
% text(time(x(2)-201)+50,ylimit(1)-75,'200 ms','FontSize',12,'HorizontalAlignment','left')
% text(time(x(2)-201)-10,ylimit(1)+25,'0 mV','FontSize',12,'HorizontalAlignment','right')
% text(time(x(2)-201)-10,ylimit(1)+175,'0.2 mV','FontSize',12,'HorizontalAlignment','right')

% ax.position = [left bottom width height]
% The left and bottom elements define the distance from the lower left corner
% of the container to the lower left corner of the position boundary
% The width and height elements are the position boundary dimensions.

%draw epoch lines 
pos(1:3) = hSub(1).Position(1:3); %get left/bottom/width of bottom axis
%top of top axis - bottom of bottom axis
% pos(4) = sum(allaxes(end).Position([2 4])) - allaxes(1).Position(2);
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
                day,trial,length(trials));
sgtitle(trial_title,'FontSize',18);  

% %add epoch labels
% text(allaxes(end),triggers(1)-200,ya1(1)+150,'Baseline','FontSize',16,'HorizontalAlignment','right')
% text(allaxes(end),triggers(2)-200,ya1(1)+150,'Sample','FontSize',16,'HorizontalAlignment','right')
% text(allaxes(end),triggers(2)+500,ya1(1)+150,'Delay','FontSize',16,'HorizontalAlignment','right')
% text(allaxes(end),triggers(3)+75,ya1(1)+150,'Match','FontSize',16,'HorizontalAlignment','left')

% 
% %draw lines delineating epochs
% %use external function located here:
% % https://github.com/michellehirsch/MATLAB-Dataspace-to-Figure-Units
% set(hf, 'currentaxes', allaxes(end));  %set current axis to top trace
% [xa1 ya1] = ds2nfu(triggers,repmat(ylimit(2),1,3)); %get fig-level plot points
% set(hf, 'currentaxes', allaxes(2));  %set current axis to bottom trace
% temp_ylim = ylim;
% [xa2 ya2] = ds2nfu(triggers,repmat(temp_ylim(1)-1,1,3)); %get fig-level plot points
% %draw lines: baseline/stimulus, stimulus/delay, delay/math
% annotation('line',[xa1(1) xa2(1)],[ya1(1) ya2(1)],'LineWidth',2,'color','k'); %draw lines
% annotation('line',[xa1(2) xa2(2)],[ya1(2) ya2(2)],'LineWidth',2,'color','k'); %draw lines
% annotation('line',[xa1(3) xa2(3)],[ya1(3) ya2(3)],'LineWidth',2,'color','k'); %draw lines

% %add channel label to each tile
% for axN=1:(chans-3)
%     text(allaxes(axN+3),time(x(1))-10,(ylimit(1)+ylimit(2))/2,areasN{length(areasN)+1-axN},'FontSize',14,'HorizontalAlignment','right')
% end


% title('Monkey 2, Day 17, Good, Correct, Rule 1, Trial 302')
hf.Color='w'; %Set background color of figure window

export_fig monkey_2_day_17_good_correct_rule_1_trial_932.png % white background



