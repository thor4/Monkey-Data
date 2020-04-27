currTrial = trials(trialN); %identify current trial
% define time, from -ms baseline prior to stimulus onset
time  = 0-trial_info.CueOnset(currTrial):size(dayN.(trialNames{trialN})(1,:),2)-trial_info.CueOnset(currTrial)-1;
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

% %temp figure for Charlie to see eye movements and comment on scale
% figure, clf
% tiledlayout(2,1); %setup tile for all subplots
% nexttile
% plot(time,h_eye,'LineWidth',2)
% title('Horizontal eye movements')
% nexttile
% plot(time,v_eye,'LineWidth',2)
% xlabel('Time (ms)')
% title('Vertical eye movements')
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
clf
cmap = colormap; %get current colormap
allColors = cmap(chans:chans:chans*chans,:); %split colormap into diff colors per chan
% ylimit = [-100,100]; %set y-limit for scaling
%legacy subplot version
for subplotN=1:chans
%     nexttile
    if subplotN<chans-2 %plot raw lfp traces
        hSub(subplotN) = subplot(chans,1,subplotN);
        x = [find(time==-500),find(time==triggers(3)+200)]; %set cutoff 500ms before sample & 200ms after match
        y = dayN.(trialNames{trialN})(subplotN,:).* 1e6; %convert to µV (1V = 10^6µV = 1,000,000µV)
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
                day,trialN,length(trials));
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