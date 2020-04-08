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
trial_info_path = strcat(path,'%s\\%s\\%s\\trial_info.mat'); %build trial_info path
trial_infoN = sprintf(trial_info_path,monkey,days_betty{day},"session01"); %create full path to trial_info.mat
load(trial_infoN,'trial_info'); %load trial_info for day's trials

%extract entire trial lfp for good/correct/rule1
[dayN,areasN] = extractDay(path,monkey,days_betty{day},1,2,1,1,'entire'); 
chans = length(areasN)+3; %total number of channels + key + 2 eye chans
length(fieldnames(dayN)) %total number of good/correct/rule1 trials
trialNames = fieldnames(dayN);
trials = str2double(strip(trials,'left','t')); %identify trials used
currTrial = trials(1); %identify current trial

% define time, from -ms baseline prior to stimulus onset
time  = 0-trial_info.CueOnset(currTrial):size(dayN.(trialNames{1})(1,:),2)-trial_info.CueOnset(currTrial)-1;
triggers = [0 ... %epoch switch base/sample
    trial_info.CueOffset(currTrial)-trial_info.CueOnset(currTrial) ... %sample/delay
    trial_info.MatchOnset(currTrial)-trial_info.CueOnset(currTrial)]; %delay/match

%load eye channels
trial_path = strcat(path,'%s\\%s\\%s\\%s%s%s.%04d.mat'); %build trial data path
trial_data = sprintf(trial_path,monkey,days_betty{day},"session01",monkey,days_betty{day},'01',currTrial);
load(trial_data,'vertical_eye'); load(trial_data,'horizontal_eye');
%subsample eye channels to get them to 10kHz
downH = downsample(horizontal_eye,30); downV = downsample(vertical_eye,30);
h_eye = [downH downH(end)]; v_eye = [downV downV(end)]; %repeat last measurement to = lfp

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
export_fig eye_movements_trial_302.png -transparent % no background

%next: work on subplots. figure out scale of eye movements to add those
%subplots as tiles

hf=figure; clf
t = tiledlayout(chans,1); %setup tile for all subplots
cmap = colormap; %get current colormap
allColors = cmap(chans:chans:chans*chans,:); %split colormap into diff colors per chan
for subplotN=1:chans-3
    nexttile
    %plot and convert to µV (1V = 10^6µV = 1,000,000µV)
    plot(time,dayN.(trialNames{1})(subplotN,:) .* 1e6,'LineWidth',2,'Color',allColors(subplotN,:)) 
end

%group the raw plots
axes = findobj(gcf,'type','axes'); %aggregate all axes from all tiles (not subplots)
linkaxes(axes,'xy') %link all tiles so axes are on same scale
y1 = get(gca,'ylim'); %get y-axis limits
%set x-axis scale + turn off box + xtick labels and rename y labels
set(axes,'Xlim',[time(1);time(end)],'Visible','off'); 

nexttile %add key for scale
line([time(end-201) time(end-201)],[y1(1) y1(1)+500],'LineWidth',2,'Color','k'); %y line
line([time(end-201) time(end-201)+200],[y1(1) y1(1)],'LineWidth',2,'Color','k'); %x line
set(gca,'ylim',y1,'xlim',[time(1);time(end)],'Visible','off');
text(time(end-201)+50,y1(1)-175,'200 ms','FontSize',12,'HorizontalAlignment','left')
text(time(end-201)-10,y1(1)+75,'0 mV','FontSize',12,'HorizontalAlignment','right')
text(time(end-201)-10,y1(1)+550,'0.5 mV','FontSize',12,'HorizontalAlignment','right')

%minimize the spacing around the perimeter of the layout & around each tile
t.Padding = 'compact'; t.TileSpacing = 'none';

%add channel label to each tile
for axN=1:chans
    text(axes(axN),time(1)-25,(y1(1)+y1(2))/2,areasN{14-axN},'FontSize',14,'HorizontalAlignment','right')
end
%add epoch labels
text(axes(end),triggers(1)-300,y1(2)+150,'Baseline','FontSize',16,'HorizontalAlignment','right')
text(axes(end),triggers(2)-200,y1(2)+150,'Sample','FontSize',16,'HorizontalAlignment','right')
text(axes(end),triggers(2)+500,y1(2)+150,'Delay','FontSize',16,'HorizontalAlignment','right')
text(axes(end),triggers(3)+100,y1(2)+150,'Match','FontSize',16,'HorizontalAlignment','left')

%add epoch lines across multiple subplots, bottom to top and top to bottom
set(gca,'Clipping','Off'), set(axes(end),'Clipping','Off'),  %turn off clipping in top & bottom plots
h1 = line([triggers(1) triggers(1)],[y1(1) (y1(2)-y1(1))*(chans+3)]); %primitive line, bottom up 
h12 = line(axes(end),[triggers(1) triggers(1)],[y1(2)+150 -(y1(2)-y1(1))*(chans+3)]); %primitive line, top down
h2 = line([triggers(2) triggers(2)],[y1(1) (y1(2)-y1(1))*(chans+3)]); %primitive line, bottom up 
h22 = line(axes(end),[triggers(2) triggers(2)],[y1(2)+150 -(y1(2)-y1(1))*(chans+3)]); %primitive line, top down
h3 = line([triggers(3) triggers(3)],[y1(1) (y1(2)-y1(1))*(chans+3)]); %primitive line, bottom up 
h32 = line(axes(end),[triggers(3) triggers(3)],[y1(2)+150 -(y1(2)-y1(1))*(chans+3)]); %primitive line, top down
%reset ylim so plot doesn't resize to accommodate line
ylim([y1(1) y1(2)]); ylim(axes(end),[y1(1) y1(2)])
set(h1,'LineWidth',2,'Color','k'); set(h12,'LineWidth',2,'Color','k')
set(h2,'LineWidth',2,'Color','k'); set(h22,'LineWidth',2,'Color','k')
set(h3,'LineWidth',2,'Color','k'); set(h32,'LineWidth',2,'Color','k')

% title('Monkey 2, Day 17, Good, Correct, Rule 1, Trial 302')

export_fig monkey_2_day_17_good_correct_rule_1_trial_302.png -transparent % no background

hf.Color='w'; %Set background color of figure window

