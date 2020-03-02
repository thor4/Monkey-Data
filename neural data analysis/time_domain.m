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
    % avg over trials & convert to �V (1V = 10^6�V = 1,000,000�V) for ERP
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
            title(erptitle); xlabel('Time (ms)'); ylabel('Voltage (�V)'); 
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
        title(erptitle); xlabel('Time (ms)'); ylabel('Voltage (�V)'); 
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

