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
%data.M1goodR1(1).d060328 %correct then .lfp, .erp or .areas
%data.M2goodR1(2).d090618 %incorrect then .lfp, .erp or .areas

%make loop to go through all days+chans both monkeys
%see about showing day's average at end of each day as well in loop
%see about saving as a movie, not a gif
%gif is here: https://www.mathworks.com/matlabcentral/answers/422908-animation-using-plot-inside-for-loop

time  = -504:size(dataM2goodCorR1.d090709.erp(2,:),2)-505; % time, from -504ms baseline
triggers = [0 505 1316]; %epoch switches base/sample, sample/delay, delay/match
monkeys = fieldnames(data)'; %extract both monkeys (2)
% monkey = 'M2goodR1';

fileID = fopen('total_erps.txt','a');
for monkey=monkeys
    behresponses = size(data.(monkey{:}),2); %extract # of responses (2)
    for response=1:behresponses
        alldays = fieldnames(data.(monkey{:})(response))'; %extract all days
        for dday=alldays
            allchans = size(data.(monkey{:})(response).(dday{:}).erp,1);
            for chan=1:allchans
                erpN = data.(monkey{:})(response).(dday{:}).erp(chan,:);
                if (response==1); resp = "Correct";
                else; resp = "Incorrect"; end
                erptitle = sprintf('%d %s Trial-averaged Monkey %d Day %d / %d Area %s Chan %d / %d',...
                    size(data.(monkey{:})(response).(dday{:}).lfp,3),resp,...
                    find(strcmp(monkeys,monkey{:})),...
                    find(strcmp(alldays,dday{:})),...
                    size(fieldnames(data.(monkey{:})(response)),1),...
                    data.(monkey{:})(response).(dday{:}).areas{chan},chan,...
                    size(data.(monkey{:})(response).(dday{:}).erp,1));
                fprintf(fileID,'%s\n',erptitle);
%                 pause(1) %pause 1 second
            end
        end
    end
end
fclose(fileID);



for dataset=data
    alldays = fieldnames((dataset{:}));
    for dday=alldays
        totalchans = size((dataset{:}).(dday).lfp,1);
        for chan=1:totalchans
            
    (i{:})
chan=3;
figure(2), clf
%correct #6DB3A5: [0.4941 0.7294 0.8000], incorrect #C9778F: [0.7882 0.4667 0.5608] 
% newcolors = {'[0.4941 0.7294 0.8000]','[0.7882 0.4667 0.5608]'};
% colororder(newcolors)
%shift correct and incorrect traces to start at 0
correct = dataM2goodCorR1.(dday).erp(chan,:) + (-dataM2goodCorR1.(dday).erp(chan,1));
incorrect = dataM2goodIncR1.(dday).erp(chan,:) + (-dataM2goodIncR1.(dday).erp(chan,1));
%day ERPs
correctD = mean(dataM2goodCorR1.(dday).erp,1) + (-mean(dataM2goodCorR1.(dday).erp(:,1)));
incorrectD = mean(dataM2goodIncR1.(dday).erp,1) + (-mean(dataM2goodIncR1.(dday).erp(:,1)));
erp = plot(time,correctD,time,incorrectD,':', 'LineWidth', 2);
% chanERP = plot(time,[dataM2goodCorR1.d090709.erp(1,:)]);
set(gca,'box','off','Xlim',[time(1);time(end)]);
y1 = get(gca,'ylim'); hold on
epochs = plot([triggers(1) triggers(1)],y1,'--', ...
    [triggers(2) triggers(2)],y1,'--',[triggers(3) triggers(3)],y1,'--'); 
epochs(1).Color = [0.5 0.5 0.5]; epochs(2).Color = [0.5 0.5 0.5];
epochs(3).Color = [0.5 0.5 0.5];
erptitleD = sprintf('%d Correct & %d Incorrect Trial-averaged Monkey %d Day %d / %d All areas/chans',...
    size(dataM2goodCorR1.(dday).lfp,3),size(dataM2goodIncR1.(dday).lfp,3),...
    monk,dayN,size(fieldnames(dataM2goodCorR1),1)); %day
erptitle = sprintf('%d Correct & %d Incorrect Trial-averaged Monkey %d Day %d / %d Area %s Chan %d / %d',...
    size(dataM2goodCorR1.(dday).lfp,3),size(dataM2goodIncR1.(dday).lfp,3),...
    monk,dayN,size(fieldnames(dataM2goodCorR1),1),...
    dataM2goodCorR1.(dday).areas{chan},chan,size(dataM2goodCorR1.(dday).erp,1));
title(erptitleD)
xlabel('Time (ms)'); ylabel('Voltage (�V)'); 
text(time(1)+100,y1(2)-1,'baseline'); text(triggers(1)+100,y1(2)-1,'sample');
text(triggers(2)+100,y1(2)-1,'delay'); text(triggers(3)+100,y1(2)-1,'match');
pause %pauses until user presses key
% chanERP(1).Color = [0.75 0.75 0.75]; chanERP(2).Color = [0.75 0.75 0.75]; %light gray

% 1-504: baseline  (505) graph: -504-1
% 505-1009: sample (505) graph: 0-504
% 1010-1820: delay (811) graph: 505-1315
% 1821-2094: match (274) graph: 1316-1589