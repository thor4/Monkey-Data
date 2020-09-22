%% Method Figure - ERP traces

% the following files are located in:
% \OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data

load('time_domain-m1.mat')
load('time_domain-m2.mat')

%package into data variable 1: correct & 2: incorrect
data.M1goodR1(1) = dataM1goodCorR1; data.M1goodR1(2) = dataM1goodIncR1; 
data.M2goodR1(1) = dataM2goodCorR1; data.M2goodR1(2) = dataM2goodIncR1; 

figure(1), clf %open fig and maximize to prepare for viewing

%init variables
time_x  = -504:size(dataM2goodCorR1.d090709.erp(2,:),2)-505; % time, from -504ms baseline
triggers = [0 505 1316]; %epoch switches base/sample, sample/delay, delay/match
monkeys = fieldnames(data)'; %extract both monkeys (2)

%% Day ERPS across chans by Monkey by behResp

for monkey=monkeys
    alldays = fieldnames(data.(monkey{:})(1))'; %extract all days
    for dday=alldays
        allchans = size(data.(monkey{:})(1).(dday{:}).erp,1);
        %make day-level erp
        correct = mean(data.(monkey{:})(1).(dday{:}).erp,1) + (-mean(data.(monkey{:})(1).(dday{:}).erp(:,1)));
        incorrect = mean(data.(monkey{:})(2).(dday{:}).erp,1) + (-mean(data.(monkey{:})(2).(dday{:}).erp(:,1)));
        if find(strcmp(monkeys,monkey{:})) == 1 monk="A";
        else monk="B"; end
        erptitle = sprintf('%d %s & %d %s Trial-averaged Monkey %s Day %d / %d All Areas & Chans',...
                size(data.(monkey{:})(1).(dday{:}).lfp,3),"Correct",...
                size(data.(monkey{:})(2).(dday{:}).lfp,3),"Incorrect",...
                monk,...
                find(strcmp(alldays,dday{:})),...
                size(fieldnames(data.(monkey{:})(1)),1));
        figure(1), clf
        %add in the (1:end-73) to ensure only 200ms of match period shows
        erpD = plot(time_x(1:end-73),correct(1:end-73),time_x(1:end-73),incorrect(1:end-73),':', 'LineWidth', 2);
        set(gca,'box','off','Xlim',[time_x(1);time_x(end-73)]);
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
        fileName = sprintf('m%s_day_%d_of_%d_erp',monk,find(strcmp(alldays,dday{:})),size(fieldnames(data.(monkey{:})(1)),1));
        export_fig(fileName,'-pdf','-transparent'); %save transparent pdf in pwd
    end
end

%% Overall ERPS across days across chans by Monkey by behResp

clf
for monkey=monkeys
    alldays = fieldnames(data.(monkey{:})(1))'; %extract all days
    i=0; %reset day counter
    for dday=alldays
        i=i+1; %iterate day count
        allchans = size(data.(monkey{:})(1).(dday{:}).erp,1);
        %make day-level erp across chans
        daysERP(1,:,i) = mean(data.(monkey{:})(1).(dday{:}).erp,1); %correct
        daysERP(2,:,i) = mean(data.(monkey{:})(2).(dday{:}).erp,1); %inc
        trialCounts(1,i) = size(data.(monkey{:})(1).(dday{:}).lfp,3);
        trialCounts(2,i) = size(data.(monkey{:})(2).(dday{:}).lfp,3);
    end
    alldaysERP = mean(daysERP,3); %overall mean across days
    correct = alldaysERP(1,:)-alldaysERP(1,1); %start at 0
    incorrect = alldaysERP(2,:)-alldaysERP(2,1); %start at 0
    alltrialCounts = sum(trialCounts,2); %overall counts for trials
    if find(strcmp(monkeys,monkey{:})) == 1 monk="A";
    else monk="B"; end
    erptitle = sprintf('%d %s & %d %s Trial-averaged Monkey %s All Days All Areas & Chans',...
            alltrialCounts(1),"Correct",alltrialCounts(2),"Incorrect",monk);
%     erptitle = sprintf('Monkey %s All Days All Areas & Chans',monk);
    figure(1), clf
    %add in the (1:end-73) to ensure only 200ms of match period shows
    erpD = plot(time_x(1:end-73),correct(1:end-73),time_x(1:end-73),incorrect(1:end-73),':', 'LineWidth', 2);
    set(gca,'box','off','Xlim',[time_x(1);time_x(end-73)]);
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
    fileName = sprintf('m%s_all_%d_days_erp',monk,size(fieldnames(data.(monkey{:})(1)),1));
    export_fig(fileName,'-pdf','-transparent'); %save transparent pdf in pwd
end

%Next up: scan to see which days to take out of the Monkey B all-days ERP
%then create a Monkey B ERP without those days, likely 4-5 near the end