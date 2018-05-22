%function to plot ERP for a recording session
%monkey should be the nested data struct from 
%mGoodStableRule1PingRej-split_by_Day_BehResp_and_Area.mat
%monk should be an integer: 1 for Monkey 1, 2 for Monkey 2
%d should be an integer (day) 1-23 for Monkey 1 or 1-24 for Monkey 2
function session_erp(monkey,monk,d)
    %build out matrices for each frontal and parietal region for the day
    %initialize variables
    resp = fieldnames(monkey(monk).day(d));
    frontal = {'a6DR','a8B','a8AD','a9L','adPFC','avPFC'};
    parietal = {'aLIP','aMIP','aPEC','aPE','aPG'};
    correctF=[]; correctP=[]; incorrectF=[]; incorrectP=[];
    for i=1:numel(resp)
        area = fieldnames(monkey(monk).day(d).(resp{i}));
        for j=1:numel(area)
            if (ismember(area(j),frontal)) && (resp{i} == "CorrectRule1")
                correctF = vertcat(correctF,monkey(monk).day(d).(resp{i}).(area{j}));
            elseif (ismember(area(j),frontal)) && (resp{i} == "IncorrectRule1")
                incorrectF = vertcat(incorrectF,monkey(monk).day(d).(resp{i}).(area{j}));
            elseif (ismember(area(j),parietal)) && (resp{i} == "CorrectRule1")
                correctP = vertcat(correctP,monkey(monk).day(d).(resp{i}).(area{j}));
            elseif (ismember(area(j),parietal)) && (resp{i} == "IncorrectRule1")
                incorrectP = vertcat(incorrectP,monkey(monk).day(d).(resp{i}).(area{j}));
            end
        end
    end
    Triggers = [0 500]; 
    %compute the average of each region across all instances per time point
    %result should be 1 x 1811 (sample points) matrix
    correctFavg = mean(correctF,1); incorrectFavg = mean(incorrectF,1);
    correctPavg = mean(correctP,1); incorrectPavg = mean(incorrectP,1);
    %shift all averages to start at zero
    correctFavgShift = correctFavg + (-correctFavg(1));
    incorrectFavgShift = incorrectFavg + (-incorrectFavg(1));
    correctPavgShift = correctPavg - (correctPavg(1));
    incorrectPavgShift = incorrectPavg + (-incorrectPavg(1));
    %figure
    time  = -500:size(correctFavg,2)-501; % time, from -500ms baseline
    figure
    subplot(2,1,1)
    traces1 = plot(time,correctFavgShift,time,incorrectFavgShift,':', 'LineWidth', 2);
    title(['Trial-averaged Monkey ' num2str(monk) ' Frontal Region Time-series']); xlabel('Time (ms)'); ylabel('Voltage (µV)'); 
    set(gca,'box','off','Xlim',[time(1);time(end)]);
    y1 = get(gca,'ylim'); hold on
    triggers1 = plot([Triggers(1) Triggers(1)],y1,'--',[Triggers(2) Triggers(2)],y1,'--'); hold off 
    legend(traces1,'Correct Trials','Incorrect Trials');
    subplot(2,1,2)
    traces2 = plot(time,correctPavgShift,time,incorrectPavgShift,':', 'LineWidth', 2);
    title(['Trial-averaged Monkey ' num2str(monk) ' Parietal Region Time-series']); xlabel('Time (ms)'); ylabel('Voltage (µV)');
    set(gca,'box','off','Xlim',[time(1);time(end)]);
    y2 = get(gca,'ylim'); hold on
    triggers2 = plot([Triggers(1) Triggers(1)],y2,'--',[Triggers(2) Triggers(2)],y2,'--'); hold off
    legend(traces2,'Correct Trials','Incorrect Trials');
end