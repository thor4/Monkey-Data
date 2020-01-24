load('mGoodStableRule1-split_by_Day_BehResp_and_Chan.mat')
%day
d = 'd%i';
%correct
for i=1:numel(monkey(1).day)
    chan = fieldnames(monkey(1).day(i).correct);
    day = sprintf(d, i);
    for j=1:numel(chan)
        [~,pings,idx] = ping_rej(monkey(1).day(i).correct.(chan{j}),-750,750);
        if (pings > 0)
            %keep track of total number of pings per area per day
            counts.day(i).correct.(day).(chan{j}) = pings;
        end
    end
end

%incorrect
for i=1:numel(monkey(1).day)
    chan = fieldnames(monkey(1).day(i).incorrect);
    day = sprintf(d, i);
    for j=1:numel(chan)
        [~,pings,idx] = ping_rej(monkey(1).day(i).incorrect.(chan{j}),-750,750);
        if (pings > 0)
            %keep track of total number of pings per area per day
            counts.day(i).incorrect.(day).(chan{j}) = pings;
        end
    end
end
