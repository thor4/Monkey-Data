%duplicate "monkey" var and rename to monkeyNew
%ping rejection for Monkey 1 Correct Trials [-750,750]
for i=1:numel(monkey(1).day)
    chan = fieldnames(monkey(1).day(i).correct);
    for j=1:numel(chan)
        [~,pings,idx] = ping_rej(monkeyNew(1).day(i).correct.(chan{j}),-750,750);
        if (pings > 0)
            %pull out all the pings from each channel
            for k=1:numel(chan)
                monkeyNew(1).day(i).correct.(chan{k})(idx,:) = [];
            end
            %keep track of total number of pings per area per day
            counts.day(i).correct.(chan{j}) = pings;
        end
    end
end

%ping rejection for Monkey 1 Incorrect Trials [-750,750]
for i=1:numel(monkey(1).day)
    chan = fieldnames(monkey(1).day(i).incorrect);
    for j=1:numel(chan)
        [~,pings,idx] = ping_rej(monkeyNew(1).day(i).incorrect.(chan{j}),-750,750);
        if (pings > 0)
            %pull out all the pings from each channel
            for k=1:numel(chan)
                monkeyNew(1).day(i).incorrect.(chan{k})(idx,:) = [];
            end
            %keep track of total number of pings per area per day
            counts.day(i).incorrect.(chan{j}) = pings;
        end
    end
end
