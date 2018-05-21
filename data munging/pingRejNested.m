%ping rejection for Monkey 1 Correct Trials [-750,750]
for i=1:numel(monkey(1).day)
    area = fieldnames(monkey(1).day(i).CorrectRule1);
    for j=1:numel(area)
        [noping,pings] = ping_rej(monkey(1).day(i).CorrectRule1.(area{j}),-750,750);
        if (pings > 0)
            %pull out all the pings and replace with artifact-free data
            monkeyNew(1).day(i).CorrectRule1.(area{j}) = noping;
            %keep track of total number of pings per area per day
            counts.day(i).(area{j}) = pings;
        end
    end
end

%ping rejection for Monkey 1 Incorrect Trials [-750,750]
for i=1:numel(monkey(1).day)
    area = fieldnames(monkey(1).day(i).IncorrectRule1);
    for j=1:numel(area)
        [noping,pings] = ping_rej(monkey(1).day(i).IncorrectRule1.(area{j}),-750,750);
        if (pings > 0)
            %pull out all the pings and replace with artifact-free data
            monkeyNew(1).day(i).IncorrectRule1.(area{j}) = noping;
            %keep track of total number of pings per area per day
            counts.day(i).(area{j}) = pings;
        end
    end
end
