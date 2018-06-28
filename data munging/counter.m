%day
d = 'd%i';

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

monkeyNew(1).day(i).correct.(chan{j});