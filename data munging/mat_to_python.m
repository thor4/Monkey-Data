%transfer large .mat file over to smaller .mat's so that python scipy can
%read them in. monkey 2 needed to be split between days 1-12 and days
%13-24 then by correct and incorrect. monkey 1 split on correct and 
%incorrect.

%Monkey 1 
%correct
for i=1:numel(monkey(1).day)
    chan = fieldnames(monkey(1).day(i).correct);
    day = sprintf(d, i);
    for j=1:numel(chan)
        monkey1.day(i).correct.(chan{j}) = monkey(1).day(i).correct.(chan{j});
    end
end

%incorrect
for i=1:numel(monkey(1).day)
    chan = fieldnames(monkey(1).day(i).incorrect);
    day = sprintf(d, i);
    for j=1:numel(chan)
        monkey1.day(i).incorrect.(chan{j}) = monkey(1).day(i).incorrect.(chan{j});
    end
end


%Monkey 2
%correct
for i=1:numel(monkey(2).day)/2
    chan = fieldnames(monkey(2).day(i).correct);
    day = sprintf(d, i);
    for j=1:numel(chan)
        monkey2C1st.day(i).correct.(chan{j}) = monkey(2).day(i).correct.(chan{j});
    end
end

for i=numel(monkey(2).day)/2+1:numel(monkey(2).day)
    chan = fieldnames(monkey(2).day(i).correct);
    day = sprintf(d, i);
    for j=1:numel(chan)
        monkey2C2nd.day(i).correct.(chan{j}) = monkey(2).day(i).correct.(chan{j});
    end
end

%incorrect
for i=1:numel(monkey(2).day)/2
    chan = fieldnames(monkey(2).day(i).incorrect);
    day = sprintf(d, i);
    for j=1:numel(chan)
        monkey2I1st.day(i).incorrect.(chan{j}) = monkey(2).day(i).incorrect.(chan{j});
    end
end

for i=numel(monkey(2).day)/2+1:numel(monkey(2).day)
    chan = fieldnames(monkey(2).day(i).incorrect);
    day = sprintf(d, i);
    for j=1:numel(chan)
        monkey2I2nd.day(i).incorrect.(chan{j}) = monkey(2).day(i).incorrect.(chan{j});
    end
end