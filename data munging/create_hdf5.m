%create an hdf5 file from monkey data
%first create the group structure in HDFview
%next loop through and add all areas for each day, behresp and chan

d = 'd%i';
%correct
for i=1:numel(monkey(2).day)
    chan = fieldnames(monkey(2).day(i).correct);
    day = sprintf(d, i);
    for j=1:numel(chan)
        mydata = monkey(2).day(i).correct.(chan{j});
        path_unfill = '/monkey2/correct/%s/%s/';
        path_fill = sprintf(path_unfill,day,chan{j});
        h5create('mGoodStableRule1PingRej.h5',path_fill,size(mydata))
        h5write('mGoodStableRule1PingRej.h5',path_fill,mydata)
    end
end

%incorrect
for i=1:numel(monkey(2).day)
    chan = fieldnames(monkey(2).day(i).incorrect);
    day = sprintf(d, i);
    for j=1:numel(chan)
        mydata = monkey(2).day(i).incorrect.(chan{j});
        path_unfill = '/monkey2/incorrect/%s/%s/';
        path_fill = sprintf(path_unfill,day,chan{j});
        h5create('mGoodStableRule1PingRej.h5',path_fill,size(mydata))
        h5write('mGoodStableRule1PingRej.h5',path_fill,mydata)
    end
end