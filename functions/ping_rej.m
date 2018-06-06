%function for ping rejection
%data should be a matrix of instances(trials) x samples(time)
%upper and lower are the bounds within which data will be kept
%they are measured in terms of the sample metric (microvolts in LFP data)
function [noping,pings,throwAwayIdx] = ping_rej(data,lower,upper)
    %find entire instances where a sample falls outside the ping range
    throwAwayIdx = any(data > upper | data < lower, 2);
    %total instances, inclusive of any pings
    total = size(data,1);
    %eliminate entire instances where a sample falls outside acceptable
    %range
    data(throwAwayIdx,:) = [];
    %report total number of pings found
    pings = total - size(data,1);
    noping = data;
end