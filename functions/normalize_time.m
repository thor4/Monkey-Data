%function to normalize temporal data to be used in time series analysis
%data & baseline should be a matrix of instances(trials) x samples(time)
function norm = normalize_time(baseline, data)
    %this will average across all baseline samples for each instance
    norm = data - mean(baseline,2);
end