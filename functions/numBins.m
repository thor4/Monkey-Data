%function for optimal number of bins for Shannon information
%by Bryan Conklin 8/19/18. ~Thanks to Mike X Cohen's "Analyzing Neural
%Time Series Data"
%data should be a vector of a single instance(channel) x samples(time)
%doesn't matter what data is measured in since Shannon's information only
%concerns itself with probabilities and distributions, not the raw data
function nBins = numBins(data)
    % optimal number of bins for histogram based on Freedman-Diaconis rule
    % (Freedman and Diaconis 1981). It states the optimal number of
    % histogram bins is related to the interquarterile range and the number
    % of data points. Does not assume data is normally distributed, unlike
    % others like Scott's rule (Scott 2010). This means it can be applied
    % to any data distribution, which is usefule for EEG/neural
    % time-frequency data which can have a power distrib. or circular
    % distrib. When data *is* normally distrib., FD and Scott provide same
    % or similar number of bins.
    n = length(data);
    maxmin_range = max(data)-min(data);
    % note: the function iqr (inter-quartile range) is in the stats toolbox. 
    % In case the toolbox is absent, write a similar function by sorting 
    % the values, finding the values that are 25% and 75% of the sorted 
    % distribution, and then subtracting the 25% number from the 75% number.
    nBins = ceil(maxmin_range/(2.0*iqr(data)*n^(-1/3))); %Freedman-Diaconis 
end

% Another point on binning:
% The procedures for selecting an appropriate number of bins outlined above are for one
% variable. When you are computing mutual information (or computing entropy for multiple
% variables), choosing an appropriate number of bins becomes more difficult because the optimal
% number of bins might be different for each variable. Because the number of bins will
% influence the estimate of entropy, the same number of bins should be used for computing
% entropy in all conditions, electrodes, and subjects. This is demonstrated in section 29.9.
% Thus, a recommended strategy is to compute the optimal number of bins for each variable,
% take the ceiling (i.e., round up) of the average of the optimal number of bins, and then apply
% this to all variables. 