%function for Shannon information
%by Bryan Conklin 8/17/18. ~Thanks to Mike X Cohen's "Analyzing Neural
%Time Series Data"
%data should be a vector of a single instance(channel) x samples(time)
%nbins should be optimal number of bins for histogram creation computed 
%from numBins function
%doesn't matter what data is measured in since Shannon's information only
%concerns itself with probabilities and distributions, not the raw data
function entro = shanInfo(data,nbins)
    %partitions data into bins, and returns the count in each bin
    %(countdata), as well as the bin edges (edgesdata).
    [countdata,edgesdata] = histcounts(data,nbins); 
    % convert counts to probability values
    countdata = countdata./sum(countdata);
    %compute entropy
    entro = -sum(countdata.*log2(countdata+eps));
end

%Why eps? It is possible that there are zero probability values for some 
%bins. This is problematic because the logarithm of zero is undefined. 
%Thus, in Matlab code, a very small number is added inside the logarithm 
%term to prevent taking the logarithm of zero. The built-in variable eps is
%used as this small number as it is the smallest numerical resolution in 
%MATLAB. Thus, instead of writing log2(p), log2(p+eps) is written. This 
%adds such a tiny number as to be insignificant to non-zero probabilities, 
%and will prevent taking the logarithm of zero. It will generate very large
%negative values when the probability is zero (e.g., –52), but this value 
%is then multiplied by zero [because of p*log2(p+eps)] and consequently 
%becomes irrelevant. Another option would be to eliminate all bins with 
%zero probabilities, but this may lead to complications in later parts of 
%the code because the number of bins could change across variables.