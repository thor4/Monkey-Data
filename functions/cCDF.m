%function for calculating the complementary cumulative distribution
%function (cCDF)
%by Bryan Conklin 8/19/18. ~Thanks to "Fundamentals of Brain Network
%Analysis"
%useful for visually checking whether a degree distribution follows a power
%law
function [deg,Pdk] = cCDF(id)
    % deg is the span of ranked degrees which match their probability of
    % being greater than some random degree k (Pdk)
    % Pdk is the P(degree > k), id is the in-degree of a directed adjacency
    % matrix
    n = length(data);
    maxmin_range = max(data)-min(data);
    % note: the function iqr (inter-quartile range) is in the stats toolbox. 
    % In case the toolbox is absent, write a similar function by sorting 
    % the values, finding the values that are 25% and 75% of the sorted 
    % distribution, and then subtracting the 25% number from the 75% number.
    nBins = ceil(maxmin_range/(2.0*iqr(data)*n^(-1/3))); %Freedman-Diaconis 
end