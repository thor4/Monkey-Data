%function for correlation & covariance matrix
%data should be a matrix of instances(channels) x samples(time)
%it's measured in terms of the sample metric (microvolts in LFP data)
function [monkeyCor, monkeyCov] = monkeyCov(data)
    %find entire instances where a sample falls outside the ping range
    monkeyData = bsxfun(@minus,data,mean(data,2));
    %unscaled covariance matrix
    monkeyCov = monkeyData*monkeyData'/(size(data,2)-1);
    imagesc(monkeyCov)
    colorbar
    %the diagonal of unscaled covariance matrix are the variances of each 
    %electrode. now, compute the "variance matrix"
    stdMat = sqrt( diag(monkeyCov)*diag(monkeyCov)' );
    %scale the (previously unscaled) covariance
    monkeyCor = monkeyCov ./ stdMat;
end