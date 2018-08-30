function [mi entropy fd_bins p_n p_z] = mutual_info(x,y,fd_bins,permtest)
% MUTUAL INFORMATION   Compute mutual information between two vectors
%  
%   Inputs:
%       x,y     :  data matrices of equal size
%
%   Optional inputs:
%       bins    :  number of bins to use for distribution discretization.
%                  be sure to run numBins function to find overall avg of 
%                  optimal number of variables first
%       permtest : perform permutation test and return mi in standard-Z 
%                  values
%
%   Outputs:
%       mi      :  mutual information in bits
%       entropy :  entropy of x, y, and joint
%       nbins   :  number of bins used for discretization
%                  (based on Freedman-Diaconis rule)
%       p_n     :  (only with permtest) p-value associated with # of 
%                  supra-threshold tests. should be similar to p_z
%       p_z     :  (only with permtest) p-value associated with Z-value's
%                  position on a Gaussian probability density (normcdf)
%
%  Bryan Conklin 08/25/18
%  Adapted from Mike X Cohen's text Neural Time Series Analysis

if nargin<2, error('Specify two inputs.'); end
if length(x)~=length(y), error('X and Y must have equal length'); end

%% determine the optimal number of bins for each variable

% vectorize in the case of matrices
% turns a 1 x time x trial matrix into a vector to compute MI over time and
% trials
x=x(:); y=y(:);

if nargin<3 || isempty(fd_bins)
    n            = length(x);
    maxmin_range = max(x)-min(x);
    fd_bins1     = ceil(maxmin_range/(2.0*iqr(x)*n^(-1/3))); % Freedman-Diaconis
    
    n            = length(y);
    maxmin_range = max(y)-min(y);
    fd_bins2     = ceil(maxmin_range/(2.0*iqr(y)*n^(-1/3)));
    
    % and use the average...
    fd_bins = ceil((fd_bins1+fd_bins2)/2);
end

%% bin data

edges = linspace(min(x),max(x),fd_bins+1);
[nPerBin1,edges1,bins1] = histcounts(x,edges); 

edges = linspace(min(y),max(y),fd_bins+1);
[nPerBin2,edges2,bins2] = histcounts(y,edges); 

%% compute entropies

% recompute entropy with optimal bins for comparison
hdat1 = nPerBin1./sum(nPerBin1);
hdat2 = nPerBin2./sum(nPerBin2);

% convert histograms to probability values
for i=1:2
    eval([ 'entropy(' num2str(i) ') = -sum(hdat' num2str(i) '.*log2(hdat' num2str(i) '+eps));' ]);
end

%% compute joint probabilities

jointprobs = zeros(fd_bins);
for i1=1:fd_bins
    for i2=1:fd_bins
        jointprobs(i1,i2) = sum(bins1==i1 & bins2==i2);
    end
end
jointprobs=jointprobs./sum(jointprobs(:));

entropy(3) = -sum(jointprobs(:).*log2(jointprobs(:)+eps));

%% mutual information

mi = sum(entropy(1:2)) - entropy(3);

%% optional permutation testing

if nargin==4
    
    npermutes = 500; %define number of permutations
    n = length(bins2);
    
    perm_mi = zeros(1,npermutes); %initialize permutation MI matrix
    
    for permi=1:npermutes
        
        jointprobs = zeros(fd_bins); %initialize joint prob distribution
        
        % shuffle bin assignments for just one of the signals
        binbreak = randsample(round(n*.8),1,1)+round(n*.1); %this is like where to cut a deck of cards
        %max binbreak can be is ~90% of n
        %now permute the bin assignments according to the cut and whether
        %the iteration modulus 4 is 0-3.
        switch mod(permi,4)
            case 0, bins2 = [ bins2(binbreak:end);    bins2(1:binbreak-1) ];
            case 1, bins2 = [ bins2(end:-1:binbreak); bins2(1:binbreak-1) ];
            case 2, bins2 = [ bins2(binbreak:end);    bins2(binbreak-1:-1:1) ];
            case 3, bins2 = [ bins2(end:-1:binbreak); bins2(binbreak-1:-1:1) ];
        end
        
        %compute new joint probability distribution for original signal 1
        %bin and shuffled signal 2 bin
        for i1=1:fd_bins
            for i2=1:fd_bins
                jointprobs(i1,i2) = sum(bins1==i1 & bins2==i2);
            end
        end
        jointprobs=jointprobs./sum(jointprobs(:));
        
        perm_jentropy = -sum(jointprobs(:).*log2(jointprobs(:)+eps));
        
        % mutual information
        perm_mi(permi) = sum(entropy(1:2)) - perm_jentropy;
    end
    
    %calculate p-value associated with z-value
    %p_n is calculated by counting the number of null-hypothesis values
    %above the observed value and dividing by the number of permutations
    p_n = sum(perm_mi>mi)/npermutes;
       
    %calculate standard z-value version of mi score
    mi = (mi-mean(perm_mi))/std(perm_mi);
    
    %p_z is calculated by converting the observed value to standard
    %deviation units and evaluating the probability of that std dev under a
    %Gaussian distribution
    p_z = 1 - normcdf(mi);
    
end

%%
% simplified replacement for randsample
function y = randsample(x,n,junk)
y=randperm(x);
y=y(1:n);
