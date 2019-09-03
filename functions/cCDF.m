%function for calculating the complementary cumulative distribution
%function (cCDF)
%by Bryan Conklin 8/19/18. ~Thanks to "Fundamentals of Brain Network
%Analysis"
%useful for visually checking whether a degree distribution follows a power
%law
function [data_sort,Pdk,idx] = cCDF(data)
    % data_sort is the ranked degrees which match their probability of
    % being greater than some random degree k (Pdk)
    % Pdk is the P(degree > k), idx is the index of the sorted id elements
    % data is a vector containing the degree counts of network nodes
    [data_sort,idx] = sort(data); %sort data in ascending order, save indices
    rank = (1:length(data)); %create rank vector for each node
    Pdk = zeros(1,length(data));  %init Pdk vector
    ii = [0, diff(data_sort(:)')==0,0]; %find where there are repetitions
    beg_idx = strfind(ii,[0 1]); %starting idx of repetitions
    end_idx = strfind(ii,[1 0]); %ending idx of repetitions
    filled_idx = 0; %init
    for i=length(end_idx):-1:1 %use highest-rank for repeated degrees
        Pdk(beg_idx(i):end_idx(i)) = 1-(rank(end_idx(i))/length(data));
        filled_idx = [filled_idx beg_idx(i):end_idx(i)]; %keep track of all 
        %indices with repetitions
    end
    filled_idx(filled_idx==0) = []; %take init 0 out
    for i=length(data):-1:1 %fill in rest of probability for non-rep's
        if ~(ismember(i,filled_idx)) %ensure didn't already calculate
            Pdk(i) = 1-(rank(i)/length(data));
        end
    end
end