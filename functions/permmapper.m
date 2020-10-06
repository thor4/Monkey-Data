%% function for creating tf condition difference maps to be used in 
%  permutation testing. 
%  assumes correct trial count > incorrect
% 
% *input*
% mData is data structure containing power (microvolts squared) for 1 monkey,
% all areas, responses and trials across all days/sessions/trials
% agnostic to which monkey.
% n_permutes is the number of permutations to run, ie: 1000
% num_frex is total number of frequencies from analysis, ie: 35 
% times2save is down-sampled time points vector, ie: <10,20,30,...>
% 
% *output*
% metaperm is a struct with the following fields: 
% diffmap is a 4D tensor: [allchans perm freqidx timeidx] avg diffmap over
% trials
% area is a 1xN array where N is total number of chans across days (allchans),
% lists which area each chan is in
% day is a 1xN array where N is total number of chans across days (allchans),
% lists which day each chan was recorded during
% power measured in  in LFP data

function metaperm = permmapper(mData,n_permutes,num_frex,times2save)
alldays = fieldnames( mData )';
metaperm = []; %init meta permutation struct
diffpermmaps = []; %init allchan across days diff perm map array
for dayN=alldays
    allchans = size(mData(1).(dayN{:}).power,1); %total # of chans
    areas = string(mData(1).(dayN{:}).areas); %all areas
    % init H0 perm map for all chans [area x permutation x freqidx x timeidx]
    chan_permmaps = zeros(length(areas),n_permutes,num_frex,length(times2save));
    for chanN = 1:allchans
        % total number of incorrect trials for chan, power: chan x freqidx x time x trials
        nitrials = size( mData(2).(dayN{:}).power,4 );
        %initialize null hypothesis permutation-level maps
        permmaps = zeros(n_permutes,num_frex,length(times2save));
        for permi = 1:n_permutes
            % randomly sample from condition 1 trials (decimation) to match
            % condition 2
            randcoridx = randperm( size( mData(1).(dayN{:}).power,4 ),nitrials );
            temp_cor = squeeze( mData(1).(dayN{:}).power(chanN,:,:,randcoridx) );
            % concatenate conditions: trials1:nitrials are from correct, 
            % trials nitrials+1:end are from incorrect
            tf3d = cat(3,temp_cor,squeeze(mData(2).(dayN{:}).power(chanN,:,:,:)));
            % randomize trials, which also randomly assigns trials to 
            % conditions, correct vs incorrect
            randtf3didx = randperm(size(tf3d,3));
            temp_tf3d = tf3d(:,:,randtf3didx); % [freqidx x time x trials]
            % compute the "difference" map under the null hypothesis: 
            % [freqidx x time]
            permmaps(permi,:,:) = squeeze( mean(temp_tf3d(:,:,1:nitrials),3) - mean(temp_tf3d(:,:,nitrials+1:end),3) );
        end
        chan_permmaps(chanN,:,:,:) = permmaps;
%         find(chan_permmaps(chanN,:,:,:)~=0); %testing
%         find(permmaps~=0); %testing
%             size(mAchan_diffpermmaps) %testing
%             find(mAchan_diffpermmaps(chanN-8,:,:,:)~=0); %testing
    end
    dayy=find(ismember(alldays,dayN{:}));
    if dayy==1 
        diffpermmaps = chan_permmaps;
        metaperm.chan = [1:allchans]; %1 x allchans vector
        metaperm.area = areas;
        metaperm.day = repmat(dayy,allchans,1)';
    else
        diffpermmaps = cat(1,diffpermmaps,chan_permmaps);
        metaperm.chan = cat(2,metaperm.chan,[1:allchans]);
        metaperm.area = cat(2,metaperm.area,areas);
        metaperm.day = cat(2,metaperm.day,repmat(dayy,allchans,1)');
    end
end
metaperm.diffmap = diffpermmaps;
end
