%% function for creating tf condition difference maps to be used in 
%  permutation testing. 
%  assumes correct trial count > incorrect
% 
% *input*
% monkey is data structure containing power for both monkeys, all areas, 
%   responses and trials across all days/sessions/trials
% monkeyN should be which monkey: '1' for A or '2' for B.
% n_permutes is the number of permutations to run, ie: 1000
% num_frex is total number of frequencies from analysis, ie: 35 
% times2save is down-sampled time points vector, ie: <10,20,30,...>
% 
% *output*
% area_permmaps is a 4-d area (length(areas)) x permmap (n_permutes) x 
% freq (num_frex) x samples (length(times2save))
% areas is which monkey's areas correspond to 1st dim of area_permmaps
% power measured in microvolts in LFP data
function [area_permmaps, areas] = permmapper(monkey,monkeyN,n_permutes,num_frex,times2save)
    % define recorded areas from each monkey
    if monkeyN==1
        areas = {'a8B', 'a9L', 'adPFC', 'avPFC', 'aLIP', 'aMIP', 'aPEC', 'aPG'}; %monkey1
    else
        areas = {'a6DR', 'a8AD', 'a8B', 'adPFC', 'aLIP', 'aPE', 'aPEC', 'aPG'}; %monkey2
    end
    % initialize null hypothesis permutation maps [permutation x freqidx x timeidx]
    permmaps = zeros(n_permutes,num_frex,length(times2save));
    % and for all areas [area x permutation x freqidx x timeidx]
    area_permmaps = zeros(length(areas),size(permmaps,1),size(permmaps,2),size(permmaps,3));
    % total number of incorrect trials for area
    for areaN = 1:numel(areas)
        ntrials = size( monkey(monkeyN).incorrect.(areas{areaN}),3 );
        for permi = 1:n_permutes
            % randomly sample from condition1 trials (decimation) to match
            % condition2
            randcoridx = randperm( size( monkey(monkeyN).correct.(areas{areaN}),3 ),ntrials );
            temp_cor = monkey(monkeyN).correct.(areas{areaN})(:,:,randcoridx);
            % concatenate conditions: trials1:ntrials are from correct, 
            % trials ntrials+1:end are from incorrect
            tf3d = cat(3,temp_cor,monkey(monkeyN).incorrect.(areas{areaN}));

            % randomize trials, which also randomly assigns trials to 
            % conditions, correct vs incorrect
            randtf3didx = randperm(size(tf3d,3));
            temp_tf3d = tf3d(:,:,randtf3didx);

            % compute the "difference" map under the null hypothesis
            permmaps(permi,:,:) = squeeze( mean(temp_tf3d(:,:,1:ntrials),3) - mean(temp_tf3d(:,:,ntrials+1:end),3) );
        end
        area_permmaps(areaN,:,:,:) = permmaps;
        % re-initialize null hypothesis maps
        permmaps = zeros(n_permutes,num_frex,length(times2save));
    end
end