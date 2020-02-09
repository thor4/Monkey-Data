% monkey: [ "betty", "clark" ]
% day:(betty) [   '090615', '090616', '090617', '090618', '090622', '090625', 
%                 '090626', '090629', '090701', '090702', '090706', '090708', 
%                 '090709', '090901', '090903', '090916', '090917', '090921', 
%                 '090923', '090924', '090928', '090929', '090930', '091001']
%     (clark) [   '060328', '060406', '060411', '060414', '060426', '060427', 
%                 '060428', '060502', '060503', '060509', '060511', '060531', 
%                 '060601', '060602', '060824', '060825', '060831', '060907', 
%                 '061212', '061213', '061214', '061215', '061221'          ]
% good: [ 0(artifacts), 1(no artifacts) ]
% stable: [ 0(transition), 1(stable performance), 2(both) ]
% behResp: [0(incorrect), 1(correct) ]
% rule: [ 1(identity), 2(location) ]
% epoch: [ 'base', 'sample', 'delay', 'match', 'all' ]

% next up, test the clark sample period to ensure it's pulling correctly
% then run the clark sample period to get counts
monk = 1; %1 = clark, 2 = betty
if monk==1
    monkey='clark';
    alldays = days_clark;
    days = { days_clark, "session02", "session03" };
else
    monkey='betty';
    alldays = days_betty;
    days = { days_betty, "session01" };
end
% 
i=1; %init counter
for lp=alldays %cycle through all days
%     doSomeOperation( mystruct( lp{:} ) );
    dayy = append('d',lp{:});
%     [lengths.(dayy),idx.(dayy)] = craw(path,monkey,lp{:},1,2,1,1,'base');
    [lengths.(dayy),numtrials] = craw(path,monkey,lp{:},1,2,0,1,'base');
    lengths.(dayy)(lengths.(dayy)==0)=[]; %get rid of placeholders
    numtrials(numtrials==0)=[]; %get rid of placeholders
    idx(i) = numtrials; %save number of trials per day in array
    mins(i)=min(lengths.(dayy)); %find the shortest length
    i=i+1;
end

idx=idx'; mins=mins'; %easier to copy-paste
clear idx mins lengths %reset each time

fn=fieldnames(lengths);

[lengths,idx] = craw(path,monkey,'090615',1,2,1,1,'base');


%unit test for craw counter, compare to 'lengths' vector
day = alldays{1}; %assign day
%switch over to craw function to load trial info for day (line 98)
length(1:floor(trial_info.CueOnset(k))-1) %epoch length for baseline
length(floor(trial_info.CueOnset(k)):floor(trial_info.CueOffset(k))-1) %epoch length for sample
length(floor(trial_info.CueOffset(k)):floor(trial_info.MatchOnset(k))-1) %epoch  length for delay
%for match
lfp_path = strcat(path,'%s\\%s\\%s\\%s%s%s.%04d.mat'); %build lfp raw data path
j=2; %betty and session 2 clark
j=3; %session 3 clark
trial_lfp = sprintf(lfp_path,monkey,day,days{j},monkey,day,days{j}{1}(8:9),k);
load(trial_lfp,'lfp_data');
length(floor(trial_info.MatchOnset(k)):length(lfp_data)) %epoch length for match
% length(floor(trial_info.MatchOnset(k)):floor(trial_info.TrialLength(k))) %epoch length for match
for k=1:1000 %don't account for stability
    if (trial_info.good_trials(k) == 1) && ...%artifacts/none
            (trial_info.BehResp(k) == 0) && ... %correct/incorrect
            (trial_info.rule(k) == 1) %identity/location
        break
    end
end


for k=1:trial_info.numTrials %don't account for stability
    if (trial_info.good_trials(k) == good) && ...%artifacts/none
            (trial_info.BehResp(k) == behResp) && ... %correct/incorrect
            (trial_info.rule(k) == rule) %identify/location
        break
    end
end
% 
% baseline tested fine for betty 090615 correct
% sample tested fine for betty 090616 incorrect
% delay tested fine for betty 090702 incorrect
% match tested fine for betty 090903 incorrect
% match tested fine for betty 090615 incorrect after craw floor modification
% trial_info.Trial_Length is *not correct* for betty day 1
% baseline tested fine for clark 060509 incorrect
% match tested fine for clark 060328 incorrect


%% test string validation
p = inputParser;
argName = 'monkey';
monkeys = { 'toejam','earl' };
% validationFcn = @(x) any(validatestring(x,monkeys));
validationFcn = @(x) any(strcmp(x,monkeys));
addRequired(p,argName,validationFcn);

parse(p,'toejam')


p = inputParser;
argName = 'monkey';
monkeys = [ "toejam","earl" ];
validationFcn = @(x) validateStringParameter(x,monkeys,mfilename,argName);
addRequired(p,argName,validationFcn);

function validateStringParameter(varargin)
    validatestring(varargin{:});
end

@(inVal)any(strcmp(inVal,possVal)) 

