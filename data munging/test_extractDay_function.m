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

% begin to build loop for ERP analysis in time_domain analysis
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
tic
for lp=alldays %cycle through all days
    dayy = append('d',lp{:});
    [data.(dayy).lfp,data.(dayy).areas] = extractDay(path,monkey,lp{:},1,2,0,1,'match');
end
toc
%94 seconds  for baseline on home pc
%29.27 seconds for sample on labpc
%151.08 seconds for delay on home pc
%27.30 seconds for match on home pc (2nd time)

clear data %reset each time


%unit test for craw, compare to data struct
day = alldays{21}; %assign day
dayy = append('d',day); %setup day for test indexing
%switch over to craw function to load j and trial_info_path, then load 
%trial info for day (line 92)

%find trials where parameters are met
for k=3:1000 %don't account for stability
    if (trial_info.good_trials(k) == 1) && ...%artifacts/none
            (trial_info.BehResp(k) == 0) && ... %correct/incorrect
            (trial_info.rule(k) == 1) %identity/location
        break
    end
end

trial3 = data.(dayy).lfp(:,:,2); %extract 1st trial to compare with raw, k is in trial name

%load raw lfp
lfp_path = strcat(path,'%s\\%s\\%s\\%s%s%s.%04d.mat'); %build lfp raw data path
j=2; %betty and session 2 clark
j=3; %session 3 clark
trial_lfp = sprintf(lfp_path,monkey,day,days{j},monkey,day,days{j}{1}(8:9),k);
load(trial_lfp,'lfp_data');

floor(trial_info.CueOnset(k))-504 %beginning of baseline period to compare
floor(trial_info.CueOnset(k)) %beginning of sample period to compare
floor(trial_info.CueOnset(k))+504 %end of sample period to compare
floor(trial_info.CueOffset(k)) %beginning of delay period to compare
floor(trial_info.CueOffset(k))+810 %end of delay period to compare
floor(trial_info.MatchOnset(k)) %beginning of match period to compare
floor(trial_info.MatchOnset(k))+273 %end of match period to compare




% baseline tested fine for betty 090615 correct
% sample tested fine for betty 090702 incorrect
% delay tested fine for betty 090709 correct
% match tested fine for betty 090921 correct
% baseline tested fine for clark 060406 correct
% sample tested fine for clark 060503 incorrect
% delay tested fine for clark 060825 incorrect
% match tested fine for clark 061214 incorrect

% trial_info.Trial_Length is *not correct* for betty day 1