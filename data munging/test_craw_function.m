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

% next up, test sample, delay and match of betty to ensure correctness
% then test clark to ensure correct epochs get pulled for her
% then begin to build loop for ERP analysis in time_domain analysis
monk = 2; %1 = clark, 2 = betty
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
    [data.(dayy).lfp,data.(dayy).areas] = craw(path,monkey,lp{:},1,2,1,1,'base');
end
toc
%94 seconds  for baseline

clear data %reset each time


%unit test for craw, compare to data struct
day = alldays{1}; %assign day
%switch over to craw function to load trial info for day (line 92)
[lfp,areas] = craw(path,monkey,day,1,2,1,1,'base'); %one day
[data.(dayy).lfp,data.(dayy).areas] = craw(path,monkey,day,1,2,1,1,'base'); %one day
trial3 = lfp(:,:,3); %extract trial to compare with raw, k will differ

%find trials where parameters are met
for k=7:1000 %don't account for stability
    if (trial_info.good_trials(k) == 1) && ...%artifacts/none
            (trial_info.BehResp(k) == 1) && ... %correct/incorrect
            (trial_info.rule(k) == 1) %identity/location
        break
    end
end

%load raw lfp
lfp_path = strcat(path,'%s\\%s\\%s\\%s%s%s.%04d.mat'); %build lfp raw data path
j=2; %betty and session 2 clark
j=3; %session 3 clark
trial_lfp = sprintf(lfp_path,monkey,day,days{j},monkey,day,days{j}{1}(8:9),k);
load(trial_lfp,'lfp_data');

floor(trial_info.CueOnset(k))-504 %beginning of sample period to compare




% baseline tested fine for betty 090615 correct

% sample tested fine for betty 090616 incorrect
% delay tested fine for betty 090702 incorrect
% match tested fine for betty 090903 incorrect
% match tested fine for betty 090615 incorrect after craw floor modification
% trial_info.Trial_Length is *not correct* for betty day 1
% baseline tested fine for clark 060509 incorrect
% sample tested fine for clark 060426 incorrect
% delay tested fine for clark 060824 correct (j=3)
% match tested fine for clark 060328 incorrect

