% New analysis 1/21/20: Trial-level power with stimulus info
%% Step 1: Pull out good, correct, stable, rule 1 + 2 trials from both monkeys %%
%initialize variables
% homepc:
path = 'H:\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\';
% labpc:
% path = 'C:\\Users\\bryan\\OneDrive\\Documents\\PhD @ FAU\\research\\High Frequency FP Activity in VWM\\';
monk = 1; %1 = clark, 2 = betty
days_betty = { '090615', '090616', '090617', '090618', '090622', '090625', '090626', '090629', '090701', '090702', '090706', '090708', '090709', '090901', '090903', '090916', '090917', '090921', '090923', '090924', '090928', '090929', '090930', '091001' };
days_clark = { '060328', '060406', '060411', '060414', '060426', '060427', '060428', '060502', '060503', '060509', '060511', '060531', '060601', '060602', '060824', '060825', '060831', '060907', '061212', '061213', '061214', '061215', '061221' };
if monk==1
    monkey='clark';
    alldays = days_clark;
    days = { days_clark, "session02", "session03" };
else
    monkey='betty';
    alldays = days_betty;
    days = { days_betty, "session01" };
end

%for testing:
% dayN = 13; %day number
% day = alldays{dayN}; %assign day
% dday = append('d',day); %setup day for test indexing

i=0; %init counter

%try saving each day's tri_info & rec_info: error in clark day where
%session 2 had 600 trials and session 3 had 579 due to a rule change in
%session 2 for trials 578-600. the trial_info.rule for session 3 goes up to
%600, while the rest end at 579. error is in extractDay:
%trial_info.(tdata{m})(tri_used)
%Index exceeds the number of array elements (579).
%figure out how to go through the sessions in clark and save trial data
%while preserving betty functionality

%change var name per monkey/rule combo (4x)
for lp=alldays %cycle through all days
    i = i+1;
    dayy = append('d',lp{:});
    [dataM1goodCorR1.(dayy).lfp,dataM1goodCorR1.(dayy).areas,...
        dataM1goodCorR1.(dayy).tri_info,dataM1goodCorR1.(dayy).rec_info] =...
        extractDay(path,monkey,lp{:},1,1,1,1,'all');
    %convert to µV (1V = 10^6µV = 1,000,000µV) 
    dataM1goodCorR1.(dayy).lfp = dataM1goodCorR1.(dayy).lfp .* 1e6;
end

%test
dataM1goodCorR1.(dayy).erp == mean(dataM1goodCorR1.(dayy).lfp,3) .* 1e6

save('time_domain-m1','dataM1goodCorR1','dataM1goodIncR1')
save('time_domain-m2','dataM2goodCorR1','dataM2goodIncR1','-v7.3') %v7.3 since filesize so big
%home pc: D:\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data