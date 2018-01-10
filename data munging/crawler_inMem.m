load('good_trials-betty.mat')
load('good_trials-clark.mat')

tic
zfc = 1;
zfi = 1;
zpc = 1;
zpi = 1;
fieldsBettyDays = fieldnames(bettyGoodTrials);
for i = 1:numel(fieldsBettyDays)
    fieldsBettyTrials = fieldnames(bettyGoodTrials.(fieldsBettyDays{i}));
    recordingRegion = bettyGoodTrials.(fieldsBettyDays{i}).recording_info.cortex;
    for j = 1:numel(fieldsBettyTrials)-2
        for k = 1:numel(recordingRegion)
            if (bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingRegion(k) == 'F') %only look at correct rule 1 frontal trials
                tempCorrectFrontalBetty(zfc,1:1001) = bettyGoodTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                tempCorrectFrontalBetty(zfc,1002:1811) = bettyGoodTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                zfc = zfc + 1;
            elseif (bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingRegion(k) == 'F') %only look at incorrect rule 1 frontal trials
                tempIncorrectFrontalBetty(zfi,1:1001) = bettyGoodTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                tempIncorrectFrontalBetty(zfi,1002:1811) = bettyGoodTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                zfi = zfi + 1;
            elseif (bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingRegion(k) == 'P') %only look at correct rule 1 parietal trials
                tempCorrectParietalBetty(zpc,1:1001) = bettyGoodTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                tempCorrectParietalBetty(zpc,1002:1811) = bettyGoodTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                zpc = zpc + 1;
            elseif (bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingRegion(k) == 'P') %only look at incorrect rule 1 parietal trials
                tempIncorrectParietalBetty(zpi,1:1001) = bettyGoodTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                tempIncorrectParietalBetty(zpi,1002:1811) = bettyGoodTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                zpi = zpi + 1;
            end
        end
    end
end
toc

numCor = sum(bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.BehResp == 1); %number of correct trials
numInc = sum(bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.BehResp == 0); %number of incorrect trials
numRule1 = sum(bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.Rule == 1); %number of rule 1 trials
numRule2 = sum(bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.Rule == 2); %number of rule 2 trials
numCorRule1 = sum((bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.BehResp == 1) & (bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.Rule == 1));