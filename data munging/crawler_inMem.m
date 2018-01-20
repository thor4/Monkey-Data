load('good_stable_trials-betty.mat')
load('good_stable_trials-clark.mat')


zfc = 1;
zfi = 1;
zpc = 1;
zpi = 1;

%betty
tic
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

%clark
tic
fieldsClarkDays = fieldnames(clarkGoodTrials);
for i = 1:numel(fieldsClarkDays)
    fieldsClarkSessions = fieldnames(clarkGoodTrials.(fieldsClarkDays{i}));
    for j = 1:numel(fieldsClarkSessions)
        fieldsClarkTrials = fieldnames(clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}));
        recordingRegion = clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).recording_info.cortex;
        for k = 1:numel(fieldsClarkTrials)-2
            for l = 1:numel(recordingRegion)
                if (clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingRegion(l) == 'F') %only look at correct rule 1 frontal trials
                    correctFrontalClark(zfc,1:1001) = clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    correctFrontalClark(zfc,1002:1811) = clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    zfc = zfc + 1;
                elseif (clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingRegion(l) == 'F') %only look at incorrect rule 1 frontal trials
                    incorrectFrontalClark(zfi,1:1001) = clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    incorrectFrontalClark(zfi,1002:1811) = clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    zfi = zfi + 1;
                elseif (clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingRegion(l) == 'P') %only look at correct rule 1 parietal trials
                    correctParietalClark(zpc,1:1001) = clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    correctParietalClark(zpc,1002:1811) = clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    zpc = zpc + 1;
                elseif (clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingRegion(l) == 'P') %only look at incorrect rule 1 parietal trials
                    incorrectParietalClark(zpi,1:1001) = clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    incorrectParietalClark(zpi,1002:1811) = clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    zpi = zpi + 1;
                end
            end
        end
    end
end
 
toc

%betty counts
numCor = sum(bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.BehResp == 1); %number of correct trials
numInc = sum(bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.BehResp == 0); %number of incorrect trials
numRule1 = sum(bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.Rule == 1); %number of rule 1 trials
numRule2 = sum(bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.Rule == 2); %number of rule 2 trials
numCorRule1 = sum((bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.BehResp == 1) & (bettyGoodTrials.(fieldsBettyDays{i}).new_trial_info.Rule == 1)); %number of correct, rule 1 trials

%clark counts
numCor = sum(clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp == 1); %number of correct trials
numInc = sum(clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp == 0); %number of incorrect trials
numRule1 = sum(clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule == 1); %number of rule 1 trials
numRule2 = sum(clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule == 2); %number of rule 2 trials
numCorRule1 = sum((clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp == 1) & (clarkGoodTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule == 1)); %number of correct, rule 1 trials