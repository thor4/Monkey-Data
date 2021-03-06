load('good_stable_trials-betty.mat')
load('good_stable_trials-clark.mat')


zfc = 1;
zfi = 1;
zpc = 1;
zpi = 1;

%betty
tic
fieldsBettyDays = fieldnames(bettyGoodStableTrials);
for i = 1:numel(fieldsBettyDays)
    fieldsBettyTrials = fieldnames(bettyGoodStableTrials.(fieldsBettyDays{i}));
    recordingRegion = bettyGoodStableTrials.(fieldsBettyDays{i}).recording_info.cortex;
    for j = 1:numel(fieldsBettyTrials)-2
        for k = 1:numel(recordingRegion)
            if (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingRegion(k) == 'F') %only look at correct rule 1 frontal trials
                correctFrontalBetty(zfc,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                correctFrontalBetty(zfc,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                zfc = zfc + 1;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingRegion(k) == 'F') %only look at incorrect rule 1 frontal trials
                incorrectFrontalBetty(zfi,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                incorrectFrontalBetty(zfi,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                zfi = zfi + 1;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingRegion(k) == 'P') %only look at correct rule 1 parietal trials
                correctParietalBetty(zpc,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                correctParietalBetty(zpc,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                zpc = zpc + 1;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingRegion(k) == 'P') %only look at incorrect rule 1 parietal trials
                incorrectParietalBetty(zpi,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                incorrectParietalBetty(zpi,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                zpi = zpi + 1;
            end
        end
    end
end
toc

%clark
tic
fieldsClarkDays = fieldnames(clarkGoodStableTrials);
for i = 1:numel(fieldsClarkDays)
    fieldsClarkSessions = fieldnames(clarkGoodStableTrials.(fieldsClarkDays{i}));
    for j = 1:numel(fieldsClarkSessions)
        fieldsClarkTrials = fieldnames(clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}));
        recordingRegion = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).recording_info.cortex;
        for k = 1:numel(fieldsClarkTrials)-2
            for l = 1:numel(recordingRegion)
                if (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingRegion(l) == 'F') %only look at correct rule 1 frontal trials
                    correctFrontalClark(zfc,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    correctFrontalClark(zfc,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    zfc = zfc + 1;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingRegion(l) == 'F') %only look at incorrect rule 1 frontal trials
                    incorrectFrontalClark(zfi,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    incorrectFrontalClark(zfi,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    zfi = zfi + 1;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingRegion(l) == 'P') %only look at correct rule 1 parietal trials
                    correctParietalClark(zpc,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    correctParietalClark(zpc,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    zpc = zpc + 1;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingRegion(l) == 'P') %only look at incorrect rule 1 parietal trials
                    incorrectParietalClark(zpi,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    incorrectParietalClark(zpi,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    zpi = zpi + 1;
                end
            end
        end
    end
end
 
toc

%betty counts
numCor = sum(bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp == 1); %number of correct trials
numInc = sum(bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp == 0); %number of incorrect trials
numRule1 = sum(bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule == 1); %number of rule 1 trials
numRule2 = sum(bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule == 2); %number of rule 2 trials
numCorRule1 = sum((bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp == 1) & (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule == 1)); %number of correct, rule 1 trials

%clark counts
numCor = sum(clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp == 1); %number of correct trials
numInc = sum(clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp == 0); %number of incorrect trials
numRule1 = sum(clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule == 1); %number of rule 1 trials
numRule2 = sum(clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule == 2); %number of rule 2 trials
numCorRule1 = sum((clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp == 1) & (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule == 1)); %number of correct, rule 1 trials