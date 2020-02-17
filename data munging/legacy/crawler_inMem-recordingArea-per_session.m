clear
load('good_stable_trials-betty.mat')
load('good_stable_trials-clark.mat')

%initialize counters for each monkey + behavioral response (behresp)
m1C=1; m1I=1; m2C=1; m2I=1; 

chan = 'chan%d%s';

tic
%betty
fieldsBettyDays = fieldnames(bettyGoodStableTrials);
for i = 1:numel(fieldsBettyDays)
    fieldsBettyTrials = fieldnames(bettyGoodStableTrials.(fieldsBettyDays{i}));
    recordingArea = bettyGoodStableTrials.(fieldsBettyDays{i}).recording_info.area;
    for k = 1:numel(recordingArea)
        channame = sprintf(chan, k, recordingArea{k});
        for j = 1:numel(fieldsBettyTrials)-2
            if (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) %look at correct rule 1
                m2CorrectRule1(m2C,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2CorrectRule1(m2C,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2C = m2C + 1;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) %look at incorrect rule 1
                m2IncorrectRule1(m2I,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2IncorrectRule1(m2I,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2I = m2I + 1;
            end            
        end
        monkey(2).day(i).correct.(channame) = m2CorrectRule1;
        monkey(2).day(i).incorrect.(channame) = m2IncorrectRule1;
        %reset the BehResp variables so they get re-created in next iteration
        clearvars m2CorrectRule1 m2IncorrectRule1
        %reset counters
        m2C=1; m2I=1;
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
        recordingArea = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).recording_info.area;
        for l = 1:numel(recordingArea)
            channame = sprintf(chan, l, recordingArea{l});
            for k = 1:numel(fieldsClarkTrials)-2
                if (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) %look at correct rule 1
                    m1CorrectRule1(m1C,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1CorrectRule1(m1C,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1C = m1C + 1;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) %look at incorrect rule 1
                    m1IncorrectRule1(m1I,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1IncorrectRule1(m1I,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1I = m1I + 1;
                end
            end
            %ensure session is Rule 1
            if (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1)
                monkey(1).day(i).correct.(channame) = m1CorrectRule1;
                monkey(1).day(i).incorrect.(channame) = m1IncorrectRule1;
            end
            %reset the BehResp variables so they get re-created in next iteration
            clearvars m1CorrectRule1 m1IncorrectRule1
            %reset counters
            m1C=1; m1I=1;
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