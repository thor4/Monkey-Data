clear
load('good_stable_trials-betty.mat')
load('good_stable_trials-clark.mat')

%M1 8B, 9L, dPFC, vPFC, LIP, MIP, PEC, PG
m1C8b=1; m1C9l=1; m1Cdpfc=1; m1Cvpfc=1; m1Clip=1; m1Cmip=1; m1Cpec=1; m1Cpg=1;
m1I8b=1; m1I9l=1; m1Idpfc=1; m1Ivpfc=1; m1Ilip=1; m1Imip=1; m1Ipec=1; m1Ipg=1;

%M2 6DR, 8AD, 8B, dPFC, LIP, PE, PEC, PG
m2C6dr=1; m2C8ad=1; m2C8b=1; m2Cdpfc=1; m2Clip=1; m2Cpe=1; m2Cpec=1; m2Cpg=1;
m2I6dr=1; m2I8ad=1; m2I8b=1; m2Idpfc=1; m2Ilip=1; m2Ipe=1; m2Ipec=1; m2Ipg=1;

tic
%betty
fieldsBettyDays = fieldnames(bettyGoodStableTrials);
for i = 1:numel(fieldsBettyDays)
    fieldsBettyTrials = fieldnames(bettyGoodStableTrials.(fieldsBettyDays{i}));
    recordingArea = bettyGoodStableTrials.(fieldsBettyDays{i}).recording_info.area;
    for j = 1:numel(fieldsBettyTrials)-2
        for k = 1:numel(recordingArea)
            if (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('6DR')) %only look at correct rule 1 6DR trials
                m2CorrectRule16DR(m2C6dr,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2CorrectRule16DR(m2C6dr,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2C6dr = m2C6dr + 1;
                monkey(2).day(i).CorrectRule1.a6DR = m2CorrectRule16DR;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('6DR')) %only look at incorrect rule 1 6DR trials
                m2IncorrectRule16DR(m2I6dr,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2IncorrectRule16DR(m2I6dr,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2I6dr = m2I6dr + 1;
                monkey(2).day(i).IncorrectRule1.a6DR = m2IncorrectRule16DR;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('8AD')) %only look at correct rule 1 8AD trials
                m2CorrectRule18AD(m2C8ad,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2CorrectRule18AD(m2C8ad,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2C8ad = m2C8ad + 1;
                monkey(2).day(i).CorrectRule1.a8AD = m2CorrectRule18AD;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('8AD')) %only look at incorrect rule 1 8AD trials
                m2IncorrectRule18AD(m2I8ad,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2IncorrectRule18AD(m2I8ad,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2I8ad = m2I8ad + 1;
                monkey(2).day(i).IncorrectRule1.a8AD = m2IncorrectRule18AD;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('8B')) %only look at correct rule 1 8B trials
                m2CorrectRule18B(m2C8b,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2CorrectRule18B(m2C8b,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2C8b = m2C8b + 1;
                monkey(2).day(i).CorrectRule1.a8B = m2CorrectRule18B;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('8B')) %only look at incorrect rule 1 8B trials
                m2IncorrectRule18B(m2I8b,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2IncorrectRule18B(m2I8b,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2I8b = m2I8b + 1;
                monkey(2).day(i).IncorrectRule1.a8B = m2IncorrectRule18B;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('dPFC')) %only look at correct rule 1 dPFC trials
                m2CorrectRule1dPFC(m2Cdpfc,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2CorrectRule1dPFC(m2Cdpfc,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2Cdpfc = m2Cdpfc + 1;
                monkey(2).day(i).CorrectRule1.adPFC = m2CorrectRule1dPFC;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('dPFC')) %only look at incorrect rule 1 dPFC trials
                m2IncorrectRule1dPFC(m2Idpfc,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2IncorrectRule1dPFC(m2Idpfc,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2Idpfc = m2Idpfc + 1;
                monkey(2).day(i).IncorrectRule1.adPFC = m2IncorrectRule1dPFC;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('LIP')) %only look at correct rule 1 LIP trials
                m2CorrectRule1LIP(m2Clip,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2CorrectRule1LIP(m2Clip,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2Clip = m2Clip + 1;
                monkey(2).day(i).CorrectRule1.aLIP = m2CorrectRule1LIP;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('LIP')) %only look at incorrect rule 1 LIP trials
                m2IncorrectRule1LIP(m2Ilip,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2IncorrectRule1LIP(m2Ilip,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2Ilip = m2Ilip + 1;
                monkey(2).day(i).IncorrectRule1.aLIP = m2IncorrectRule1LIP;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('PE')) %only look at correct rule 1 PE trials
                m2CorrectRule1PE(m2Cpe,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2CorrectRule1PE(m2Cpe,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2Cpe = m2Cpe + 1;
                monkey(2).day(i).CorrectRule1.aPE = m2CorrectRule1PE;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('PE')) %only look at incorrect rule 1 PE trials
                m2IncorrectRule1PE(m2Ipe,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2IncorrectRule1PE(m2Ipe,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2Ipe = m2Ipe + 1;
                monkey(2).day(i).IncorrectRule1.aPE = m2IncorrectRule1PE;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('PEC')) %only look at correct rule 1 PEC trials
                m2CorrectRule1PEC(m2Cpec,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2CorrectRule1PEC(m2Cpec,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2Cpec = m2Cpec + 1;
                monkey(2).day(i).CorrectRule1.aPEC = m2CorrectRule1PEC;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('PEC')) %only look at incorrect rule 1 PEC trials
                m2IncorrectRule1PEC(m2Ipec,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2IncorrectRule1PEC(m2Ipec,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2Ipec = m2Ipec + 1;
                monkey(2).day(i).IncorrectRule1.aPEC = m2IncorrectRule1PEC;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 1) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('PG')) %only look at correct rule 1 PG trials
                m2CorrectRule1PG(m2Cpg,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2CorrectRule1PG(m2Cpg,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2Cpg = m2Cpg + 1;
                monkey(2).day(i).CorrectRule1.aPG = m2CorrectRule1PG;
            elseif (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == 0) && ...
                    (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == 1) && ...
                    (recordingArea{k} == string('PG')) %only look at incorrect rule 1 PG trials
                m2IncorrectRule1PG(m2Ipg,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                m2IncorrectRule1PG(m2Ipg,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                m2Ipg = m2Ipg + 1;
                monkey(2).day(i).IncorrectRule1.aPG = m2IncorrectRule1PG;
            end
        end
    end
    %reset the area variables so they get re-created in next iteration
    clearvars m2CorrectRule16DR m2CorrectRule18B m2CorrectRule1dPFC ...
        m2CorrectRule18AD m2CorrectRule1LIP m2CorrectRule1PE ...
        m2CorrectRule1PEC m2CorrectRule1PG  m2IncorrectRule16DR ...
        m2IncorrectRule18B m2IncorrectRule1dPFC m2IncorrectRule18AD ...
        m2IncorrectRule1LIP m2IncorrectRule1PE m2IncorrectRule1PEC ...
        m2IncorrectRule1PG
    %reset M2 6DR, 8AD, 8B, dPFC, LIP, PE, PEC, PG
    m2C6dr=1; m2C8ad=1; m2C8b=1; m2Cdpfc=1; m2Clip=1; m2Cpe=1; m2Cpec=1; m2Cpg=1;
    m2I6dr=1; m2I8ad=1; m2I8b=1; m2Idpfc=1; m2Ilip=1; m2Ipe=1; m2Ipec=1; m2Ipg=1;
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
        for k = 1:numel(fieldsClarkTrials)-2
            for l = 1:numel(recordingArea)
                if (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('8B')) %only look at correct rule 1 8B trials
                    m1CorrectRule18B(m1C8b,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1CorrectRule18B(m1C8b,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1C8b = m1C8b + 1;
                    monkey(1).day(i).CorrectRule1.a8B = m1CorrectRule18B;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('8B')) %only look at incorrect rule 1 8B trials
                    m1IncorrectRule18B(m1I8b,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1IncorrectRule18B(m1I8b,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1I8b = m1I8b + 1;
                    monkey(1).day(i).IncorrectRule1.a8B = m1IncorrectRule18B;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('9L')) %only look at correct rule 1 9L trials
                    m1CorrectRule19L(m1C9l,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1CorrectRule19L(m1C9l,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1C9l = m1C9l + 1;
                    monkey(1).day(i).CorrectRule1.a9L = m1CorrectRule19L;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('9L')) %only look at incorrect rule 1 9L trials
                    m1IncorrectRule19L(m1I9l,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1IncorrectRule19L(m1I9l,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1I9l = m1I9l + 1;
                    monkey(1).day(i).IncorrectRule1.a9L = m1IncorrectRule19L;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('dPFC')) %only look at correct rule 1 dPFC trials
                    m1CorrectRule1dPFC(m1Cdpfc,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1CorrectRule1dPFC(m1Cdpfc,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Cdpfc = m1Cdpfc + 1;
                    monkey(1).day(i).CorrectRule1.adPFC = m1CorrectRule1dPFC;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('dPFC')) %only look at incorrect rule 1 dPFC trials
                    m1IncorrectRule1dPFC(m1Idpfc,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1IncorrectRule1dPFC(m1Idpfc,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Idpfc = m1Idpfc + 1;
                    monkey(1).day(i).IncorrectRule1.adPFC = m1IncorrectRule1dPFC;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('vPFC')) %only look at correct rule 1 vPFC trials
                    m1CorrectRule1vPFC(m1Cvpfc,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1CorrectRule1vPFC(m1Cvpfc,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Cvpfc = m1Cvpfc + 1;
                    monkey(1).day(i).CorrectRule1.avPFC = m1CorrectRule1vPFC;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('vPFC')) %only look at incorrect rule 1 vPFC trials
                    m1IncorrectRule1vPFC(m1Ivpfc,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1IncorrectRule1vPFC(m1Ivpfc,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Ivpfc = m1Ivpfc + 1;
                    monkey(1).day(i).IncorrectRule1.avPFC = m1IncorrectRule1vPFC;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('LIP')) %only look at correct rule 1 LIP trials
                    m1CorrectRule1LIP(m1Clip,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1CorrectRule1LIP(m1Clip,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Clip = m1Clip + 1;
                    monkey(1).day(i).CorrectRule1.aLIP = m1CorrectRule1LIP;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('LIP')) %only look at incorrect rule 1 LIP trials
                    m1IncorrectRule1LIP(m1Ilip,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1IncorrectRule1LIP(m1Ilip,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Ilip = m1Ilip + 1;
                    monkey(1).day(i).IncorrectRule1.aLIP = m1IncorrectRule1LIP;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('MIP')) %only look at correct rule 1 MIP trials
                    m1CorrectRule1MIP(m1Cmip,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1CorrectRule1MIP(m1Cmip,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Cmip = m1Cmip + 1;
                    monkey(1).day(i).CorrectRule1.aMIP = m1CorrectRule1MIP;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('MIP')) %only look at incorrect rule 1 MIP trials
                    m1IncorrectRule1MIP(m1Imip,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1IncorrectRule1MIP(m1Imip,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Imip = m1Imip + 1;
                    monkey(1).day(i).IncorrectRule1.aMIP = m1IncorrectRule1MIP;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('PEC')) %only look at correct rule 1 PEC trials
                    m1CorrectRule1PEC(m1Cpec,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1CorrectRule1PEC(m1Cpec,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Cpec = m1Cpec + 1;
                    monkey(1).day(i).CorrectRule1.aPEC = m1CorrectRule1PEC;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('PEC')) %only look at incorrect rule 1 PEC trials
                    m1IncorrectRule1PEC(m1Ipec,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1IncorrectRule1PEC(m1Ipec,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Ipec = m1Ipec + 1;
                    monkey(1).day(i).IncorrectRule1.aPEC = m1IncorrectRule1PEC;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('PG')) %only look at correct rule 1 PG trials
                    m1CorrectRule1PG(m1Cpg,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1CorrectRule1PG(m1Cpg,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Cpg = m1Cpg + 1;
                    monkey(1).day(i).CorrectRule1.aPG = m1CorrectRule1PG;
                elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
                        (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
                        (recordingArea{l} == string('PG')) %only look at incorrect rule 1 PG trials
                    m1IncorrectRule1PG(m1Ipg,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
                    m1IncorrectRule1PG(m1Ipg,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
                    m1Ipg = m1Ipg + 1;
                    monkey(1).day(i).IncorrectRule1.aPG = m1IncorrectRule1PG;
                end
            end
        end
    end
    %reset the area variables so they get re-created in next iteration
    clearvars m1CorrectRule19L m1CorrectRule18B m1CorrectRule1dPFC ...
        m1CorrectRule1vPFC m1CorrectRule1LIP m1CorrectRule1MIP ...
        m1CorrectRule1PEC m1CorrectRule1PG  m1IncorrectRule19L ...
        m1IncorrectRule18B m1IncorrectRule1dPFC m1IncorrectRule1vPFC ...
        m1IncorrectRule1LIP m1IncorrectRule1MIP m1IncorrectRule1PEC ...
        m1IncorrectRule1PG
    %reset M1 8B, 9L, dPFC, vPFC, LIP, MIP, PEC, PG
    m1C8b=1; m1C9l=1; m1Cdpfc=1; m1Cvpfc=1; m1Clip=1; m1Cmip=1; m1Cpec=1; m1Cpg=1;
    m1I8b=1; m1I9l=1; m1Idpfc=1; m1Ivpfc=1; m1Ilip=1; m1Imip=1; m1Ipec=1; m1Ipg=1;
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