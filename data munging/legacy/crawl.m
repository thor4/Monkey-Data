%BehResp can only be 0 for incorrect or 1 for correct
%Rule can only be 1 for identity rule or 2 for location rule
%region can only be 'F' for frontal or 'P' for parietal
%area can only be '9L','8B','6DR','8AD','vPFC','dPFC','LIP','MIP','PE','PG'
%, or 'PEC'

%doc: https://www.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html

function [M1,M2] = crawl(data,varargin)
    p = inputParser; %create inputParser object to check inputs
    %define default optional parameter values
    defaultBehResp = 1; %correct trials
    defaultRule = 1; %identity rule
    defaultRegion = 'no'; %no region specified
    defaultArea = 'no'; %no area specified
    %add required input + optional parameter values & verify datatype
    addRequired(p,'data',@isstruct);
    addParameter(p,'BehResp',defaultBehResp,@isnumeric)
    addParameter(p,'Rule',defaultRule,@isnumeric)
    addParameter(p,'Region',defaultRegion,@ischar)
    addParameter(p,'Area',defaultArea,@ischar)
    parse(p,data,varargin{:}) %parse all the inputs
    %report parsed inputs back to user for confirmation
    disp(['Crawling: ',p.Results.data])
    if ~isempty(p.UsingDefaults) %check if using any default values
        disp('Using defaults: ')
        disp(p.UsingDefaults)
    end
    z = 1;
    M1=0;
    %betty
    fieldsBettyDays = fieldnames(bettyGoodStableTrials);
    for i = 1:numel(fieldsBettyDays)
        fieldsBettyTrials = fieldnames(bettyGoodStableTrials.(fieldsBettyDays{i}));
        recordingRegion = bettyGoodStableTrials.(fieldsBettyDays{i}).recording_info.cortex;
        recordingArea = bettyGoodStableTrials.(fieldsBettyDays{i}).recording_info.area;
        for j = 1:numel(fieldsBettyTrials)-2
            if BehResp && Rule && region %provide BehResp,Rule,region slice
                for k = 1:numel(recordingRegion)                
                    if (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == BehResp) && ...
                            (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == Rule) && ...
                            (recordingRegion(k) == region) %only look at trials from requested region
                        M2(z,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                        M2(z,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                        z = z + 1;
                    end
                end
            elseif BehResp && Rule && area %provide BehResp,Rule,area slice
                for k = 1:numel(recordingArea)
                   if (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.BehResp(j) == BehResp) && ...
                            (bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.Rule(j) == Rule) && ...
                            (recordingArea(k) == area) %only look at trials from requested area
                        M2(z,1:1001) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,1:1001);
                        M2(z,1002:1811) = bettyGoodStableTrials.(fieldsBettyDays{i}).(fieldsBettyTrials{j})(k,bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+1:bettyGoodStableTrials.(fieldsBettyDays{i}).new_trial_info.CueOffset(j)+810);
                        z = z + 1;
                   end 
                end
            else
                print("Please enter either BehResp,Rule,Region combo or BehResp,Rule,area");
            end
        end
    end
end
    %clark
%     tic
%     fieldsClarkDays = fieldnames(clarkGoodStableTrials);
%     for i = 1:numel(fieldsClarkDays)
%         fieldsClarkSessions = fieldnames(clarkGoodStableTrials.(fieldsClarkDays{i}));
%         for j = 1:numel(fieldsClarkSessions)
%             fieldsClarkTrials = fieldnames(clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}));
%             recordingRegion = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).recording_info.cortex;
%             for k = 1:numel(fieldsClarkTrials)-2
%                 for l = 1:numel(recordingRegion)
%                     if (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
%                             (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
%                             (recordingRegion(l) == 'F') %only look at correct rule 1 frontal trials
%                         correctFrontalClark(zfc,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
%                         correctFrontalClark(zfc,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
%                         zfc = zfc + 1;
%                     elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
%                             (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
%                             (recordingRegion(l) == 'F') %only look at incorrect rule 1 frontal trials
%                         incorrectFrontalClark(zfi,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
%                         incorrectFrontalClark(zfi,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
%                         zfi = zfi + 1;
%                     elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 1) && ...
%                             (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
%                             (recordingRegion(l) == 'P') %only look at correct rule 1 parietal trials
%                         correctParietalClark(zpc,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
%                         correctParietalClark(zpc,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
%                         zpc = zpc + 1;
%                     elseif (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.BehResp(k) == 0) && ...
%                             (clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.Rule(k) == 1) && ...
%                             (recordingRegion(l) == 'P') %only look at incorrect rule 1 parietal trials
%                         incorrectParietalClark(zpi,1:1001) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,1:1001);
%                         incorrectParietalClark(zpi,1002:1811) = clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).(fieldsClarkTrials{k})(l,clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+1:clarkGoodStableTrials.(fieldsClarkDays{i}).(fieldsClarkSessions{j}).new_trial_info.CueOffset(k)+810);
%                         zpi = zpi + 1;
%                     end
%                 end
%             end
%         end
%     end
% end
% toc