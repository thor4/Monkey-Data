data=rand(1,10)
 
r=eye(8)
 
r(1,:) %channel loc 6DR
 
 
% lfpData-sessionID-day-trialNum-chanLoc-chan-rule-resp [double floating point array]
 
 
v=[data,r(1,:)]
 
num_samples=500;
 
% data=zeros(61399,1200+3+1+1+2+11+2+2); %(N,1222)
data=rand(61399,1200+3+1+1+2+11+2+2); %(N,1222)
 
imagesc(data)
pause
 
% 
% for i = 1:num_samples
%     
%     [r,c]=size(voltage)
%    
%     
%     data(i,1:c)=voltage
%     
%     
%     
%     data(i,1201:1203)=session %[0,1,0]
%     data(i,1204) = day
%     data(i,1205) = trial num
%     data(i,1206:1207)=chanloc
%     data(i,1208:1218)=chan
%     data(i,1219:1220) =rule
%     data(i,1221:1222) =resp
%     
%     
% end
% 
 
chan=eye(11)
 
% {8B, 9L, 6DR, 8AD, vPFC, dPFC, LIP, MIP, PE, PG, PEC}

% crawling baseline period
baseline_period = (trial_info.CueOnset(29) - 50) - (trial_info.CueOnset(29) - 450);

lfp_data(1,trial_info.CueOnset(29)-450:trial_info.CueOnset(29)-51);

find(bad_trials ~= 0)

lfp_data(1,floor(trial_info.CueOnset(29))-450:floor(trial_info.CueOnset(29))-51);

data1 = zeros(400,10000);
data1 = data1(400,1:7531);

data=data(1:400,1:80166);

inc = incorrect(1:810,1:7531);

dataPAD = zeros(810,7531);
dataPAD(1:400,1:7531) = data(1:400,1:7531);

% horizontal concatenation of areas into regions

incorrectFrontal = horzcat(area8B, area9L, area6DR, area8AD, areavPFC, areadPFC);
incorrectParietal = horzcat(areaLIP, areaMIP, areaPE, areaPG, areaPEC);
incorrect = horzcat(incorrectFrontal, incorrectParietal);

correctFrontal = horzcat(area8B, area9L, area6DR, area8AD, areavPFC, areadPFC);
correctParietal = horzcat(areaLIP, areaMIP, areaPE, areaPG, areaPEC);
correct = horzcat(correctFrontal, correctParietal);
