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
 
% {8B, 9L, LIP, MIP, PG, vPFC 6DR, 8AD, PE, PEC, dPFC}
