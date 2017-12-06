clear
%load fft norm power
%replace location and BehResp acIncdingly
%create indices of required bands
deltaIdx = find(bdPFCfreqIncR1>=0 & bdPFCfreqIncR1<4);
thetaIdx = find(bdPFCfreqIncR1>=4 & bdPFCfreqIncR1<8);
alphaIdx = find(bdPFCfreqIncR1>=8 & bdPFCfreqIncR1<12);
betaIdx = find(bdPFCfreqIncR1>=12 & bdPFCfreqIncR1<30);
lowgammaIdx = find(bdPFCfreqIncR1>=30 & bdPFCfreqIncR1<=55);
highgamma1Idx = find(bdPFCfreqIncR1>=65 & bdPFCfreqIncR1<=115);
highgamma2Idx = find(bdPFCfreqIncR1>=125 & bdPFCfreqIncR1<=175);
highgamma3Idx = find(bdPFCfreqIncR1>=185 & bdPFCfreqIncR1<=235);

%pull out the frequencies in each band then average over band
bdPFCpowIncR1delta = bdPFCpowIncR1(:,deltaIdx(1):deltaIdx(end));
bdPFCbandsIncR1(:,1) = mean(bdPFCpowIncR1delta(:,1:end),2);
bdPFCpowIncR1theta = bdPFCpowIncR1(:,thetaIdx(1):thetaIdx(end));
bdPFCbandsIncR1(:,2) = mean(bdPFCpowIncR1theta(:,1:end),2);
bdPFCpowIncR1alpha = bdPFCpowIncR1(:,alphaIdx(1):alphaIdx(end));
bdPFCbandsIncR1(:,3) = mean(bdPFCpowIncR1alpha(:,1:end),2);
bdPFCpowIncR1beta = bdPFCpowIncR1(:,betaIdx(1):betaIdx(end));
bdPFCbandsIncR1(:,4) = mean(bdPFCpowIncR1beta(:,1:end),2);
bdPFCpowIncR1lowgamma = bdPFCpowIncR1(:,lowgammaIdx(1):lowgammaIdx(end));
bdPFCbandsIncR1(:,5) = mean(bdPFCpowIncR1lowgamma(:,1:end),2);
bdPFCpowIncR1highgamma1 = bdPFCpowIncR1(:,highgamma1Idx(1):highgamma1Idx(end));
bdPFCbandsIncR1(:,6) = mean(bdPFCpowIncR1highgamma1(:,1:end),2);
bdPFCpowIncR1highgamma2 = bdPFCpowIncR1(:,highgamma2Idx(1):highgamma2Idx(end));
bdPFCbandsIncR1(:,7) = mean(bdPFCpowIncR1highgamma2(:,1:end),2);
bdPFCpowIncR1highgamma3 = bdPFCpowIncR1(:,highgamma3Idx(1):highgamma3Idx(end));
bdPFCbandsIncR1(:,8) = mean(bdPFCpowIncR1highgamma3(:,1:end),2);