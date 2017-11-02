clear
%load fft norm power
%replace location and BehResp acCordingly
%create indices of required bands
deltaIdx = find(cvPFCfreqCorR1>=0 & cvPFCfreqCorR1<4);
thetaIdx = find(cvPFCfreqCorR1>=4 & cvPFCfreqCorR1<8);
alphaIdx = find(cvPFCfreqCorR1>=8 & cvPFCfreqCorR1<12);
betaIdx = find(cvPFCfreqCorR1>=12 & cvPFCfreqCorR1<30);
lowgammaIdx = find(cvPFCfreqCorR1>=30 & cvPFCfreqCorR1<=55);
highgamma1Idx = find(cvPFCfreqCorR1>=65 & cvPFCfreqCorR1<=115);
highgamma2Idx = find(cvPFCfreqCorR1>=125 & cvPFCfreqCorR1<=175);
highgamma3Idx = find(cvPFCfreqCorR1>=185 & cvPFCfreqCorR1<=235);

%pull out the frequencies in each band then average over band
cvPFCpowCorR1delta = cvPFCpowCorR1(:,deltaIdx(1):deltaIdx(end));
cvPFCbandsCorR1(:,1) = mean(cvPFCpowCorR1delta(:,1:end),2);
cvPFCpowCorR1theta = cvPFCpowCorR1(:,thetaIdx(1):thetaIdx(end));
cvPFCbandsCorR1(:,2) = mean(cvPFCpowCorR1theta(:,1:end),2);
cvPFCpowCorR1alpha = cvPFCpowCorR1(:,alphaIdx(1):alphaIdx(end));
cvPFCbandsCorR1(:,3) = mean(cvPFCpowCorR1alpha(:,1:end),2);
cvPFCpowCorR1beta = cvPFCpowCorR1(:,betaIdx(1):betaIdx(end));
cvPFCbandsCorR1(:,4) = mean(cvPFCpowCorR1beta(:,1:end),2);
cvPFCpowCorR1lowgamma = cvPFCpowCorR1(:,lowgammaIdx(1):lowgammaIdx(end));
cvPFCbandsCorR1(:,5) = mean(cvPFCpowCorR1lowgamma(:,1:end),2);
cvPFCpowCorR1highgamma1 = cvPFCpowCorR1(:,highgamma1Idx(1):highgamma1Idx(end));
cvPFCbandsCorR1(:,6) = mean(cvPFCpowCorR1highgamma1(:,1:end),2);
cvPFCpowCorR1highgamma2 = cvPFCpowCorR1(:,highgamma2Idx(1):highgamma2Idx(end));
cvPFCbandsCorR1(:,7) = mean(cvPFCpowCorR1highgamma2(:,1:end),2);
cvPFCpowCorR1highgamma3 = cvPFCpowCorR1(:,highgamma3Idx(1):highgamma3Idx(end));
cvPFCbandsCorR1(:,8) = mean(cvPFCpowCorR1highgamma3(:,1:end),2);