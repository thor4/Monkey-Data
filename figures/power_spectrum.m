%power spectrum
% doc plot_vector
% figure
% subplot(2,1,1)
% plot_vector(ScDelay(:,10),fcDelay,'n',[],'b')
% hold on
% fill(theta,ScDelay(thetaIdx(1):thetaIdx(end),10),'r');
% hold on
% fill(fcDelay(1):fcDelay(5),ScDelay(1,10):ScDelay(length(fcDelay(1):fcDelay(5)),10),'r')
% subplot(2,1,2)
% plot_vector(SiDelay(:,10),fiDelay,'n',[],'g')
deltaIdx = find(f>=0 & f<4);
thetaIdx = find(f>=4 & f<8);
alphaIdx = find(f>=8 & f<12);
betaIdx = find(f>=12 & f<30);
lowgammaIdx = find(f>=30 & f<=50);
highgamma1Idx = find(f>=70 & f<=110);
highgamma2Idx = find(f>=130 & f<=170);
highgamma3Idx = find(f>=190 & f<=230);

%correct power spectrum
figure
subplot(1,2,1);
area(f(deltaIdx(1):deltaIdx(end)),ScNorm(deltaIdx(1):deltaIdx(end)), 'FaceColor', [.75 0 0]);
hold on
area(f(deltaIdx(end):thetaIdx(end)),ScNorm(deltaIdx(end):thetaIdx(end)), 'FaceColor', [0 .75 0]);
hold on
area(f(thetaIdx(end):alphaIdx(end)),ScNorm(thetaIdx(end):alphaIdx(end)), 'FaceColor', [0 0 .75]);
hold on
area(f(alphaIdx(end):betaIdx(end)),ScNorm(alphaIdx(end):betaIdx(end)), 'FaceColor', [0.5 0 0]);
hold on
area(f(betaIdx(end):lowgammaIdx(end)),ScNorm(betaIdx(end):lowgammaIdx(end)), 'FaceColor', [0 .5 0]);
hold on
area(f(highgamma1Idx(1):highgamma1Idx(end)),ScNorm(highgamma1Idx(1):highgamma1Idx(end)), 'FaceColor', [0 0 .5]);
hold on
area(f(highgamma2Idx(1):highgamma2Idx(end)),ScNorm(highgamma2Idx(1):highgamma2Idx(end)), 'FaceColor', [.25 0 0]);
hold on
area(f(highgamma3Idx(1):highgamma3Idx(end)),ScNorm(highgamma3Idx(1):highgamma3Idx(end)), 'FaceColor', [0 .25 0]);

%incorrect power spectrum
subplot(1,2,2);
area(f(deltaIdx(1):deltaIdx(end)),SiNorm(deltaIdx(1):deltaIdx(end)), 'FaceColor', [.75 0 0]);
hold on
area(f(deltaIdx(end):thetaIdx(end)),SiNorm(deltaIdx(end):thetaIdx(end)), 'FaceColor', [0 .75 0]);
hold on
area(f(thetaIdx(end):alphaIdx(end)),SiNorm(thetaIdx(end):alphaIdx(end)), 'FaceColor', [0 0 .75]);
hold on
area(f(alphaIdx(end):betaIdx(end)),SiNorm(alphaIdx(end):betaIdx(end)), 'FaceColor', [0.5 0 0]);
hold on
area(f(betaIdx(end):lowgammaIdx(end)),SiNorm(betaIdx(end):lowgammaIdx(end)), 'FaceColor', [0 .5 0]);
hold on
area(f(highgamma1Idx(1):highgamma1Idx(end)),SiNorm(highgamma1Idx(1):highgamma1Idx(end)), 'FaceColor', [0 0 .5]);
hold on
area(f(highgamma2Idx(1):highgamma2Idx(end)),SiNorm(highgamma2Idx(1):highgamma2Idx(end)), 'FaceColor', [.25 0 0]);
hold on
area(f(highgamma3Idx(1):highgamma3Idx(end)),SiNorm(highgamma3Idx(1):highgamma3Idx(end)), 'FaceColor', [0 .25 0]);

