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

figure
subplot(2,1,1);
plot(f, ScNormdB, f, SiNormdB, ':', 'LineWidth', 2);
title('LFP Slepian/Chronux Spectrum During Delay');
xlabel('Frequency (Hz)');
ylabel('Baseline-normalized Power (dB)');
legend('Correct Trials','Incorrect Trials');

subplot(2,1,2);
time = (0:250);
plot(time, 10*(log10(arpower_cAvgNorm)), time, 10*(log10(arpower_iAvgNorm)), ':', 'LineWidth', 2);
title('LFP AR/BSMART Spectrum During Delay');
xlabel('Frequency (Hz)');
ylabel('Baseline-normalized Power (dB)');
legend('Correct Trials','Incorrect Trials');


%various plots of power
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
area(f(deltaIdx(1):deltaIdx(end)),ScNormdB(deltaIdx(1):deltaIdx(end)), 'FaceColor', [.75 0 0]);
hold on
area(f(deltaIdx(end):thetaIdx(end)),ScNormdB(deltaIdx(end):thetaIdx(end)), 'FaceColor', [0 .75 0]);
hold on
area(f(thetaIdx(end):alphaIdx(end)),ScNormdB(thetaIdx(end):alphaIdx(end)), 'FaceColor', [0 0 .75]);
hold on
area(f(alphaIdx(end):betaIdx(end)),ScNormdB(alphaIdx(end):betaIdx(end)), 'FaceColor', [0.5 0 0]);
hold on
area(f(betaIdx(end):lowgammaIdx(end)),ScNormdB(betaIdx(end):lowgammaIdx(end)), 'FaceColor', [0 .5 0]);
hold on
area(f(highgamma1Idx(1):highgamma1Idx(end)),ScNormdB(highgamma1Idx(1):highgamma1Idx(end)), 'FaceColor', [0 0 .5]);
hold on
area(f(highgamma2Idx(1):highgamma2Idx(end)),ScNormdB(highgamma2Idx(1):highgamma2Idx(end)), 'FaceColor', [.25 0 0]);
hold on
area(f(highgamma3Idx(1):highgamma3Idx(end)),ScNormdB(highgamma3Idx(1):highgamma3Idx(end)), 'FaceColor', [0 .25 0]);

%incorrect power spectrum
subplot(1,2,2);
area(f(deltaIdx(1):deltaIdx(end)),SiNormdB(deltaIdx(1):deltaIdx(end)), 'FaceColor', [.75 0 0]);
hold on
area(f(deltaIdx(end):thetaIdx(end)),SiNormdB(deltaIdx(end):thetaIdx(end)), 'FaceColor', [0 .75 0]);
hold on
area(f(thetaIdx(end):alphaIdx(end)),SiNormdB(thetaIdx(end):alphaIdx(end)), 'FaceColor', [0 0 .75]);
hold on
area(f(alphaIdx(end):betaIdx(end)),SiNormdB(alphaIdx(end):betaIdx(end)), 'FaceColor', [0.5 0 0]);
hold on
area(f(betaIdx(end):lowgammaIdx(end)),SiNormdB(betaIdx(end):lowgammaIdx(end)), 'FaceColor', [0 .5 0]);
hold on
area(f(highgamma1Idx(1):highgamma1Idx(end)),SiNormdB(highgamma1Idx(1):highgamma1Idx(end)), 'FaceColor', [0 0 .5]);
hold on
area(f(highgamma2Idx(1):highgamma2Idx(end)),SiNormdB(highgamma2Idx(1):highgamma2Idx(end)), 'FaceColor', [.25 0 0]);
hold on
area(f(highgamma3Idx(1):highgamma3Idx(end)),SiNormdB(highgamma3Idx(1):highgamma3Idx(end)), 'FaceColor', [0 .25 0]);

%sum total energy in each power band
%Incorrect
chronuxDeltaI = sum((SiNormdB(deltaIdx(1):deltaIdx(end))));
chronuxThetaI = sum((SiNormdB(thetaIdx(1):thetaIdx(end))));
chronuxAlphaI = sum((SiNormdB(alphaIdx(1):alphaIdx(end))));
chronuxBetaI = sum((SiNormdB(betaIdx(1):betaIdx(end))));
chronuxLowGammaI = sum((SiNormdB(lowgammaIdx(1):lowgammaIdx(end))));
chronuxHighGamma1I = sum((SiNormdB(highgamma1Idx(1):highgamma1Idx(end))));
chronuxHighGamma2I = sum((SiNormdB(highgamma2Idx(1):highgamma2Idx(end))));
chronuxHighGamma3I = sum((SiNormdB(highgamma3Idx(1):highgamma3Idx(end))));

%Correct
chronuxDeltaC = sum((ScNormdB(deltaIdx(1):deltaIdx(end))));
chronuxThetaC = sum((ScNormdB(thetaIdx(1):thetaIdx(end))));
chronuxAlphaC = sum((ScNormdB(alphaIdx(1):alphaIdx(end))));
chronuxBetaC = sum((ScNormdB(betaIdx(1):betaIdx(end))));
chronuxLowGammaC = sum((ScNormdB(lowgammaIdx(1):lowgammaIdx(end))));
chronuxHighGamma1C = sum((ScNormdB(highgamma1Idx(1):highgamma1Idx(end))));
chronuxHighGamma2C = sum((ScNormdB(highgamma2Idx(1):highgamma2Idx(end))));
chronuxHighGamma3C = sum((ScNormdB(highgamma3Idx(1):highgamma3Idx(end))));

%bar graph
freq_cat = categorical({'delta', 'theta', 'alpha', 'beta', 'low gammma', 'high gamma 1', 'high gamma 2', 'high gamma 3'});
energy = [chronuxDeltaC chronuxDeltaI; chronuxThetaC chronuxThetaI; chronuxAlphaC chronuxAlphaI; chronuxBetaC chronuxBetaI; chronuxLowGammaC chronuxLowGammaI; chronuxHighGamma1C chronuxHighGamma1I; chronuxHighGamma2C chronuxHighGamma2I; chronuxHighGamma3C chronuxHighGamma3I];
bar(energy,freq_cat)
bar(energy)