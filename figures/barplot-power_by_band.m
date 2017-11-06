%bar plot with error bars and significance markers, rename band each time
clear
%mean of frequency band across all trials (rows) in each band (columns)
%correct in first column, incorrect in second column
cDeltaRegionsR1mean(1,1) = mean(c8BbandsCorR1sampled(:,1));
cDeltaRegionsR1mean(2,1) = mean(c9LbandsCorR1sampled(:,1));
cDeltaRegionsR1mean(3,1) = mean(cdPFCbandsCorR1sampled(:,1));
cDeltaRegionsR1mean(4,1) = mean(cvPFCbandsCorR1sampled(:,1));
cDeltaRegionsR1mean(5,1) = mean(cLIPbandsCorR1sampled(:,1));
cDeltaRegionsR1mean(6,1) = mean(cMIPbandsCorR1sampled(:,1));
cDeltaRegionsR1mean(7,1) = mean(cPECbandsCorR1sampled(:,1));
cDeltaRegionsR1mean(8,1) = mean(cPGbandsCorR1sampled(:,1));
cDeltaRegionsR1mean(1,2) = mean(c8BbandsIncR1(:,1));
cDeltaRegionsR1mean(2,2) = mean(c9LbandsIncR1(:,1));
cDeltaRegionsR1mean(3,2) = mean(cdPFCbandsIncR1(:,1));
cDeltaRegionsR1mean(4,2) = mean(cvPFCbandsIncR1(:,1));
cDeltaRegionsR1mean(5,2) = mean(cLIPbandsIncR1(:,1));
cDeltaRegionsR1mean(6,2) = mean(cMIPbandsIncR1(:,1));
cDeltaRegionsR1mean(7,2) = mean(cPECbandsIncR1(:,1));
cDeltaRegionsR1mean(8,2) = mean(cPGbandsIncR1(:,1));

%confidence intervals across all trials (rows) in each band (columns)
%correct in first column, incorrect in second column
ci = 0.95 ; alpha = 1 - ci;
n = size(c8BbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(1,1) = T_multiplier.*std(c8BbandsCorR1sampled(:,1))./sqrt(n);
n = size(c9LbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(2,1) = T_multiplier.*std(c9LbandsCorR1sampled(:,1))./sqrt(n);
n = size(cdPFCbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(3,1) = T_multiplier.*std(cdPFCbandsCorR1sampled(:,1))./sqrt(n);
n = size(cvPFCbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(4,1) = T_multiplier.*std(cvPFCbandsCorR1sampled(:,1))./sqrt(n);
n = size(cLIPbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(5,1) = T_multiplier.*std(cLIPbandsCorR1sampled(:,1))./sqrt(n);
n = size(cMIPbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(6,1) = T_multiplier.*std(cMIPbandsCorR1sampled(:,1))./sqrt(n);
n = size(cPECbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(7,1) = T_multiplier.*std(cPECbandsCorR1sampled(:,1))./sqrt(n);
n = size(cPGbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(8,1) = T_multiplier.*std(cPGbandsCorR1sampled(:,1))./sqrt(n);
n = size(c8BbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(1,2) = T_multiplier.*std(c8BbandsIncR1(:,1))./sqrt(n);
n = size(c9LbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(2,2) = T_multiplier.*std(c9LbandsIncR1(:,1))./sqrt(n);
n = size(cdPFCbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(3,2) = T_multiplier.*std(cdPFCbandsIncR1(:,1))./sqrt(n);
n = size(cvPFCbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(4,2) = T_multiplier.*std(cvPFCbandsIncR1(:,1))./sqrt(n);
n = size(cLIPbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(5,2) = T_multiplier.*std(cLIPbandsIncR1(:,1))./sqrt(n);
n = size(cMIPbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(6,2) = T_multiplier.*std(cMIPbandsIncR1(:,1))./sqrt(n);
n = size(cPECbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(7,2) = T_multiplier.*std(cPECbandsIncR1(:,1))./sqrt(n);
n = size(cPGbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
cDeltaRegionsR1ci95(8,2) = T_multiplier.*std(cPGbandsIncR1(:,1))./sqrt(n);

figure
hb=barwitherr(cDeltaRegionsR1ci95,cDeltaRegionsR1mean);
%sigstar({'8B'{1},'8B'{2}})
% sigstar({[1,2], [1,3]})
newXticklabel = {'8B','9L','dPFC','vPFC','LIP','MIP','PEC','PG'};
set(gca,'XtickLabel',newXticklabel,'FontSize',25);
title('Clark Delta Band Delay-period .95 CI');
%xlabel('Frequency Bands');
ylabel('Power');
legend('Correct Trials','Incorrect Trials');

ctr2 = bsxfun(@plus, hb(2).XData, [hb(2).XOffset]');
hold on
plot(ctr2(1:2), [1 1]*cDeltaRegionsR1mean(1,2)*1.1, '-k', 'LineWidth',2)
plot(mean(ctr2(1:2)), cDeltaRegionsR1mean(1,2)*1.15, '*k')
hold off


https://www.mathworks.com/matlabcentral/answers/203782-how-to-determine-the-x-axis-positions-of-the-bars-in-a-grouped-bar-chart
https://www.mathworks.com/matlabcentral/answers/175193-creating-sigstar-in-bar-graph