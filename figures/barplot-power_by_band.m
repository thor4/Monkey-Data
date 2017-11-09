%bar plot with error bars and significance markers, rename band each time
clear
%mean of frequency band across all trials (rows) in each band (columns)
%correct in first column, incorrect in second column
bDeltaRegionsR1mean(1,1) = mean(b6DRbandsCorR1sampled(:,1));
bDeltaRegionsR1mean(2,1) = mean(b8ADbandsCorR1sampled(:,1));
bDeltaRegionsR1mean(3,1) = mean(b8BbandsCorR1sampled(:,1));
bDeltaRegionsR1mean(4,1) = mean(bdPFCbandsCorR1sampled(:,1));
bDeltaRegionsR1mean(5,1) = mean(bLIPbandsCorR1sampled(:,1));
bDeltaRegionsR1mean(6,1) = mean(bPEbandsCorR1sampled(:,1));
bDeltaRegionsR1mean(7,1) = mean(bPECbandsCorR1sampled(:,1));
bDeltaRegionsR1mean(8,1) = mean(bPGbandsCorR1sampled(:,1));
bDeltaRegionsR1mean(1,2) = mean(b6DRbandsIncR1(:,1));
bDeltaRegionsR1mean(2,2) = mean(b8ADbandsIncR1(:,1));
bDeltaRegionsR1mean(3,2) = mean(b8BbandsIncR1(:,1));
bDeltaRegionsR1mean(4,2) = mean(bdPFCbandsIncR1(:,1));
bDeltaRegionsR1mean(5,2) = mean(bLIPbandsIncR1(:,1));
bDeltaRegionsR1mean(6,2) = mean(bPEbandsIncR1(:,1));
bDeltaRegionsR1mean(7,2) = mean(bPECbandsIncR1(:,1));
bDeltaRegionsR1mean(8,2) = mean(bPGbandsIncR1(:,1));

%confidence intervals across all trials (rows) in each band (columns)
%correct in first column, incorrect in second column
ci = 0.95 ; alpha = 1 - ci;
n = size(b6DRbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(1,1) = T_multiplier.*std(b6DRbandsCorR1sampled(:,1))./sqrt(n);
n = size(b8ADbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(2,1) = T_multiplier.*std(b8ADbandsCorR1sampled(:,1))./sqrt(n);
n = size(b8BbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(3,1) = T_multiplier.*std(b8BbandsCorR1sampled(:,1))./sqrt(n);
n = size(bdPFCbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(4,1) = T_multiplier.*std(bdPFCbandsCorR1sampled(:,1))./sqrt(n);
n = size(bLIPbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(5,1) = T_multiplier.*std(bLIPbandsCorR1sampled(:,1))./sqrt(n);
n = size(bPEbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(6,1) = T_multiplier.*std(bPEbandsCorR1sampled(:,1))./sqrt(n);
n = size(bPECbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(7,1) = T_multiplier.*std(bPECbandsCorR1sampled(:,1))./sqrt(n);
n = size(bPGbandsCorR1sampled,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(8,1) = T_multiplier.*std(bPGbandsCorR1sampled(:,1))./sqrt(n);
n = size(b6DRbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(1,2) = T_multiplier.*std(b6DRbandsIncR1(:,1))./sqrt(n);
n = size(b8ADbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(2,2) = T_multiplier.*std(b8ADbandsIncR1(:,1))./sqrt(n);
n = size(b8BbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(3,2) = T_multiplier.*std(b8BbandsIncR1(:,1))./sqrt(n);
n = size(bdPFCbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(4,2) = T_multiplier.*std(bdPFCbandsIncR1(:,1))./sqrt(n);
n = size(bLIPbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(5,2) = T_multiplier.*std(bLIPbandsIncR1(:,1))./sqrt(n);
n = size(bPEbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(6,2) = T_multiplier.*std(bPEbandsIncR1(:,1))./sqrt(n);
n = size(bPECbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(7,2) = T_multiplier.*std(bPECbandsIncR1(:,1))./sqrt(n);
n = size(bPGbandsIncR1,1); T_multiplier = tinv(1-alpha/2, n-1);
bDeltaRegionsR1ci95(8,2) = T_multiplier.*std(bPGbandsIncR1(:,1))./sqrt(n);

% ci95(1,:) = bDeltaRegionsR1ci95(2,:);
% ci95(2,:) = bDeltaRegionsR1ci95(8,:);
% y(1,:) = bDeltaRegionsR1mean(2,:);
% y(2,:) = bDeltaRegionsR1mean(8,:);

figure
[hBar hErrorbar] = barwitherr(bDeltaRegionsR1ci95,bDeltaRegionsR1mean);
%[hBar hErrorbar] = barwitherr(ci95,y);
%sigstar({'8B'{1},'8B'{2}})
% sigstar({[1,2], [1,3]})
% for k = 1:size(y,2)
%     b(k).CData = k;
% end
hBar(1).FaceColor = [0 .271 .937]; hBar(2).FaceColor = 'red';
hBar(1).LineStyle = 'none'; hBar(2).LineStyle = 'none';
xticks([])
% ax = gca
% ax.YGrid = 'off'
axis square
set(gca,'box','off')
set(gca,'FontSize',25);

newXticklabel = {'8B','9L','dPFC','vPFC','LIP','MIP','PEC','PG'};
%set(gca,'XtickLabel',newXticklabel,'FontSize',25);
title('Clark Delta Band Delay-period .95 CI');
%xlabel('Frequency Bands');
ylabel('Power');
legend('Correct Trials','Incorrect Trials');

ctr2 = bsxfun(@plus, hb(2).XData, [hb(2).XOffset]');
hold on
plot(ctr2(1:2), [1 1]*bDeltaRegionsR1mean(1,2)*1.1, '-k', 'LineWidth',2)
plot(mean(ctr2(1:2)), bDeltaRegionsR1mean(1,2)*1.15, '*k')
hold off
