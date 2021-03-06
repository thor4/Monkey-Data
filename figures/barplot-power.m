clear
%means
cLIPbandsCorR1samplemeans = mean(cLIPbandsCorR1samples,1);
cLIPbandsIncR1samplemeans = mean(cLIPbandsIncR1samples,1);
y = [cLIPbandsCorR1samplemeans(1) cLIPbandsIncR1samplemeans(1); ...
    cLIPbandsCorR1samplemeans(2) cLIPbandsIncR1samplemeans(2); ...
    cLIPbandsCorR1samplemeans(3) cLIPbandsIncR1samplemeans(3); ...
    cLIPbandsCorR1samplemeans(4) cLIPbandsIncR1samplemeans(4); ...
    cLIPbandsCorR1samplemeans(5) cLIPbandsIncR1samplemeans(5); ...
    cLIPbandsCorR1samplemeans(6) cLIPbandsIncR1samplemeans(6); ...
    cLIPbandsCorR1samplemeans(7) cLIPbandsIncR1samplemeans(7); ...
    cLIPbandsCorR1samplemeans(8) cLIPbandsIncR1samplemeans(8)];
%standard deviation then standard error of the mean
n = size(cLIPbandsCorR1samples,1);
cLIPbandsCorR1samplesem = std(cLIPbandsCorR1samples)/sqrt(n);
cLIPbandsIncR1samplesem = std(cLIPbandsIncR1samples)/sqrt(n);
sem = [cLIPbandsCorR1samplesem(1) cLIPbandsIncR1samplesem(1); ...
    cLIPbandsCorR1samplesem(2) cLIPbandsIncR1samplesem(2); ...
    cLIPbandsCorR1samplesem(3) cLIPbandsIncR1samplesem(3); ...
    cLIPbandsCorR1samplesem(4) cLIPbandsIncR1samplesem(4); ...
    cLIPbandsCorR1samplesem(5) cLIPbandsIncR1samplesem(5); ...
    cLIPbandsCorR1samplesem(6) cLIPbandsIncR1samplesem(6); ...
    cLIPbandsCorR1samplesem(7) cLIPbandsIncR1samplesem(7); ...
    cLIPbandsCorR1samplesem(8) cLIPbandsIncR1samplesem(8)];
%db 10.*(log10(y))
figure
hold on
bar(y)
errorbar(y,sem,'LineStyle','none')
% change the xticklabel with characters
newXticklabel = {'?','?','?','?','Low ?','High ?','High ?','High ?'};
set(gca,'XtickLabel',newXticklabel,'FontSize',25);
% xAx = get(gca,'XAxis');
% set(xAx,'FontSize',25);
% yAx = get(gca,'YAxis');
% set(yAx,'FontSize',25);
title('Clark LIP Delay-period');
xlabel('Frequency Bands');
ylabel('Baseline-normalized Power (dB)');
legend('Correct Trials','Incorrect Trials');

%conf intervals
ci = 0.9999 ;
alpha = 1 - ci;
T_multiplier = tinv(1-alpha/2, n-1);
cLIPbandsCorR1sampleci95 = T_multiplier.*std(cLIPbandsCorR1samples)./sqrt(n);
cLIPbandsIncR1sampleci95 = T_multiplier.*std(cLIPbandsIncR1samples)./sqrt(n);
ci95 = [cLIPbandsCorR1sampleci95(1) cLIPbandsIncR1sampleci95(1); ...
    cLIPbandsCorR1sampleci95(2) cLIPbandsIncR1sampleci95(2); ...
    cLIPbandsCorR1sampleci95(3) cLIPbandsIncR1sampleci95(3); ...
    cLIPbandsCorR1sampleci95(4) cLIPbandsIncR1sampleci95(4); ...
    cLIPbandsCorR1sampleci95(5) cLIPbandsIncR1sampleci95(5); ...
    cLIPbandsCorR1sampleci95(6) cLIPbandsIncR1sampleci95(6); ...
    cLIPbandsCorR1sampleci95(7) cLIPbandsIncR1sampleci95(7); ...
    cLIPbandsCorR1sampleci95(8) cLIPbandsIncR1sampleci95(8)];

figure
b = barwitherr(ci95(1,:),y(1,:));
ylim([3.9e8,5.5e8])
newXticklabel = {'Correct','Incorrect'};
set(gca,'XtickLabel',newXticklabel,'FontSize',25);
% xAx = get(gca,'XAxis');
% set(xAx,'FontSize',25);
% yAx = get(gca,'YAxis');
% set(yAx,'FontSize',25);
title('Clark LIP Delay-period ? Band 99.99% CI');
%xlabel('Frequency Bands');
ylabel('Power');
legend('Correct Trials','Incorrect Trials');