%build out matrix for n-way ANOVA, n number of group. here it's 2: region
%and response
regionIdx = cell(16000,1);
respIdx = cell(16000,1);
idx = 1; %index counter

cBetaSamples(1:1000,1) = c8BbandsCorR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'8B'};
    respIdx(idx) = {'Correct'};
    idx = idx + 1;
end
cBetaSamples(1001:2000,1) = c9LbandsCorR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'9L'};
    respIdx(idx) = {'Correct'};
    idx = idx + 1;
end
cBetaSamples(2001:3000,1) = cdPFCbandsCorR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'dPFC'};
    respIdx(idx) = {'Correct'};
    idx = idx + 1;
end
cBetaSamples(3001:4000,1) = cvPFCbandsCorR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'vPFC'};
    respIdx(idx) = {'Correct'};
    idx = idx + 1;
end
cBetaSamples(4001:5000,1) = cLIPbandsCorR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'LIP'};
    respIdx(idx) = {'Correct'};
    idx = idx + 1;
end
cBetaSamples(5001:6000,1) = cMIPbandsCorR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'MIP'};
    respIdx(idx) = {'Correct'};
    idx = idx + 1;
end
cBetaSamples(6001:7000,1) = cPECbandsCorR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'PEC'};
    respIdx(idx) = {'Correct'};
    idx = idx + 1;
end
cBetaSamples(7001:8000,1) = cPGbandsCorR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'PG'};
    respIdx(idx) = {'Correct'};
    idx = idx + 1;
end
cBetaSamples(8001:9000,1) = c8BbandsIncR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'8B'};
    respIdx(idx) = {'Incorrect'};
    idx = idx + 1;
end
cBetaSamples(9001:10000,1) = c9LbandsIncR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'9L'};
    respIdx(idx) = {'Incorrect'};
    idx = idx + 1;
end
cBetaSamples(10001:11000,1) = cdPFCbandsIncR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'dPFC'};
    respIdx(idx) = {'Incorrect'};
    idx = idx + 1;
end
cBetaSamples(11001:12000,1) = cvPFCbandsIncR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'vPFC'};
    respIdx(idx) = {'Incorrect'};
    idx = idx + 1;
end
cBetaSamples(12001:13000,1) = cLIPbandsIncR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'LIP'};
    respIdx(idx) = {'Incorrect'};
    idx = idx + 1;
end
cBetaSamples(13001:14000,1) = cMIPbandsIncR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'MIP'};
    respIdx(idx) = {'Incorrect'};
    idx = idx + 1;
end
cBetaSamples(14001:15000,1) = cPECbandsIncR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'PEC'};
    respIdx(idx) = {'Incorrect'};
    idx = idx + 1;
end
cBetaSamples(15001:16000,1) = cPGbandsIncR1samples(:,4);
for i=1:1000
    regionIdx(idx) = {'PG'};
    respIdx(idx) = {'Incorrect'};
    idx = idx + 1;
end