%prepare matrix for 2-way ANOVA (independent measures)
%build observation matrix: BehResp in rows and Brodmann area in columns
%Rows: 1000xCorrect, 1000xIncorrect
%Columns: 1=8B, 2=9L, 3=dPFC, 4=vPFC, 5=LIP, 6=MIP, 7=PEC, 8=PG
cBetaSamples = zeros(2000,8);
cBetaSamples(1:1000,1) = c8BbandsCorR1samples(:,4);
cBetaSamples(1:1000,2) = c9LbandsCorR1samples(:,4);
cBetaSamples(1:1000,3) = cdPFCbandsCorR1samples(:,4);
cBetaSamples(1:1000,4) = cvPFCbandsCorR1samples(:,4);
cBetaSamples(1:1000,5) = cLIPbandsCorR1samples(:,4);
cBetaSamples(1:1000,6) = cMIPbandsCorR1samples(:,4);
cBetaSamples(1:1000,7) = cPECbandsCorR1samples(:,4);
cBetaSamples(1:1000,8) = cPGbandsCorR1samples(:,4);
cBetaSamples(1001:2000,1) = c8BbandsIncR1samples(:,4);
cBetaSamples(1001:2000,2) = c9LbandsIncR1samples(:,4);
cBetaSamples(1001:2000,3) = cdPFCbandsIncR1samples(:,4);
cBetaSamples(1001:2000,4) = cvPFCbandsIncR1samples(:,4);
cBetaSamples(1001:2000,5) = cLIPbandsIncR1samples(:,4);
cBetaSamples(1001:2000,6) = cMIPbandsIncR1samples(:,4);
cBetaSamples(1001:2000,7) = cPECbandsIncR1samples(:,4);
cBetaSamples(1001:2000,8) = cPGbandsIncR1samples(:,4);

%group 1
for i=1:8
    cBetaSamples(1:1000,1) = c8BbandsCorR1samples(:,4);
    cBetaSamples(1:1000,2) = c9LbandsCorR1samples(:,4);
cBetaSamples(1:1000,3) = cdPFCbandsCorR1samples(:,4);
cBetaSamples(1:1000,4) = cvPFCbandsCorR1samples(:,4);
cBetaSamples(1:1000,5) = cLIPbandsCorR1samples(:,4);
cBetaSamples(1:1000,6) = cMIPbandsCorR1samples(:,4);
cBetaSamples(1:1000,7) = cPECbandsCorR1samples(:,4);
cBetaSamples(1:1000,8) = cPGbandsCorR1samples(:,4);
cBetaSamples(1001:2000,1) = c8BbandsIncR1samples(:,4);
cBetaSamples(1001:2000,2) = c9LbandsIncR1samples(:,4);
cBetaSamples(1001:2000,3) = cdPFCbandsIncR1samples(:,4);
cBetaSamples(1001:2000,4) = cvPFCbandsIncR1samples(:,4);
cBetaSamples(1001:2000,5) = cLIPbandsIncR1samples(:,4);
cBetaSamples(1001:2000,6) = cMIPbandsIncR1samples(:,4);
cBetaSamples(1001:2000,7) = cPECbandsIncR1samples(:,4);
cBetaSamples(1001:2000,8) = cPGbandsIncR1samples(:,4);
resp = {'hi';'hi';'lo'};