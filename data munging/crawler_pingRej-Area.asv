load('m1GoodStableRule1-split_by_BehResp_and_Area.mat')
load('m2GoodStableRule1-split_by_BehResp_and_Area.mat')

%specify conditions and choose only those row indices that have at least 
%one column instance satisfying conditions
%correct, monkey 1
m1CorrectRule18BThrowAway = any(m1CorrectRule18B > 750 | m1CorrectRule18B < -750, 2);
numArtifactM1CorrectRule18B = size(m1CorrectRule18B,1);
m1CorrectRule19LThrowAway = any(m1CorrectRule19L > 750 | m1CorrectRule19L < -750, 2);
numArtifactM1CorrectRule19L = size(m1CorrectRule19L,1);
m1CorrectRule1dPFCThrowAway = any(m1CorrectRule1dPFC > 750 | m1CorrectRule1dPFC < -750, 2);
numArtifactM1CorrectRule1dPFC = size(m1CorrectRule1dPFC,1);
m1CorrectRule1LIPThrowAway = any(m1CorrectRule1LIP > 750 | m1CorrectRule1LIP < -750, 2);
numArtifactM1CorrectRule1LIP = size(m1CorrectRule1LIP,1);
m1CorrectRule1MIPThrowAway = any(m1CorrectRule1MIP > 750 | m1CorrectRule1MIP < -750, 2);
numArtifactM1CorrectRule1MIP = size(m1CorrectRule1MIP,1);
m1CorrectRule1PECThrowAway = any(m1CorrectRule1PEC > 750 | m1CorrectRule1PEC < -750, 2);
numArtifactM1CorrectRule1PEC = size(m1CorrectRule1PEC,1);
m1CorrectRule1PGThrowAway = any(m1CorrectRule1PG > 750 | m1CorrectRule1PG < -750, 2);
numArtifactM1CorrectRule1PG = size(m1CorrectRule1PG,1);
m1CorrectRule1vPFCThrowAway = any(m1CorrectRule1vPFC > 750 | m1CorrectRule1vPFC < -750, 2);
numArtifactM1CorrectRule1vPFC = size(m1CorrectRule1vPFC,1);

%incorrect, monkey 1
m1IncorrectRule18BThrowAway = any(m1IncorrectRule18B > 750 | m1IncorrectRule18B < -750, 2);
numArtifactM1IncorrectRule18B = size(m1IncorrectRule18B,1);
m1IncorrectRule19LThrowAway = any(m1IncorrectRule19L > 750 | m1IncorrectRule19L < -750, 2);
numArtifactM1IncorrectRule19L = size(m1IncorrectRule19L,1);
m1IncorrectRule1dPFCThrowAway = any(m1IncorrectRule1dPFC > 750 | m1IncorrectRule1dPFC < -750, 2);
numArtifactM1IncorrectRule1dPFC = size(m1IncorrectRule1dPFC,1);
m1IncorrectRule1LIPThrowAway = any(m1IncorrectRule1LIP > 750 | m1IncorrectRule1LIP < -750, 2);
numArtifactM1IncorrectRule1LIP = size(m1IncorrectRule1LIP,1);
m1IncorrectRule1MIPThrowAway = any(m1IncorrectRule1MIP > 750 | m1IncorrectRule1MIP < -750, 2);
numArtifactM1IncorrectRule1MIP = size(m1IncorrectRule1MIP,1);
m1IncorrectRule1PECThrowAway = any(m1IncorrectRule1PEC > 750 | m1IncorrectRule1PEC < -750, 2);
numArtifactM1IncorrectRule1PEC = size(m1IncorrectRule1PEC,1);
m1IncorrectRule1PGThrowAway = any(m1IncorrectRule1PG > 750 | m1IncorrectRule1PG < -750, 2);
numArtifactM1IncorrectRule1PG = size(m1IncorrectRule1PG,1);
m1IncorrectRule1vPFCThrowAway = any(m1IncorrectRule1vPFC > 750 | m1IncorrectRule1vPFC < -750, 2);
numArtifactM1IncorrectRule1vPFC = size(m1IncorrectRule1vPFC,1);

%remove entire rows where at least one instance out of bounds occurs
%correct, monkey 1
m1CorrectRule18B(m1CorrectRule18BThrowAway,:) = [];
numNoArtifactM1CorrectRule18B = size(m1CorrectRule18B,1);
m1CorrectRule19L(m1CorrectRule19LThrowAway,:) = [];
numNoArtifactM1CorrectRule19L = size(m1CorrectRule19L,1);
m1CorrectRule1dPFC(m1CorrectRule1dPFCThrowAway,:) = [];
numNoArtifactM1CorrectRule1dPFC = size(m1CorrectRule1dPFC,1);
m1CorrectRule1vPFC(m1CorrectRule1vPFCThrowAway,:) = [];
numNoArtifactM1CorrectRule1vPFC = size(m1CorrectRule1vPFC,1);
m1CorrectRule1LIP(m1CorrectRule1LIPThrowAway,:) = [];
numNoArtifactM1CorrectRule1LIP = size(m1CorrectRule1LIP,1);
m1CorrectRule1MIP(m1CorrectRule1MIPThrowAway,:) = [];
numNoArtifactM1CorrectRule1MIP = size(m1CorrectRule1MIP,1);
m1CorrectRule1PEC(m1CorrectRule1PECThrowAway,:) = [];
numNoArtifactM1CorrectRule1PEC = size(m1CorrectRule1PEC,1);
m1CorrectRule1PG(m1CorrectRule1PGThrowAway,:) = [];
numNoArtifactM1CorrectRule1PG = size(m1CorrectRule1PG,1);
%incorrect, monkey 1
m1IncorrectRule18B(m1IncorrectRule18BThrowAway,:) = [];
numNoArtifactM1IncorrectRule18B = size(m1IncorrectRule18B,1);
m1IncorrectRule19L(m1IncorrectRule19LThrowAway,:) = [];
numNoArtifactM1IncorrectRule19L = size(m1IncorrectRule19L,1);
m1IncorrectRule1dPFC(m1IncorrectRule1dPFCThrowAway,:) = [];
numNoArtifactM1IncorrectRule1dPFC = size(m1IncorrectRule1dPFC,1);
m1IncorrectRule1vPFC(m1IncorrectRule1vPFCThrowAway,:) = [];
numNoArtifactM1IncorrectRule1vPFC = size(m1IncorrectRule1vPFC,1);
m1IncorrectRule1LIP(m1IncorrectRule1LIPThrowAway,:) = [];
numNoArtifactM1IncorrectRule1LIP = size(m1IncorrectRule1LIP,1);
m1IncorrectRule1MIP(m1IncorrectRule1MIPThrowAway,:) = [];
numNoArtifactM1IncorrectRule1MIP = size(m1IncorrectRule1MIP,1);
m1IncorrectRule1PEC(m1IncorrectRule1PECThrowAway,:) = [];
numNoArtifactM1IncorrectRule1PEC = size(m1IncorrectRule1PEC,1);
m1IncorrectRule1PG(m1IncorrectRule1PGThrowAway,:) = [];
numNoArtifactM1IncorrectRule1PG = size(m1IncorrectRule1PG,1);

%count total number of artifacts removed
%correct
artifactM1CorrectRule18B = numArtifactM1CorrectRule18B - numNoArtifactM1CorrectRule18B;
artifactM1CorrectRule19L = numArtifactM1CorrectRule19L - numNoArtifactM1CorrectRule19L;
artifactM1CorrectRule1dPFC = numArtifactM1CorrectRule1dPFC - numNoArtifactM1CorrectRule1dPFC;
artifactM1CorrectRule1vPFC = numArtifactM1CorrectRule1vPFC - numNoArtifactM1CorrectRule1vPFC;
artifactM1CorrectRule1LIP = numArtifactM1CorrectRule1LIP - numNoArtifactM1CorrectRule1LIP;
artifactM1CorrectRule1MIP = numArtifactM1CorrectRule1MIP - numNoArtifactM1CorrectRule1MIP;
artifactM1CorrectRule1PEC = numArtifactM1CorrectRule1PEC - numNoArtifactM1CorrectRule1PEC;
artifactM1CorrectRule1PG = numArtifactM1CorrectRule1PG - numNoArtifactM1CorrectRule1PG;

%incorrect
artifactM1IncorrectRule18B = numArtifactM1IncorrectRule18B - numNoArtifactM1IncorrectRule18B;
artifactM1IncorrectRule19L = numArtifactM1IncorrectRule19L - numNoArtifactM1IncorrectRule19L;
artifactM1IncorrectRule1dPFC = numArtifactM1IncorrectRule1dPFC - numNoArtifactM1IncorrectRule1dPFC;
artifactM1IncorrectRule1vPFC = numArtifactM1IncorrectRule1vPFC - numNoArtifactM1IncorrectRule1vPFC;
artifactM1IncorrectRule1LIP = numArtifactM1IncorrectRule1LIP - numNoArtifactM1IncorrectRule1LIP;
artifactM1IncorrectRule1MIP = numArtifactM1IncorrectRule1MIP - numNoArtifactM1IncorrectRule1MIP;
artifactM1IncorrectRule1PEC = numArtifactM1IncorrectRule1PEC - numNoArtifactM1IncorrectRule1PEC;
artifactM1IncorrectRule1PG = numArtifactM1IncorrectRule1PG - numNoArtifactM1IncorrectRule1PG;