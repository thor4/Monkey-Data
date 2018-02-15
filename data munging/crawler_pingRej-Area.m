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
%correct, monkey 2
m2CorrectRule16DRThrowAway = any(m2CorrectRule16DR > 350 | m2CorrectRule16DR < -350, 2);
numArtifactM2CorrectRule16DR = size(m2CorrectRule16DR,1);
m2CorrectRule18ADThrowAway = any(m2CorrectRule18AD > 350 | m2CorrectRule18AD < -350, 2);
numArtifactM2CorrectRule18AD = size(m2CorrectRule18AD,1);
m2CorrectRule18BThrowAway = any(m2CorrectRule18B > 350 | m2CorrectRule18B < -350, 2);
numArtifactM2CorrectRule18B = size(m2CorrectRule18B,1);
m2CorrectRule1dPFCThrowAway = any(m2CorrectRule1dPFC > 350 | m2CorrectRule1dPFC < -350, 2);
numArtifactM2CorrectRule1dPFC = size(m2CorrectRule1dPFC,1);
m2CorrectRule1LIPThrowAway = any(m2CorrectRule1LIP > 350 | m2CorrectRule1LIP < -350, 2);
numArtifactM2CorrectRule1LIP = size(m2CorrectRule1LIP,1);
m2CorrectRule1PEThrowAway = any(m2CorrectRule1PE > 350 | m2CorrectRule1PE < -350, 2);
numArtifactM2CorrectRule1PE = size(m2CorrectRule1PE,1);
m2CorrectRule1PECThrowAway = any(m2CorrectRule1PEC > 350 | m2CorrectRule1PEC < -350, 2);
numArtifactM2CorrectRule1PEC = size(m2CorrectRule1PEC,1);
m2CorrectRule1PGThrowAway = any(m2CorrectRule1PG > 350 | m2CorrectRule1PG < -350, 2);
numArtifactM2CorrectRule1PG = size(m2CorrectRule1PG,1);
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
%incorrect, monkey 2
m2IncorrectRule16DRThrowAway = any(m2IncorrectRule16DR > 350 | m2IncorrectRule16DR < -350, 2);
numArtifactM2IncorrectRule16DR = size(m2IncorrectRule16DR,1);
m2IncorrectRule18ADThrowAway = any(m2IncorrectRule18AD > 350 | m2IncorrectRule18AD < -350, 2);
numArtifactM2IncorrectRule18AD = size(m2IncorrectRule18AD,1);
m2IncorrectRule18BThrowAway = any(m2IncorrectRule18B > 350 | m2IncorrectRule18B < -350, 2);
numArtifactM2IncorrectRule18B = size(m2IncorrectRule18B,1);
m2IncorrectRule1dPFCThrowAway = any(m2IncorrectRule1dPFC > 350 | m2IncorrectRule1dPFC < -350, 2);
numArtifactM2IncorrectRule1dPFC = size(m2IncorrectRule1dPFC,1);
m2IncorrectRule1LIPThrowAway = any(m2IncorrectRule1LIP > 350 | m2IncorrectRule1LIP < -350, 2);
numArtifactM2IncorrectRule1LIP = size(m2IncorrectRule1LIP,1);
m2IncorrectRule1PEThrowAway = any(m2IncorrectRule1PE > 350 | m2IncorrectRule1PE < -350, 2);
numArtifactM2IncorrectRule1PE = size(m2IncorrectRule1PE,1);
m2IncorrectRule1PECThrowAway = any(m2IncorrectRule1PEC > 350 | m2IncorrectRule1PEC < -350, 2);
numArtifactM2IncorrectRule1PEC = size(m2IncorrectRule1PEC,1);
m2IncorrectRule1PGThrowAway = any(m2IncorrectRule1PG > 350 | m2IncorrectRule1PG < -350, 2);
numArtifactM2IncorrectRule1PG = size(m2IncorrectRule1PG,1);

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
%correct, monkey 2
m2CorrectRule16DR(m2CorrectRule16DRThrowAway,:) = [];
numNoArtifactM2CorrectRule16DR = size(m2CorrectRule16DR,1);
m2CorrectRule18AD(m2CorrectRule18ADThrowAway,:) = [];
numNoArtifactM2CorrectRule18AD = size(m2CorrectRule18AD,1);
m2CorrectRule18B(m2CorrectRule18BThrowAway,:) = [];
numNoArtifactM2CorrectRule18B = size(m2CorrectRule18B,1);
m2CorrectRule1dPFC(m2CorrectRule1dPFCThrowAway,:) = [];
numNoArtifactM2CorrectRule1dPFC = size(m2CorrectRule1dPFC,1);
m2CorrectRule1LIP(m2CorrectRule1LIPThrowAway,:) = [];
numNoArtifactM2CorrectRule1LIP = size(m2CorrectRule1LIP,1);
m2CorrectRule1PE(m2CorrectRule1PEThrowAway,:) = [];
numNoArtifactM2CorrectRule1PE = size(m2CorrectRule1PE,1);
m2CorrectRule1PEC(m2CorrectRule1PECThrowAway,:) = [];
numNoArtifactM2CorrectRule1PEC = size(m2CorrectRule1PEC,1);
m2CorrectRule1PG(m2CorrectRule1PGThrowAway,:) = [];
numNoArtifactM2CorrectRule1PG = size(m2CorrectRule1PG,1);
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
%incorrect, monkey 2
m2IncorrectRule16DR(m2IncorrectRule16DRThrowAway,:) = [];
numNoArtifactM2IncorrectRule16DR = size(m2IncorrectRule16DR,1);
m2IncorrectRule18AD(m2IncorrectRule18ADThrowAway,:) = [];
numNoArtifactM2IncorrectRule18AD = size(m2IncorrectRule18AD,1);
m2IncorrectRule18B(m2IncorrectRule18BThrowAway,:) = [];
numNoArtifactM2IncorrectRule18B = size(m2IncorrectRule18B,1);
m2IncorrectRule1dPFC(m2IncorrectRule1dPFCThrowAway,:) = [];
numNoArtifactM2IncorrectRule1dPFC = size(m2IncorrectRule1dPFC,1);
m2IncorrectRule1LIP(m2IncorrectRule1LIPThrowAway,:) = [];
numNoArtifactM2IncorrectRule1LIP = size(m2IncorrectRule1LIP,1);
m2IncorrectRule1PE(m2IncorrectRule1PEThrowAway,:) = [];
numNoArtifactM2IncorrectRule1PE = size(m2IncorrectRule1PE,1);
m2IncorrectRule1PEC(m2IncorrectRule1PECThrowAway,:) = [];
numNoArtifactM2IncorrectRule1PEC = size(m2IncorrectRule1PEC,1);
m2IncorrectRule1PG(m2IncorrectRule1PGThrowAway,:) = [];
numNoArtifactM2IncorrectRule1PG = size(m2IncorrectRule1PG,1);

%count total number of artifacts removed
%correct, monkey 1
artifactM1CorrectRule18B = numArtifactM1CorrectRule18B - numNoArtifactM1CorrectRule18B;
artifactM1CorrectRule19L = numArtifactM1CorrectRule19L - numNoArtifactM1CorrectRule19L;
artifactM1CorrectRule1dPFC = numArtifactM1CorrectRule1dPFC - numNoArtifactM1CorrectRule1dPFC;
artifactM1CorrectRule1vPFC = numArtifactM1CorrectRule1vPFC - numNoArtifactM1CorrectRule1vPFC;
artifactM1CorrectRule1LIP = numArtifactM1CorrectRule1LIP - numNoArtifactM1CorrectRule1LIP;
artifactM1CorrectRule1MIP = numArtifactM1CorrectRule1MIP - numNoArtifactM1CorrectRule1MIP;
artifactM1CorrectRule1PEC = numArtifactM1CorrectRule1PEC - numNoArtifactM1CorrectRule1PEC;
artifactM1CorrectRule1PG = numArtifactM1CorrectRule1PG - numNoArtifactM1CorrectRule1PG;
%correct, monkey 2
artifactM2CorrectRule16DR = numArtifactM2CorrectRule16DR - numNoArtifactM2CorrectRule16DR;
artifactM2CorrectRule18AD = numArtifactM2CorrectRule18AD - numNoArtifactM2CorrectRule18AD;
artifactM2CorrectRule18B = numArtifactM2CorrectRule18B - numNoArtifactM2CorrectRule18B;
artifactM2CorrectRule1dPFC = numArtifactM2CorrectRule1dPFC - numNoArtifactM2CorrectRule1dPFC;
artifactM2CorrectRule1LIP = numArtifactM2CorrectRule1LIP - numNoArtifactM2CorrectRule1LIP;
artifactM2CorrectRule1PE = numArtifactM2CorrectRule1PE - numNoArtifactM2CorrectRule1PE;
artifactM2CorrectRule1PEC = numArtifactM2CorrectRule1PEC - numNoArtifactM2CorrectRule1PEC;
artifactM2CorrectRule1PG = numArtifactM2CorrectRule1PG - numNoArtifactM2CorrectRule1PG;
%incorrect, monkey 1
artifactM1IncorrectRule18B = numArtifactM1IncorrectRule18B - numNoArtifactM1IncorrectRule18B;
artifactM1IncorrectRule19L = numArtifactM1IncorrectRule19L - numNoArtifactM1IncorrectRule19L;
artifactM1IncorrectRule1dPFC = numArtifactM1IncorrectRule1dPFC - numNoArtifactM1IncorrectRule1dPFC;
artifactM1IncorrectRule1vPFC = numArtifactM1IncorrectRule1vPFC - numNoArtifactM1IncorrectRule1vPFC;
artifactM1IncorrectRule1LIP = numArtifactM1IncorrectRule1LIP - numNoArtifactM1IncorrectRule1LIP;
artifactM1IncorrectRule1MIP = numArtifactM1IncorrectRule1MIP - numNoArtifactM1IncorrectRule1MIP;
artifactM1IncorrectRule1PEC = numArtifactM1IncorrectRule1PEC - numNoArtifactM1IncorrectRule1PEC;
artifactM1IncorrectRule1PG = numArtifactM1IncorrectRule1PG - numNoArtifactM1IncorrectRule1PG;
%incorrect, monkey 2
artifactM2IncorrectRule16DR = numArtifactM2IncorrectRule16DR - numNoArtifactM2IncorrectRule16DR;
artifactM2IncorrectRule18AD = numArtifactM2IncorrectRule18AD - numNoArtifactM2IncorrectRule18AD;
artifactM2IncorrectRule18B = numArtifactM2IncorrectRule18B - numNoArtifactM2IncorrectRule18B;
artifactM2IncorrectRule1dPFC = numArtifactM2IncorrectRule1dPFC - numNoArtifactM2IncorrectRule1dPFC;
artifactM2IncorrectRule1LIP = numArtifactM2IncorrectRule1LIP - numNoArtifactM2IncorrectRule1LIP;
artifactM2IncorrectRule1PE = numArtifactM2IncorrectRule1PE - numNoArtifactM2IncorrectRule1PE;
artifactM2IncorrectRule1PEC = numArtifactM2IncorrectRule1PEC - numNoArtifactM2IncorrectRule1PEC;
artifactM2IncorrectRule1PG = numArtifactM2IncorrectRule1PG - numNoArtifactM2IncorrectRule1PG;