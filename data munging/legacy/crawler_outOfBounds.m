load('bettyGoodStableRule1-split_by_BehResp_and_Region.mat')
load('clarkGoodStableRule1-split_by_BehResp_and_Region.mat')

%specify conditions and choose only those row indices that have at least 
%one column instance satisfying conditions
%correct
correctFrontalM1ThrowAway = any(correctFrontalClark > 750 | correctFrontalClark < -750, 2);
numArtifactCorrectFrontalM1 = size(correctFrontalClark,1);
correctParietalM1ThrowAway = any(correctParietalClark > 750 | correctParietalClark < -750, 2);
numArtifactCorrectParietalM1 = size(correctParietalClark,1);
correctFrontalM2ThrowAway = any(correctFrontalBetty > 350 | correctFrontalBetty < -350, 2);
numArtifactCorrectFrontalM2 = size(correctFrontalBetty,1);
correctParietalM2ThrowAway = any(correctParietalBetty > 350 | correctParietalBetty < -350, 2);
numArtifactCorrectParietalM2 = size(correctParietalBetty,1);
%incorrect
incorrectFrontalM1ThrowAway = any(incorrectFrontalClark > 750 | incorrectFrontalClark < -750, 2);
numArtifactIncorrectFrontalM1 = size(incorrectFrontalClark,1);
incorrectParietalM1ThrowAway = any(incorrectParietalClark > 750 | incorrectParietalClark < -750, 2);
numArtifactIncorrectParietalM1 = size(incorrectParietalClark,1);
incorrectFrontalM2ThrowAway = any(incorrectFrontalBetty > 350 | incorrectFrontalBetty < -350, 2);
numArtifactIncorrectFrontalM2 = size(incorrectFrontalBetty,1);
incorrectParietalM2ThrowAway = any(incorrectParietalBetty > 350 | incorrectParietalBetty < -350, 2);
numArtifactIncorrectParietalM2 = size(incorrectParietalBetty,1);

%remove entire rows where at least one instance out of bounds occurs
%correct
correctFrontalClark(correctFrontalM1ThrowAway,:) = [];
numNoArtifactCorrectFrontalM1 = size(correctFrontalClark,1);
correctParietalClark(correctParietalM1ThrowAway,:) = [];
numNoArtifactCorrectParietalM1 = size(correctParietalClark,1);
correctFrontalBetty(correctFrontalM2ThrowAway,:) = [];
numNoArtifactCorrectFrontalM2 = size(correctFrontalBetty,1);
correctParietalBetty(correctParietalM2ThrowAway,:) = [];
numNoArtifactCorrectParietalM2 = size(correctParietalBetty,1);
%incorrect
incorrectFrontalClark(incorrectFrontalM1ThrowAway,:) = [];
numNoArtifactIncorrectFrontalM1 = size(incorrectFrontalClark,1);
incorrectParietalClark(incorrectParietalM1ThrowAway,:) = [];
numNoArtifactIncorrectParietalM1 = size(incorrectParietalClark,1);
incorrectFrontalBetty(incorrectFrontalM2ThrowAway,:) = [];
numNoArtifactIncorrectFrontalM2 = size(incorrectFrontalBetty,1);
incorrectParietalBetty(incorrectParietalM2ThrowAway,:) = [];
numNoArtifactIncorrectParietalM2 = size(incorrectParietalBetty,1);

%count total number of artifacts removed
%correct
artifactTrialsCorrectFrontalM1 = numArtifactCorrectFrontalM1 - numNoArtifactCorrectFrontalM1;
artifactTrialsCorrectParietalM1 = numArtifactCorrectParietalM1 - numNoArtifactCorrectParietalM1;
artifactTrialsCorrectFrontalM2 = numArtifactCorrectFrontalM2 - numNoArtifactCorrectFrontalM2;
artifactTrialsCorrectParietalM2 = numArtifactCorrectParietalM2 - numNoArtifactCorrectParietalM2;
%incorrect
artifactTrialsIncorrectFrontalM1 = numArtifactIncorrectFrontalM1 - numNoArtifactIncorrectFrontalM1;
artifactTrialsIncorrectParietalM1 = numArtifactIncorrectParietalM1 - numNoArtifactIncorrectParietalM1;
artifactTrialsIncorrectFrontalM2 = numArtifactIncorrectFrontalM2 - numNoArtifactIncorrectFrontalM2;
artifactTrialsIncorrectParietalM2 = numArtifactIncorrectParietalM2 - numNoArtifactIncorrectParietalM2;