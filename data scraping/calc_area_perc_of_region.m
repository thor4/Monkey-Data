
liarea9L = length(Incorrect.area9L);
liarea8B = length(Incorrect.area8B);
liarea6DR = length(Incorrect.area6DR);
liarea8AD = length(Incorrect.area8AD);
liareavPFC = length(Incorrect.areavPFC);
liareadPFC = length(Incorrect.areadPFC);
liareaLIP = length(Incorrect.areaLIP);
liareaMIP = length(Incorrect.areaMIP);
liareaPE = length(Incorrect.areaPE);
liareaPG = length(Incorrect.areaPG);
liareaPEC = length(Incorrect.areaPEC);
lif = liarea9L + liarea8B + liarea6DR + liarea8AD + liareavPFC + liareadPFC;
lip = liareaLIP + liareaMIP + liareaPE + liareaPG + liareaPEC;

lcarea9L = length(Correct.area9L);
lcarea8B = length(Correct.area8B);
lcarea6DR = length(Correct.area6DR);
lcarea8AD = length(Correct.area8AD);
lcareavPFC = length(Correct.areavPFC);
lcareadPFC = length(Correct.areadPFC);
lcareaLIP = length(Correct.areaLIP);
lcareaMIP = length(Correct.areaMIP);
lcareaPE = length(Correct.areaPE);
lcareaPG = length(Correct.areaPG);
lcareaPEC = length(Correct.areaPEC);
lcf = lcarea9L + lcarea8B + lcarea6DR + lcarea8AD + lcareavPFC + lcareadPFC;
lcp = lcareaLIP + lcareaMIP + lcareaPE + lcareaPG + lcareaPEC;

% calculate percentages for each area of incorrect & correct
% frontal
fi6DR = liarea6DR / lif; fc6DR = lcarea6DR / lcf;
fi8AD = liarea8AD / lif; fc8AD = lcarea8AD / lcf;
fi8B = liarea8B / lif; fc8B = lcarea8B / lcf;
fi9L = liarea9L / lif; fc9L = lcarea9L / lcf;
fidPFC = liareadPFC / lif; fcdPFC = lcareadPFC / lcf;
fivPFC = liareavPFC / lif; fcvPFC = lcareavPFC / lcf;
% parietal
piLIP = liareaLIP / lip; pcLIP = lcareaLIP / lcp;
piMIP = liareaMIP / lip; pcMIP = lcareaMIP / lcp;
piPE = liareaPE / lip; pcPE = lcareaPE / lcp;
piPEC = liareaPEC / lip; pcPEC = lcareaPEC / lcp;
piPG = liareaPG / lip; pcPG = lcareaPG / lcp;