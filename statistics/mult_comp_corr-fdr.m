clear
%load uncorrected p values for all regions within one monkey
%false discovery rate correction
pvals = zeros(8,8); %initialize
%combine p values from all regions (rows) by bands (columns)
% rows={'8B','9L','dPFC','vPFC','LIP','MIP','PEC','PG'}; %clark
rows={'6DR','8AD','8B','dPFC','LIP','PE','PEC','PG'}; %clark
columns={'delta','theta','alpha','beta','low gamma','high gamma1','high gamma2','high gamma3'};
pvals(1,:) = b6DRbandsR1p;
pvals(2,:) = b8ADbandsR1p;
pvals(3,:) = b8BbandsR1p;
pvals(4,:) = bdPFCbandsR1p;
pvals(5,:) = bLIPbandsR1p;
pvals(6,:) = bPEbandsR1p;
pvals(7,:) = bPECbandsR1p;
pvals(8,:) = bPGbandsR1p;
%run FDR with report 'yes' and .05 alpha 
%return new critical p (new alpha) and p-values adjusted for mult comp
[h crit_p adj_p]=fdr_bh(pvals,.0001,'pdep','yes');

