load('raw-0_20_filtered-base_norm.mat');
cNorm = cNorm'; iNorm = iNorm';
%data should be trials x samples
cNorm = cNorm'; iNorm = iNorm';

%ftest to test whether two time points in a time series have significantly
%different variances. if so, don't assume equal variance in the ttest
[hf,pf,cif,statsf] = vartest2(cNorm,iNorm,'Alpha',0.01);
find(hf==0) %if empty, shows all were rejected with h=1
%all were rejected, which means we can't assume equal variance in the ttest

%ttest to test whether two time points in a time series are significantly 
%different, 

[h,p,ci,stats] = ttest2(cNorm,iNorm,'Vartype','unequal','Alpha',0.001);
find(h==0) 
%shows only 21 time points were not significant: 1 & [519-538]