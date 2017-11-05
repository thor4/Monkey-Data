%determine statistical significance using nonparametric permutation null
%hypothesis distribution created earlier
%
%load null hypothesis distribution along with the normalized power band 
%sample distributions from same area

%run ttest and get t statistic from stats struct
[h,p,ci,stats] = ttest(c8BbandsCorR1samples,c8BbandsIncR1samples);
%standardize t statistic to z score
z = (stats.tstat-mean(c8BbandsR1nullHypoDistrib,1))./std(c8BbandsR1nullHypoDistrib,1);
%transform z score to p value by evaluating its position on a Gaussian
%probability density
p = normcdf(z); %this isn't right, work on this