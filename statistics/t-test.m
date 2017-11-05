%determine statistical significance using nonparametric permutation null
%hypothesis distribution created earlier
%
%load null hypothesis distribution along with the normalized power band 
%sample distributions from same area

%run ttest and get t statistic from stats struct
[h,p,ci,stats] = ttest(bPGbandsCorR1sampled,bPGbandsIncR1);
%standardize t statistic to z score
bPGbandsR1z = (stats.tstat-mean(bPGbandsR1nullHypoDistrib,1))./std(bPGbandsR1nullHypoDistrib,1);
%transform z score to p value by evaluating its position on a Gaussian
%probability density
bPGbandsR1p = zeros(1,8);
format long
for i=1:8
    if bPGbandsR1z(i) > 0
        bPGbandsR1p(i) = 1-normcdf(bPGbandsR1z(i),mean(bPGbandsR1nullHypoDistrib(:,i),1),std(bPGbandsR1nullHypoDistrib(:,i),1));
    else
        bPGbandsR1p(i) = normcdf(bPGbandsR1z(i),mean(bPGbandsR1nullHypoDistrib(:,i),1),std(bPGbandsR1nullHypoDistrib(:,i),1));
    end
end