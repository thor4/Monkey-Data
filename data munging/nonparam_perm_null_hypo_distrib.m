%creating a nonparametric permutation testing null hypothesis distribution
%load distribution samples from relevant area
clear
rng(42); %seed rng for reproducibility
c8BbandsR1nullHypoDistrib = zeros(1000,8); %initialize null hypothesis distribution
for i=1:1000
    %combine and shuffle twice
    mixitup = vertcat(c8BbandsCorR1samples,c8BbandsIncR1samples);
    mixitup2 =  mixitup(randperm(end),:); %once
    mixitup3 =  mixitup(randperm(end),:); %twice
    %null hypothesis says it shouldn't matter which were correct and which were
    %incorrect so arbitrarily split in two
    [~,~,~,stats] = ttest(mixitup3(1:1000,:),mixitup3(1001:2000,:));
    c8BbandsR1nullHypoDistrib(i,:) = stats.tstat;
end

%sanity check for normal distribution
c8Bhgamma1R1nullHypoDistrib = c8BbandsR1nullHypoDistrib(:,6);
c8BdeltaR1nullHypoDistrib = c8BbandsR1nullHypoDistrib(:,1);
figure
subplot(1,2,1)
histogram(c8Bhgamma1R1nullHypoDistrib);
title('High Gamma 1');
subplot(1,2,2)
histogram(c8BdeltaR1nullHypoDistrib);
title('Delta');