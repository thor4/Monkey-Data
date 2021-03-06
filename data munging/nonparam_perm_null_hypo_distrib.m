%creating a nonparametric permutation testing null hypothesis distribution
%load distribution samples from relevant area
clear
rng(42); %seed rng for reproducibility
bPGbandsR1nullHypoDistrib = zeros(1000,8); %initialize null hypothesis distribution
%combine
mixitup = vertcat(bPGbandsCorR1sampled,bPGbandsIncR1);
for i=1:1000
    %shuffle
    mixitup =  mixitup(randperm(end),:); %once
    %null hypothesis says it shouldn't matter which were correct and which were
    %incorrect so arbitrarily split in two
    [~,~,~,stats] = ttest(mixitup(1:100,:),mixitup(101:200,:));
    bPGbandsR1nullHypoDistrib(i,:) = stats.tstat;
end

%sanity check for normal distribution
bPGhgamma1R1nullHypoDistrib = bPGbandsR1nullHypoDistrib(:,6);
bPGdeltaR1nullHypoDistrib = bPGbandsR1nullHypoDistrib(:,1);
figure
subplot(1,2,1)
histogram(bPGhgamma1R1nullHypoDistrib);
title('High Gamma 1');
subplot(1,2,2)
histogram(bPGdeltaR1nullHypoDistrib);
title('Delta');