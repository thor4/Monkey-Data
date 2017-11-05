%equalize correct to incorrect
clear
%load correct and incorrect for each region, replacing the region as you go
%seed the random number generator for reproducibility
rng(42);
%randomly sample with replacement and take the mean across trials for all
%bands. 
bPGbandsidxCorR1 = randi(size(bPGbandsCorR1,1),size(bPGbandsIncR1,1),1);
bPGbandsCorR1sampled = bPGbandsCorR1(bPGbandsidxCorR1,:);