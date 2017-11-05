clear
%load correct and incorrect for each region, replacing the region as you go
%seed the random number generator for reproducibility
rng(42);
%randomly sample with replacement and take the mean across trials for all
%bands. do this 1,000 times to build sample distribution
c9LbandsCorR1samples = zeros(100,8);
c9LbandsIncR1samples = zeros(100,8);
for i=1:100
    c9LbandsidxCorR1 = randi(size(c9LbandsCorR1,1),size(c9LbandsIncR1,1),1);
    c9LbandsidxIncR1 = randi(size(c9LbandsIncR1,1),size(c9LbandsIncR1,1),1);
    c9LbandsCorR1samples(i,:) = mean(c9LbandsCorR1(c9LbandsidxCorR1,:),1);
    c9LbandsIncR1samples(i,:) = mean(c9LbandsIncR1(c9LbandsidxIncR1,:),1);
end

%sanity check for normal distribution
c9Lhgamma3CorR1 = c9LbandsCorR1samples(:,6);
c9Lhgamma3IncR1 = c9LbandsIncR1samples(:,6);
figure
subplot(1,2,1)
histogram(c9Lhgamma3CorR1);
subplot(1,2,2)
histogram(c9Lhgamma3IncR1);