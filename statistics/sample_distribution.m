clear
%load correct and incorrect for each region, replacing the region as you go
%seed the random number generator for reproducibility
rng(42);
%randomly sample with replacement and take the mean across trials for all
%bands. do this 1,000 times to build sample distribution
c8BbandsCorR1samples = zeros(1000,8);
c8BbandsIncR1samples = zeros(1000,8);
for i=1:1000
    c8BbandsidxCorR1 = randi(size(c8BbandsCorR1,1),size(c8BbandsIncR1,1),1);
    c8BbandsidxIncR1 = randi(size(c8BbandsIncR1,1),size(c8BbandsIncR1,1),1);
    c8BbandsCorR1samples(i,:) = mean(c8BbandsCorR1(c8BbandsidxCorR1,:),1);
    c8BbandsIncR1samples(i,:) = mean(c8BbandsIncR1(c8BbandsidxIncR1,:),1);
end

%sanity check for normal distribution
c8Bhgamma3CorR1 = c8BbandsCorR1samples(:,6);
c8Bhgamma3IncR1 = c8BbandsIncR1samples(:,6);
figure
subplot(1,2,1)
histogram(c8Bhgamma3CorR1);
subplot(1,2,2)
histogram(c8Bhgamma3IncR1);