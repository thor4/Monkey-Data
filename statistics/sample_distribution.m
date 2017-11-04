clear
%load correct and incorrect for each region, replacing the region as you go
%seed the random number generator for reproducibility
rng(42);
%randomly sample with replacement and take the mean. do this 1,000 times to
%build sample distribution
cvPFCbandsCorR1samples = zeros(1000,8);
cvPFCbandsIncR1samples = zeros(1000,8);
for i=1:1000
    cvPFCbandsidxCorR1 = randi(size(cvPFCbandsCorR1,1),size(cvPFCbandsIncR1,1),1);
    cvPFCbandsidxIncR1 = randi(size(cvPFCbandsIncR1,1),size(cvPFCbandsIncR1,1),1);
    cvPFCbandsCorR1samples(i,:) = mean(cvPFCbandsCorR1(cvPFCbandsidxCorR1,:),1);
    cvPFCbandsIncR1samples(i,:) = mean(cvPFCbandsIncR1(cvPFCbandsidxIncR1,:),1);
end

%sanity check for normal distribution
cvPFChgamma3CorR1 = cvPFCbandsCorR1samples(:,6);
cvPFChgamma3IncR1 = cvPFCbandsIncR1samples(:,6);
figure
subplot(1,2,1)
histogram(cvPFChgamma3CorR1);
subplot(1,2,2)
histogram(cvPFChgamma3IncR1);