%initialize descriptive statistics matrices
c1stat = zeros(length(correct1),7);
c2stat = zeros(length(correct2),7);
c3stat = zeros(length(correct3),7);

i1stat = zeros(length(incorrect1),7);
i2stat = zeros(length(incorrect2),7);
i3stat = zeros(length(incorrect3),7);

test = [1, 7, 3; 4, 5, 6];
var(test(:,1))

for i=1:length(correct1)
    c1stat(i,1) = min(correct1(:,i));
    c1stat(i,2) = max(correct1(:,i));
    c1stat(i,3) = mean(correct1(:,i));
    c1stat(i,4) = median(correct1(:,i));
    c1stat(i,5) = mode(correct1(:,i));
    c1stat(i,6) = std(correct1(:,i));
    c1stat(i,7) = var(correct1(:,i));
end

for i=1:length(correct2)
    c2stat(i,1) = min(correct2(:,i));
    c2stat(i,2) = max(correct2(:,i));
    c2stat(i,3) = mean(correct2(:,i));
    c2stat(i,4) = median(correct2(:,i));
    c2stat(i,5) = mode(correct2(:,i));
    c2stat(i,6) = std(correct2(:,i));
    c2stat(i,7) = var(correct2(:,i));
end

for i=1:length(correct3)
    c3stat(i,1) = min(correct3(:,i));
    c3stat(i,2) = max(correct3(:,i));
    c3stat(i,3) = mean(correct3(:,i));
    c3stat(i,4) = median(correct3(:,i));
    c3stat(i,5) = mode(correct3(:,i));
    c3stat(i,6) = std(correct3(:,i));
    c3stat(i,7) = var(correct3(:,i));
end

for i=1:length(incorrect1)
    i1stat(i,1) = min(incorrect1(:,i));
    i1stat(i,2) = max(incorrect1(:,i));
    i1stat(i,3) = mean(incorrect1(:,i));
    i1stat(i,4) = median(incorrect1(:,i));
    i1stat(i,5) = mode(incorrect1(:,i));
    i1stat(i,6) = std(incorrect1(:,i));
    i1stat(i,7) = var(incorrect1(:,i));
end

for i=1:length(incorrect2)
    i2stat(i,1) = min(incorrect2(:,i));
    i2stat(i,2) = max(incorrect2(:,i));
    i2stat(i,3) = mean(incorrect2(:,i));
    i2stat(i,4) = median(incorrect2(:,i));
    i2stat(i,5) = mode(incorrect2(:,i));
    i2stat(i,6) = std(incorrect2(:,i));
    i2stat(i,7) = var(incorrect2(:,i));
end

for i=1:length(incorrect3)
    i3stat(i,1) = min(incorrect3(:,i));
    i3stat(i,2) = max(incorrect3(:,i));
    i3stat(i,3) = mean(incorrect3(:,i));
    i3stat(i,4) = median(incorrect3(:,i));
    i3stat(i,5) = mode(incorrect3(:,i));
    i3stat(i,6) = std(incorrect3(:,i));
    i3stat(i,7) = var(incorrect3(:,i));
end