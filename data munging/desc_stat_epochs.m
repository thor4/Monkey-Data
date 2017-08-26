%initialize descriptive statistics matrices
c1stat = zeros(length(cSub1),7);
c2stat = zeros(length(cSub2),7);
c3stat = zeros(length(cSub3),7);

i1stat = zeros(length(iSub1),7);
i2stat = zeros(length(iSub2),7);
i3stat = zeros(length(iSub3),7);

test = [1, -7, 3; -4, 5, 6];
max(abs(test(1,:)))
max(test(1,:))

same = c1stat(:,1) == c1stat(:,5);
find(same==0)


for i=1:length(cSub1)
    c1stat(i,1) = min(cSub1(i,:));
    c1stat(i,2) = max(cSub1(i,:));
    c1stat(i,3) = mean(cSub1(i,:));
    c1stat(i,4) = median(cSub1(i,:));
    c1stat(i,5) = std(cSub1(i,:));
    c1stat(i,6) = var(cSub1(i,:));
    c1stat(i,7) = max(abs(cSub1(i,:)));
end

for i=1:length(cSub2)
    c2stat(i,1) = min(cSub2(i,:));
    c2stat(i,2) = max(cSub2(i,:));
    c2stat(i,3) = mean(cSub2(i,:));
    c2stat(i,4) = median(cSub2(i,:));
    c2stat(i,5) = std(cSub2(i,:));
    c2stat(i,6) = var(cSub2(i,:));
    c2stat(i,7) = max(abs(cSub2(i,:)));
end

for i=1:length(cSub3)
    c3stat(i,1) = min(cSub3(i,:));
    c3stat(i,2) = max(cSub3(i,:));
    c3stat(i,3) = mean(cSub3(i,:));
    c3stat(i,4) = median(cSub3(i,:));
    c3stat(i,5) = std(cSub3(i,:));
    c3stat(i,6) = var(cSub3(i,:));
    c3stat(i,7) = max(abs(cSub3(i,:)));
end

for i=1:length(iSub1)
    i1stat(i,1) = min(iSub1(i,:));
    i1stat(i,2) = max(iSub1(i,:));
    i1stat(i,3) = mean(iSub1(i,:));
    i1stat(i,4) = median(iSub1(i,:));
    i1stat(i,5) = std(iSub1(i,:));
    i1stat(i,6) = var(iSub1(i,:));
    i1stat(i,7) = max(abs(iSub1(i,:)));
end

for i=1:length(iSub2)
    i2stat(i,1) = min(iSub2(i,:));
    i2stat(i,2) = max(iSub2(i,:));
    i2stat(i,3) = mean(iSub2(i,:));
    i2stat(i,4) = median(iSub2(i,:));
    i2stat(i,5) = std(iSub2(i,:));
    i2stat(i,6) = var(iSub2(i,:));
    i2stat(i,7) = max(abs(iSub2(i,:)));
end

for i=1:length(iSub3)
    i3stat(i,1) = min(iSub3(i,:));
    i3stat(i,2) = max(iSub3(i,:));
    i3stat(i,3) = mean(iSub3(i,:));
    i3stat(i,4) = median(iSub3(i,:));
    i3stat(i,5) = std(iSub3(i,:));
    i3stat(i,6) = var(iSub3(i,:));
    i3stat(i,7) = max(abs(iSub3(i,:)));
end