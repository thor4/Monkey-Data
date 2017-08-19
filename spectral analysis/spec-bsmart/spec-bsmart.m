load('incorrect.mat')
load('correct.mat')
porder = 10;
fs = 1000;
arpower_cBase = zeros(352967,501);
tic
for i = 1:size(correctBase,2)
    x = correctBase(:,i);
    x = transpose(x);
    Nr = size(x,1);
    Nl = size(x,2);
    pp = arpower(x,Nr,Nl,porder,fs);
    arpower_cBase(i,:) = pp(1,:);
end
toc

%only interested in 0-250Hz
arpower_iBase = arpower_iBase(:,1:251);
arpower_cBase = arpower_cBase(:,1:251);

%baseline normalization
arpower_cNorm = arpower_cDelay ./ arpower_cBase;
arpower_iNorm = arpower_iDelay ./ arpower_iBase;

%average across trials
arpower_iNormAvg = mean(arpower_iNorm);
arpower_cNormAvg = mean(arpower_cNorm);