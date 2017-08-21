load('incorrect.mat')
load('correct.mat')
porder = 10;
fs = 1000;
arpower_iBase = zeros(137078,501);
tic
for i = 1:size(incorrectBase,2)
    x = incorrectBase(:,i);
    x = transpose(x);
    Nr = size(x,1);
    Nl = size(x,2);
    pp = arpower(x,Nr,Nl,porder,fs);
    arpower_iBase(i,:) = pp(1,:);
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

%second option: average then normalize
arpower_cBaseAvg = mean(arpower_cBase);
arpower_cDelayAvg = mean(arpower_cDelay);
arpower_cAvgNorm = arpower_cDelayAvg ./ arpower_cBaseAvg;