load('C:\Users\bryan\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\incorrect.mat')
load('C:\Users\bryan\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\correct.mat')
porder = 10;
fs = 1000;
arpower_cBase = zeros(352967,501);
tic
for i = 1:size(correct,2)
    x = correct(:,i);
    x = transpose(x);
    Nr = size(x,1);
    Nl = size(x,2);
    pp = arpower(x,Nr,Nl,porder,fs);
    arpower_cBase(i,:) = pp(1,:);
end
toc

%only interested in 0-250Hz
arpower_i = arpower_i(:,1:251);
arpower_c = arpower_c(:,1:251);
arpower_iBase = arpower_iBase(:,1:251);
arpower_cBase = arpower_cBase(:,1:251);