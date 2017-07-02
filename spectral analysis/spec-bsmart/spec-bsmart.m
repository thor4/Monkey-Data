load('C:\Users\bryan\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\incorrect.mat')
load('C:\Users\bryan\OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\correct.mat')
porder = 10;
fs = 1000;
arpower_i = zeros(137078,501);
tic
for i = 1:size(incorrect,2)
    x = incorrect(:,i);
    x = transpose(x);
    Nr = size(x,1);
    Nl = size(x,2);
    pp = arpower(x,Nr,Nl,porder,fs);
    arpower_i(i,:) = pp(1,:);
end
toc

