% New analysis 12/10/20
% Work on raw time-frequency power from wavelets_power.m script
% Files: mAgoodR1_pow_erp_lfp.mat & mBgoodR1_pow_erp_lfp located in
% OneDrive\Documents\PhD @ FAU\research\High Frequency FP Activity in VWM\data\step 2 TF - power

%init vars
min_freq = 4; %in Hz (need several cycles in an epoch, these epochs are 500ms min so 4Hz = 2 cycles)
max_freq = 100; %nothing above 100
num_frex = 35; %50 for 200Hz, 35 for 100Hz, better for statistics mult comp corr, less smooth spectrogram
frex = logspace(log10(min_freq),log10(max_freq),num_frex); %total num of freq's
% vector of time points to save in post-analysis downsampling
times2save = -400:10:1466; % in ms, 1466 = 505 (sample) + 811 (delay) + 150 (match)

%load non-normalized power files
monkeys = {'mA','mB'}; %setup monkeys array
mi = 2; %choose which monkey: mA=1, mB=2

monkey=monkeys{mi}; %load data for chosen monkey
if monkey=='mA'
    load('mAgoodR1_pow_erp_lfp.mat'); %monkey A data w raw power
    mData = mAgoodR1; clear mAgoodR1 
else
    load('mBgoodR1_pow_erp_lfp.mat'); %monkey B data w raw power
    mData = mBgoodR1; clear mBgoodR1
end

