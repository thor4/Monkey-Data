%power spectrum
doc plot_vector
figure
subplot(2,1,1)
plot_vector(ScDelay(:,10),fcDelay,'n',[],'b')
hold on
fill(theta,ScDelay(thetaIdx(1):thetaIdx(end),10),'r');
hold on
fill(fcDelay(1):fcDelay(5),ScDelay(1,10):ScDelay(length(fcDelay(1):fcDelay(5)),10),'r')
subplot(2,1,2)
plot_vector(SiDelay(:,10),fiDelay,'n',[],'g')
delta = fcDelay(fcDelay>0 & fcDelay<4);
theta = fcDelay(fcDelay>4 & fcDelay<8);
alpha = 8:12;
beta = 12:30;
low_gamma = 30:55;
high_gamma1 = 65:115;
high_gamma2 = 125:175;
high_gamma3 = 185:235;

freq = [ScDelay(deltaIdx(1):deltaIdx(end),10); ScDelay(deltaIdx(end):thetaIdx(end),10)];
h = area(fcDelay(deltaIdx(1):deltaIdx(end)),freq);

figure
h = area(fcDelay(deltaIdx(1):deltaIdx(end)),ScDelay(deltaIdx(1):deltaIdx(end),10));
hold on
area(fcDelay(deltaIdx(end):thetaIdx(end)),ScDelay(deltaIdx(end):thetaIdx(end),10));
plot(fcDelay,ScDelay(:,10));

Y = [1, 5, 3;
     3, 2, 7;
     1, 5, 3;
     2, 6, 1];
z = area(Y,'LineStyle',':');

figure
%fill(theta,ScDelay(thetaIdx(1):thetaIdx(end),10),'r');
area(theta,ScDelay(thetaIdx(1):thetaIdx(end),10),'r');
patch(theta,ScDelay(thetaIdx(1):thetaIdx(end),10),'r');

deltaIdx = find(fcDelay>=0 & fcDelay<4);
thetaIdx = find(fcDelay>=4 & fcDelay<8);

x = (-10:0.1:10);
xs = x(x>-4 & x<4);
figure;
hold on;
area(xs,normpdf(xs,0,3));
plot(x,normpdf(x,0,3));