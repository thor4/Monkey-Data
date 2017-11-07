time = 0:0.001:1;
frequency = 6;
%phase = 97;
%phase_in_rad = degtorad(phase);
%y = sin(2 * pi * frequency * time + phase_in_rad);
y = sin(2 * pi * frequency * time);
%plot(time, y), xlabel('Samples'), ylabel('Sine wave')
plot(time, y,'LineWidth',5,'Color',[0 .271 .937]);
axis off