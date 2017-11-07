figure
pcolor(adj_p)
colorbar('Direction','reverse')

cmin = 0;
cmax = 0.25;
figure
imagesc(adj_p,[cmin cmax])
colorbar('Ticks',[.05 .1 .15,.2,.25],...
         'TickLabels',{'.05 (Significant)','.1','.15','.2','.25'},...
         'Direction','reverse')
     
figure
contourf(adj_p)
colorbar('Ticks',[.05 .1 .15,.2,.25],...
         'TickLabels',{'.05 (Significant)','.1','.15','.2','.25'},...
         'Direction','reverse')

figure
heatmap(adj_p)
colorbar('Direction','reverse')
     