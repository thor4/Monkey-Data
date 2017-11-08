figure
pcolor(adj_p)
colorbar('Direction','reverse')

cmin = 0;
cmax = max(max(adj_p));
figure
im = imagesc(adj_p,[cmin cmax])
colorbar('Direction','reverse')
% colorbar('Ticks',[.05:.05:cmax],...
%          'TickLabels',{'0','.05 (Significant)','.1','.15','.2','.25','.3','.35','.4','.45','.5'},...
%          'Direction','reverse')
colorbar('Ticks',[.05:.05:cmax],...
         'TickLabels',{'0','.05 (Sig)','.1','.15','.2','.25','.3','.35','.4','.45','.5'},...
         'Direction','reverse')
axis off
map = [0, .14, 0
    0, .28 , 0
    0, .42 , 0
    0, .56, 0
    0, .7, 0
    0, .84, 0
    0, .98,0]
colormap(map)
%     0, 0, 0.6
%     0, 0, 0.8
%     0, 0, 1.0];
% blue 0045EF
%  [red=0, green=69, blue=239]
%  red	FF0000
%  [red=255, green=0, blue=0]
%  orange F79256
%  [red=247, green=146, blue=86]

     
figure
contourf(adj_p)
colorbar('Ticks',[.05 .1 .15,.2,.25],...
         'TickLabels',{'.05 (Significant)','.1','.15','.2','.25'},...
         'Direction','reverse')

figure
hm = heatmap(adj_p);
ax  hm.plot;

hm = HeatMap(data);
ax = hm.plot; % 'ax' will be a handle to a standard MATLAB axes.
colorbar(ax); % Turn the colorbar on
colorbar('Direction','reverse')
caxis(ax, [-3 10]); % Adjust the color limits
     