clear

cmin = 0;
cmax = 0.001;
figure
im = imagesc(adj_p,[cmin cmax])
%im = imagesc(adj_p)
colorbar('Direction','reverse')
% colorbar('Ticks',[.05:.05:cmax],...
%          'TickLabels',{'0','.05 (Significant)','.1','.15','.2','.25','.3','.35','.4','.45','.5'},...
%          'Direction','reverse')
c = colorbar('Ticks',[0,0.0001,.00025,0.0005,.00075,.001],...
         'TickLabels',{'0','.0001 (Sig)','.00025','.0005','.00075','.001'},...
         'Direction','reverse','FontSize',24,'Box','off','TickLength',0.000001)
axis off
colormap(summer)
brighten(.25)
set(c, 'ylim', [0 .1])
%summer brighten(1)
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
     