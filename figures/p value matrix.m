clear

cmin = 0;
% cmax = 0.001;
cmax = max(max(adj_p));
figure
im = imagesc(adj_p,[cmin cmax])
%im = imagesc(adj_p)
colorbar('Direction','reverse')
colorbar('Ticks',[0:.05:cmax],...
         'TickLabels',{'0','.05 (Sig)','.1','.15','.2','.25','.3','.35','.4','.45','.5'},...
         'Direction','reverse','FontSize',24,'Box','off','TickLength',0.000001)
% c = colorbar('Ticks',[0,0.0001,.00025,0.0005,.00075,.001],...
%          'TickLabels',{'0','.0001 (Sig)','.00025','.0005','.00075','.001'},...
%          'Direction','reverse','FontSize',24,'Box','off','TickLength',0.000001)
axis off
colormap(summer)
brighten(.25)

% blue 0045EF
%  [red=0, green=69, blue=239]
%  red	FF0000
%  [red=255, green=0, blue=0]
%  orange F79256
%  [red=247, green=146, blue=86]