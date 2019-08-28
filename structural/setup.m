%%setup adjacency matrix for analysis & visualizations

%step 1: import adjacency matrix (AM) spreadsheet and save as variable
%step 2: setup the nodes to be in exact order as in 
%prefrontal_parietal_areas spreadsheet, format: source x target

nodes = ["10","9","32","14","25","8B","8Ad","9/46d","46d","46v","9/46v",...
    "8Av","45","47/12","13","11","6DR","PE","PEci","PEc","PEa","PF",...
    "PFop","PFG","IPd","POa","PG","PGm","PGop","Opt"]; %string variable

%step 3: save AM and nodes as AM.mat data file
load('AM.mat')

%identify in-degree, out-degree and in+out=total degree per node
[id,od,deg] = degrees_dir(AM);

%identify joint degree distribution, J: joint degree distribution matrix 
%(shifted by one) J_od: number of vertices with od>id. J_id: number of 
%vertices with id>od. J_bl: number of vertices with id=od.
[J,J_od,J_id,J_bl] = jdegree(AM);

%identify density (fraction of present connections to possible connections)
kden = density_dir(AM);

%identify clutering coeff (the fraction of triangles around a node (equiv. 
%the fraction of node's neighbors that are neighbors of each other)
C = clustering_coef_bd(AM);

%identify optimal community structure (a subdivision of the network into
%nonoverlapping groups of nodes which maximizes the number of within-
%group edges, and minimizes the number of between-group edges. q is 
%optimized community-structure statistic
[M,Q] = community_louvain(AM);

%calculate distance matrix (contains lengths of shortest paths between all
%pairs of nodes. An entry (u,v) represents the length of shortest path 
%from node u to node v. The average shortest path length is the 
%characteristic path length of the network. Lengths between disconnected 
%nodes are set to Inf.
D = distance_bin(AM);

%calculate characteristic path length (average shortest path length between
%all pairs of nodes in the network) & global efficiency (average inverse 
%shortest path length in the network. lambda=inf, efficiency=0.5686
[lambda,efficiency] = charpath(D); %lambda=char path length

%calculate betweenness centrality (the fraction of all shortest paths in 
%the network that contain a given node. Nodes with high values of 
%betweenness centrality participate in a large number of shortest paths.
BC = betweenness_bin(AM);

%% visualize association matrix
imagesc(AM)
% colormap(cool)
colormap(summer)
yticks(1:30); xticks(1:30); yticklabels(nodes(:)); xticklabels(nodes(:));
set(gca,'XAxisLocation','top') %move x-axis to top
xtickangle(90); set(gca,'TickLength',[0 0],'FontSize',24) %remove ticks
cbar = colorbar; lim = get(cbar,'Limits'); cbar.TickLabels=(["No","Yes"]); % top/bottom
set(get(cbar,'label'),'string','Raw power (\muV^2)');
cbar.Label.String = 'Connection?'; pos = cbar.Label.Position; 
cbar.Label.Position=[pos(1)-2.5 pos(2)];
export_fig AM.png -transparent % no background
%%make font size bigger then export, in powerpoint make axis labels "from"
%%and "to"

% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',(1:30),'ytick',(1:30),'TickDir','out');
% set original axes as active
axes(a)
% link axes in case of zooming
linkaxes([a b])
yticks('auto')

%% visualize in & out degree distributions
[id_sort,id_idx] = sort(id); %sort elements of id in ascending order and save indices in idx
od_id_sort = od(id_idx); %od sorted by the sorted id vector
deg=horzcat(id_sort',od_id_sort'); %concatenate id + od to have single deg mat
figure(2), clf
bh = barh(1:30,deg,'grouped','FaceColor','flat'); 
% bh(1).FaceColor = '#C9778F'; %in-deg
% bh(2).FaceColor = '#7EBACC'; %out-deg
% bh(1).CData(1,:) = '#C9778F';
% stem(1:30,id_sort)
yticks(1:30); yticklabels(nodes(id_idx));
xticks([10 mean(id) 20]); 
h=gca; h.XAxis.TickLength = [0.005 0]; % shorten tick marks on only x-axis
h.YAxis.TickLength = [0 0]; % del tick marks on only y-axis
for i=1:length(id_idx) %color each tick label & bar's area acc to region
    if (id_idx(i) < 18) %frontal area #7EBACC
        h.YTickLabel{i} = ['\color[rgb]{0.4941 0.7294 0.8000}' h.YTickLabel{i}];
        bh(1).CData(i,:) = [0.4941 0.7294 0.8000]; %color both bars
        bh(2).CData(i,:) = [0.4941 0.7294 0.8000]; %color both bars
    else %parietal area #C9778F
        h.YTickLabel{i} = ['\color[rgb]{0.7882 0.4667 0.5608}' h.YTickLabel{i}];
        bh(1).CData(i,:) = [0.7882 0.4667 0.5608]; %color both bars
        bh(2).CData(i,:) = [0.7882 0.4667 0.5608]; %color both bars
    end
end
bh(1).LineStyle = '-'; %set in-degree line style to solid
bh(2).LineStyle = ':'; %set out-degree line style to dotted
bh(1).LineWidth = 0.25; bh(2).LineWidth = 0.25; %set line size

xline(mean(id),'--','Mean In-Degree','LabelVerticalAlignment','bottom',...
    'Color','#C9778F');
xline(mean(od),'--','Mean Out-Degree','LabelVerticalAlignment','bottom',...
    'LabelHorizontalAlignment','left','Color','#7EBACC');
hold on %add phantom plots to get legend to present correctly
sol=plot(21,2,'k-'); dot=plot(21,3,'k:'); hold off
legend([sol dot],'In-degree','Out-degree'); title('In & Out-Degree Distributions')
legend('boxoff')
export_fig deg_dist.pdf -transparent % no background
export_fig deg_dist.png -transparent % no background

%% in degree distribution visualization v2
[id_sort,id_idx] = sort(id); %sort elements of id in ascending order and save indices in idx
figure(3), clf
bh = barh(1:30,id_sort,'FaceColor','flat'); 
yticks(1:30); yticklabels(nodes(id_idx));
xticks([10 mean(id) 20]); 
h=gca; h.XAxis.TickLength = [0.005 0]; % shorten tick marks on only x-axis
h.YAxis.TickLength = [0 0]; % del tick marks on only y-axis
h.YAxisLocation = 'Right'; h.XDir = 'reverse'; %reflect about y-axis
for i=1:length(id_idx) %color each tick label & bar's area acc to region
    if (id_idx(i) < 18) %frontal area #6DB3A5
        h.YTickLabel{i} = ['\color[rgb]{0.4941 0.7294 0.8000}' h.YTickLabel{i}];
        bh.CData(i,:) = [0.4941 0.7294 0.8000]; %color bars
    else %parietal area #C9778F
        h.YTickLabel{i} = ['\color[rgb]{0.7882 0.4667 0.5608}' h.YTickLabel{i}];
        bh.CData(i,:) = [0.7882 0.4667 0.5608]; %color both bars
    end
end
xline(mean(id),'--','Mean','LabelVerticalAlignment','bottom',...
    'Color',[0.5 0.5 0.5]);
h.Color = 'none'; % turn off background color
box off
export_fig in_deg_rt.eps -transparent % no background

figure(4), clf %histogram
id_hist = histogram(id_sort,6,'FaceColor',[0.5 0.5 0.5]);
xline(mean(id),'--','Color',[0.5 0.5 0.5]); 
axis off
h.Color = 'none'; % turn off background color


%% out degree distribution visualization v2
[od_sort,od_idx] = sort(od); %sort elements of id in ascending order and save indices in idx
figure(5), clf
bh = barh(1:30,od_sort,'FaceColor','flat'); 
yticks(1:30); yticklabels(nodes(od_idx));
xticks([10 mean(od) 20]); 
h=gca; h.XAxis.TickLength = [0.005 0]; % shorten tick marks on only x-axis
h.YAxis.TickLength = [0 0]; % del tick marks on only y-axis
for i=1:length(od_idx) %color each tick label & bar's area acc to region
    if (od_idx(i) < 18) %frontal area #6DB3A5
        h.YTickLabel{i} = ['\color[rgb]{0.4941 0.7294 0.8000}' h.YTickLabel{i}];
        bh.CData(i,:) = [0.4941 0.7294 0.8000]; %color bars
    else %parietal area #C9778F
        h.YTickLabel{i} = ['\color[rgb]{0.7882 0.4667 0.5608}' h.YTickLabel{i}];
        bh.CData(i,:) = [0.7882 0.4667 0.5608]; %color both bars
    end
end
xline(mean(od),'--','Mean','LabelVerticalAlignment','bottom',...
    'Color',[0.5 0.5 0.5]);
h.Color = 'none'; % turn off background color
box off

figure(6), clf %histogram
id_hist = histogram(od_sort,6,'FaceColor',[0.5 0.5 0.5]);
xline(mean(od),'--','Color',[0.5 0.5 0.5]); 
h.Color = 'none'; % turn off background color
axis off

