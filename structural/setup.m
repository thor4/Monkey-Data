%%setup adjacency matrix for analysis & visualizations

%step 1: import adjacency matrix (AM) spreadsheet and save as variable
%step 2: setup the nodes to be in exact order as in 
%prefrontal_parietal_areas spreadsheet, format: source x target

nodes = {'10','9','32','14','25','8B','8Ad','9/46d','46d','46v','9/46v',...
    '8Av','45','47/12','13','11','6DR','PE','PEci','PEc','PEa','PF',...
    'PFop','PFG','IPd','POa','PG','PGm','PGop','Opt'}; %cell variable

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