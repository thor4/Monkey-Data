%%setup adjacency matrix for analysis & visualizations

%step 1: import adjacency matrix (AM) spreadsheet and save as variable
%step 2: setup the nodes to be in exact order as in 
%prefrontal_parietal_areas spreadsheet, format: source x target

nodes = ["10","9","32","14","25","8B","8Ad","9/46d","46d","46v","9/46v",...
    "8Av","45","47/12","13","11","6DR","PE","PEci","PEc","PEa","PF",...
    "PFop","PFG","IPd","POa","PG","PGm","PGop","Opt"]; %string variable

%step 3: save AM and nodes as AM.mat data file
load('AM.mat')
% load('AMv.mat') %has 0.5 along diagonal for AM visualization
% following lines are attempts to get AM into format for import into gephi
% EdgeL = adj2gephilab(filename,AM);
% [result] = Matrix2GraphML(AM,'C:\Users\bryan\Documents\GitHub\Monkey-Data\structural\fpn.graphml',nodes);

%identify in-degree, out-degree and in+out=total degree per node
[id,od,deg] = degrees_dir(AM);

% %identify joint degree distribution, J: joint degree distribution matrix 
% %(shifted by one) J_od: number of vertices with od>id. J_id: number of 
% %vertices with id>od. J_bl: number of vertices with id=od.
% [J,J_od,J_id,J_bl] = jdegree(AM);
% 
% %identify density (fraction of present connections to possible connections)
% kden = density_dir(AM);
% 
% %identify clutering coeff (the fraction of triangles around a node (equiv. 
% %the fraction of node's neighbors that are neighbors of each other)
% C = clustering_coef_bd(AM);
% 
% %identify optimal community structure (a subdivision of the network into
% %nonoverlapping groups of nodes which maximizes the number of within-
% %group edges, and minimizes the number of between-group edges. q is 
% %optimized community-structure statistic
% [M,Q] = community_louvain(AM);

%calculate distance matrix (contains lengths of shortest paths between all
%pairs of nodes. An entry (u,v) represents the length of shortest path 
%from node u to node v. The average shortest path length is the 
%characteristic path length of the network. Lengths between disconnected 
%nodes are set to Inf.
D = distance_bin(AM);

%calculate characteristic path length (average shortest path length between
%all pairs of nodes in the network) & global efficiency (average inverse 
%shortest path length in the network. lambda=inf, efficiency=0.5686
[char_path_length,efficiency] = charpath(D); 
% 
% %calculate betweenness centrality (the fraction of all shortest paths in 
% %the network that contain a given node. Nodes with high values of 
% %betweenness centrality participate in a large number of shortest paths.
% BC = betweenness_bin(AM);

%% visualize association matrix
figure(1), clf
% AMv = imagesc(AM);
AMv = imagesc((1:size(AM,2))-0.5, (1:size(AM,1))-0.5, AM);
% G = digraph(AM)
% plot(G)
colormap([1 1 1; 0.75 0.75 0.75; 0.25 0.25 0.25;]) % black 1, white 0
grid on
ax = gca; ax.GridColor = [0 0 0]; ax.LineWidth = 2; % make grid show up
ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
% white [1 1 1] , black [0 0 0], grey [0.5 0.5 0.5]
% colormap(cool)
% colormap(summer)
yticks(1:30); xticks(1:30); yticklabels([]); xticklabels([]);
set(gca,'XAxisLocation','top') %move x-axis to top
ytk=get(gca,'ytick').'; xtik=get(gca,'xtick'); ypos=-0.1; xpos=-.2;
% place y-labels
text(repmat(ypos,size(ytk(1:17))),ytk(1:17),nodes(1:17),'FontSize',16,'Color',[0.4941 0.7294 0.8000],'vertical','bottom','HorizontalAlignment','right') %frontal
text(repmat(ypos,size(ytk(18:30))),ytk(18:30),nodes(18:30),'FontSize',16,'Color',[0.7882 0.4667 0.5608],'vertical','bottom','HorizontalAlignment','right') %parietal
% place x-labels
text(ytk(1:17)-0.85,repmat(xpos,size(ytk(1:17))),nodes(1:17),'FontSize',16,'Color',[0.4941 0.7294 0.8000],'vertical','top','HorizontalAlignment','left','Rotation',90) %frontal
text(ytk(18:30)-0.85,repmat(xpos,size(ytk(18:30))),nodes(18:30),'FontSize',16,'Color',[0.7882 0.4667 0.5608],'vertical','top','HorizontalAlignment','left','Rotation',90) %parietal
% xtickangle(90); 
% set(gca,'TickLength',[0 0],'FontSize',24) %remove ticks
% cbar = colorbar; lim = get(cbar,'Limits'); cbar.TickLabels=(["No","Yes"]); % top/bottom
% cbar.Label.String = 'Connection?'; pos = cbar.Label.Position; 
% cbar.Label.Position=[pos(1)-2.5 pos(2)];
export_fig AM.png -transparent % no background
export_fig AM.eps -transparent % no background
%%in powerpoint make axis labels "from" and "to"

%%change colors of different regional areas to match frontal and parietal
%%color scheme



%% visualize in & out degree distributions - grouped v1
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
%need to turn FaceAlpha to 1 to turn transparency off to deal with
%export_fig bug only for histograms
id_hist = histogram(id_sort,'BinEdges',(2:3:23),'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
h=gca; h.Color = 'none'; % turn off background color
ylim([0,10]); % same-scale y-axis as fig 6
xticks([10 mean(od) 20]); % same tick marks
xline(mean(id),'--','Color',[0.5 0.5 0.5]); 
h.XAxis.TickLength = [0.005 0]; 
h.YAxisLocation = 'Right'; h.XDir = 'reverse'; %reflect about y-axis
box off; h.YAxis.Visible = 'off'; % turn off y-axis
% axis off
export_fig in_deg_hist.eps -transparent % no background

%% testing best cCDF plot
% [xi_deg,yi_pdk,idxi] = cCDF(id); %in-degree *legacy*
% [xo_deg,yo_pdk,idxo] = cCDF(od); %out-degree *legacy*
% get xi/o's and yi/o's from next section's c values, first in-deg
xi_deg = c(:,1); yi_pdk = c(:,2); %in-deg cCDF
% now re-run to get new c for out-degree
xo_deg = c(:,1); yo_pdk = c(:,2); %out-deg cCDF
% xi=c(:,1); yi=c(:,2); %pull out deg & prob from cCDF for curve fitting
%use curve fitting app to find best fit for cCDFs and make fitFPN function
[fpnresult,fpngof] = fitFPN(xi_deg, yi_pdk, xo_deg, yo_pdk); %idx 4 & 8 are best
save('fitted_models','fpnresult'); save('gof_models','fpngof') %save models and their gof

% now pass the xi/o_deg into each model to get its prediction
% yi_poly = fpnresult{1}(xi_deg); %Polynomial best fit of cCDF in-deg fpn
yi_pow = fpnresult{2}(xi_deg); %Power best fit of cCDF in-deg fpn
yi_exp = fpnresult{3}(xi_deg); %Exponential best fit of cCDF in-deg fpn
yi_gauss = fpnresult{4}(xi_deg); %Gaussian best fit of cCDF in-deg fpn
% yo_poly = fpnresult{1}(xo_deg); %Polynomial best fit of cCDF out-deg fpn
yo_pow = fpnresult{2}(xo_deg); %Power best fit of cCDF out-deg fpn
yo_exp = fpnresult{3}(xo_deg); %Exponential best fit of cCDF out-deg fpn
yo_gauss = fpnresult{4}(xo_deg); %Gaussian best fit of cCDF out-deg fpn


figure(9), clf %cCDF in-deg fits
loglog(xi_deg,yi_pdk,'ko','MarkerSize',8,'MarkerFaceColor',[1 1 1]); %in
hold on
loglog(xi_deg,yi_exp,'b:','LineWidth',2); %exp fit
loglog(xi_deg,yi_pow,'r--','LineWidth',2); %power law fit
loglog(xi_deg,yi_gauss,'m-.','LineWidth',2); %gaussian fit
hold off; 
set(gca,'XLim',[1,30],'XTick',10.^xrt);
set(gca,'YLim',[0.03,1.25],'YTick',10.^yrt,'FontSize',16);
ylabel('P(degree \geq x)','FontSize',18); xlabel('x','FontSize',18)
title('In-degree Distribution cCDF','FontSize',20); %update for id/od/deg accordingly
export_fig id_ccdf_matlab.eps -transparent % no background
export_fig id_ccdf_matlab.png -transparent % no background

figure(10), clf %cCDF out-deg fits
loglog(xo_deg,yo_pdk,'ko','MarkerSize',8,'MarkerFaceColor',[1 1 1]); %out
hold on
loglog(xo_deg,yo_exp,'b:','LineWidth',2); %exp fit
loglog(xo_deg,yo_pow,'r--','LineWidth',2); %power law fit
loglog(xo_deg,yo_gauss,'m-.','LineWidth',2); %gaussian fit
hold off; 
set(gca,'XLim',[2,30],'XTick',10.^xrt);
set(gca,'YLim',[0.03,1.25],'YTick',10.^yrt,'FontSize',16);
ylabel('P(degree \geq x)','FontSize',18); xlabel('x','FontSize',18)
title('Out-degree Distribution cCDF','FontSize',20); %update for id/od/deg accordingly
export_fig od_ccdf_matlab.eps -transparent % no background
export_fig od_ccdf_matlab.png -transparent % no background


% %use this if you want to find degrees where its like a power law dist,
% %between degrees [6,14]
% fit_xdeg = xi_deg(3:15)'; %remove outliers for in-deg with no power dist
% fit_ypdk = yi_pdk(3:15)';
% [fp,gofp,outputp] = fit(fit_xdeg,fit_ypdk,'power1'); % y = a*x^b
% %in sse: 0.0143 r^2: 0.9149 adjr^2: 0.9063 rmse: 0.0379 [lin reg]
% %fp.b = -0.5433 (not in barabasi's range)
% yp = fp(fit_xdeg); %power law best fit
% hold on
% loglog(fit_xdeg,yp,'g-')
% plot(f_pow,'predobs')

% %frontal & parietal in+out degree analyses (worry about this for SI)
% frontal = AM(1:17,1:17); %pull out all frontal areas
% parietal = AM(18:30,18:30); %pull out all parietal areas
% %identify in-degree, out-degree and in+out=total degree per node
% [fid,fod,fdeg] = degrees_dir(frontal); %frontal in+out degrees
% [pid,pod,pdeg] = degrees_dir(parietal); %parietal in+out degrees
% [fxi_deg,fyi_pdk,~] = cCDF(fid); [fxo_deg,fyo_pdk,~] = cCDF(fod); %cCDFs frontal
% [pxi_deg,pyi_pdk,~] = cCDF(pid); [pxo_deg,pyo_pdk,~] = cCDF(pod); %cCDFs frontal
% %use curve fitting app to explore cCDFs and find best fit
% [fpresult,fpgof] = fitFrontalParietal(fxi_deg, fyi_pdk, fxo_deg, fyo_pdk,...
%     pxi_deg, pyi_pdk, pxo_deg, pyo_pdk); %idx 2, 3, 5 & 7 are best
% yif = fpresult{2}(fxi_deg); %Gaussian best fit of cCDF in-deg frontal
% 
% figure(9), clf %cCDF of in-deg frontal
% loglog(fxi_deg,fyi_pdk,'ro'); %in
% hold on
% loglog(fxi_deg,yif,'k-'); %in-deg frontal Gaussian best-fit

%see about other metrics to explore
%make structural figures and write up methods + discussion surrounding it

%% test for power law distribution
% files located here: http://tuvalu.santafe.edu/~aaronc/powerlaws/
[alpha, xmin, L] = plfit(id','finite');
% The output 'alpha' is the maximum likelihood estimate of the scaling
% exponent, 'xmin' is the estimate of the lower bound of the power-law
% behavior, and L is the log-likelihood of the data x>=xmin under the
% fitted power law.
[p, gof] = plpva(id',xmin,'reps',5000); %p=0.0158, reject power law hypo for id deg, not drawn from a power-law dist

[alpha, xmin, L] = plfit(od','finite');
[p, gof] = plpva(od',xmin,'reps',5000); %p=0.1716, >0.1 so power law is plausible hypo for data

% now go to python for curve fitting comparisons
% find in-degree is better approximated by an exponential than a power law

% find out-degree

%visualize od dist along with fitted power-law dist on log-log axes
% h = plplot(od,xmin,alpha); 

%% visualize id+od dist along with fitted power-law dist on log-log axes
% 99% of this taken directly from plplot
% reshape input vector
x = reshape(id,numel(id),1); %state whether looking at in/out/total deg
% initialize storage for output handles
h = zeros(2,1);

% select method (discrete or continuous) for plotting
if     isempty(setdiff(x,floor(x))), f_dattype = 'INTS';
elseif isreal(x),    f_dattype = 'REAL';
else                 f_dattype = 'UNKN';
end;
if strcmp(f_dattype,'INTS') && min(x) > 50,
    f_dattype = 'REAL';
end;

% calc cCDF (c) and best-fit line (cf)
switch f_dattype,
    
    case 'REAL',
        n = length(x);
        c = [sort(x) (n:-1:1)'./n];
        q = sort(x(x>=xmin));
        cf = [q (q./xmin).^(1-alpha)];
        cf(:,2) = cf(:,2) .* c(find(c(:,1)>=xmin,1,'first'),2);
        xi= c(9:end,1); %start from xmin (12)
        yi= [1 0.790343 0.624642 0.493681 0.390177 0.308374 0.243721 ...
            0.192623 0.120321]; %exp fit from python from xmin (12) to 21
%         xo= c(9:end,1); %start from xmin (10)
%         yo= [1 0.821421 0.674732 0.554239 0.455264 0.373963 0.307181 ...
%             0.252325 0.207265 0.170252 0.139848 0.114874 0.0775095]; 

        figure;
        h(1) = loglog(c(:,1),c(:,2),'ko','MarkerSize',8,'MarkerFaceColor',[1 1 1]); hold on;
        h(2) = loglog(cf(:,1),cf(:,2),'r--','LineWidth',2); %power law fit
        h(3) = loglog(xi,yi,'b:','LineWidth',2); hold off; %exp fit
%         h(3) = loglog(xo,yo,'b:','LineWidth',2); hold off; %exp fit
        xr  = [10.^floor(log10(min(x))) 10.^ceil(log10(max(x)))];
        xrt = (round(log10(xr(1))):2:round(log10(xr(2))));
        if length(xrt)<4, xrt = (round(log10(xr(1))):1:round(log10(xr(2)))); end;
        yr  = [10.^floor(log10(1/n)) 1];
        yrt = (round(log10(yr(1))):2:round(log10(yr(2))));
        if length(yrt)<4, yrt = (round(log10(yr(1))):1:round(log10(yr(2)))); end;
        set(gca,'XLim',xr,'XTick',10.^xrt);
        set(gca,'YLim',yr,'YTick',10.^yrt,'FontSize',16);
        ylabel('P(degree \geq x)','FontSize',18);
        xlabel('x','FontSize',18)

    case 'INTS',
        n = length(x);        
        q = unique(x);
        c = hist(x,q)'./n; %c is p(degree=x)
        c = [[q; q(end)+1] 1-[0; cumsum(c)]]; c(c(:,2)<10^-10,:) = []; %now c is p(degree>=x)
        cf = ((xmin:q(end))'.^-alpha)./(zeta(alpha) - sum((1:xmin-1).^-alpha));
        cf = [(xmin:q(end)+1)' 1-[0; cumsum(cf)]];
        cf(:,2) = cf(:,2) .* c(c(:,1)==xmin,2);
        xi= c(9:end,1); %start from xmin (12)
        yi= [1 0.790343 0.624642 0.493681 0.390177 0.308374 0.243721 ...
            0.192623 0.120321]; %exp fit from python from xmin (12) to 21
        %exp fit from python from xmin (10) to 23
%         xo= c(6:end,1); %start from xmin (10) 
%         yo= [1 0.821421 0.674732 0.554239 0.455264 0.373963 0.307181 ...
%             0.252325 0.207265 0.170252 0.139848 0.114874 0.0775095]; 


        figure;
        h(1) = loglog(c(:,1),c(:,2),'ko','MarkerSize',8,'MarkerFaceColor',[1 1 1]); hold on;
        h(2) = loglog(cf(:,1),cf(:,2),'r--','LineWidth',2); %power law fit
        h(3) = loglog(xi,yi,'b:','LineWidth',2); hold off; %exp fit
%         h(3) = loglog(xo,yo,'b:','LineWidth',2); hold off; %exp fit
        xr  = [10.^floor(log10(min(x))) 10.^ceil(log10(max(x)))];
        xrt = (round(log10(xr(1))):2:round(log10(xr(2))));
        if length(xrt)<4, xrt = (round(log10(xr(1))):1:round(log10(xr(2)))); end;
        yr  = [10.^floor(log10(1/n)) 1];
        yrt = (round(log10(yr(1))):2:round(log10(yr(2))));
        if length(yrt)<4, yrt = (round(log10(yr(1))):1:round(log10(yr(2)))); end;
        set(gca,'XLim',xr,'XTick',10.^xrt);
        set(gca,'YLim',yr,'YTick',10.^yrt,'FontSize',16);
        ylabel('P(degree \geq x)','FontSize',18);
        xlabel('x','FontSize',18)

    otherwise,
        fprintf('(PLPLOT) Error: x must contain only reals or only integers.\n');
        h = [];
        return;
end;



set(gca,'XLim',[1,30],'XTick',10.^xrt);
set(gca,'YLim',[0.03,1.25],'YTick',10.^yrt,'FontSize',16);
ylabel('P(degree \geq x)','FontSize',18); xlabel('x','FontSize',18)
title('In-degree Distribution cCDF','FontSize',20); %update for id/od/deg accordingly
% legend('empirical data','power-law fit','exponential fit')
export_fig id_ccdf.eps -transparent % no background
export_fig id_ccdf.png -transparent % no background




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
export_fig out_deg_lt.eps -transparent % no background

figure(6), clf %histogram
od_hist = histogram(od_sort,'BinEdges',(2:3:23),'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
h=gca; h.Color = 'none'; % turn off background color
ylim([0,10]); % same-scale y-axis as fig 6
xticks([10 mean(od) 20]); % same tick marks
xline(mean(id),'--','Color',[0.5 0.5 0.5]); 
h.XAxis.TickLength = [0.005 0]; 
box off; h.YAxis.Visible = 'off'; % turn off y-axis
% axis off
export_fig out_deg_hist.eps -transparent % no background


%% surrogate networks

% %be sure not to include the diagonal in the mean distance calculation
% rowMean = sum(D,2) ./ sum(D~=0,2); %avg shortest path length
% char_path_length = mean(rowMean); %characteristic path length

%random Maslov-Sneppen re-wiring
% ..the number of iterations should exceed at least 100 times the number of 
% connections in the network (Milo et al., 2004). 100*30 = 3000
iter = 100; % # of iterations
networks = 100; % # of surrogate networks
ensemble = zeros(size(AM,1),size(AM,2),networks); %init ensemble
C_ensemble = zeros(networks,1); %init clust coef ensemble
L_ensemble = zeros(networks,1); %init char path length ensemble
f_ensemble = zeros(13,networks); %init network motif freq fingerprint ensemble
% doc randmio_dir
tic
for i=1:networks
    % R: randomized network, eff: number of actual rewirings carried out
    [R eff] = randmio_dir(AM, iter); 
    ensemble(:,:,i) = R; %build ensemble of surrogate networks
    C_ensemble(i) = mean(clustering_coef_bd(R)); 
    D_rand = distance_bin(R); %shortest path length for each node
    [L_ensemble(i),~] = charpath(D_rand); 
    [f_ensemble(:,i),~]=motif3struct_bin(R); %network motif freq fingerprint
    clear R
end
toc
% 19.14 seconds = 100 surrogate networks
% [id_test,od_test,deg_test] = degrees_dir(ensemble(:,:,57)); %test deg dist
% id==id_test
% od==od_test
% all deg dist are maintained

%lattice

%Rlatt,  latticized network in original node ordering
%Rrp, latticized network in node ordering used for latticization
%ind_rp, node ordering used for latticization
%eff, number of actual rewirings carried out
latt_ensemble = zeros(size(AM,1),size(AM,2),networks); %init ensemble
C_latt_ensemble = zeros(networks,1); %init clust coef ensemble

tic
for i=1:networks
    % R: randomized network, eff: number of actual rewirings carried out
    [Rlatt,Rrp,ind_rp,eff] = latmio_dir_connected(AM, iter); 
    latt_ensemble(:,:,i) = Rrp; %build ensemble of surrogate networks
    C_latt_ensemble(i) = mean(clustering_coef_bd(Rrp)); 
    D_latt = distance_bin(Rrp); %shortest path length for each node
    clear Rrp
end
toc
%25 seconds = 100 surrogate networks

%% Motif analysis
%first must make_motif34lib.m which generates the required motif34lib.mat
%which all motif releated functions require
make_motif34lib

%next be sure to generate null hypothesis distributions through random and
%lattice-based network ensemble generation

[f,F]=motif3struct_bin(AM); % real data
% The motif frequency of occurrence around an individual node is known as
% the motif fingerprint of that node. The total motif frequency of
% occurrence in the entire network is known as the motif fingerprint of the network. 
%F is node motif frequency fingerprint
%f is network motif frequency fingerprint
f_rand = mean(f_ensemble,2); %compute the mean motif fingerprint over all randomly generated networks
f_std = std(f_ensemble,0,2); %compute the std dev over all rand gen networks
f_zscores = (f - f_rand) > (0.1 * f_rand); % Milo, 2002 method, only 9 & 13 sig

%% Small-world analysis
%run surrogate networks first to generate ensembles

%Humphries et al. (2006& 2008) method of calculating sigma:
%only uses random networks
C = mean(clustering_coef_bd(AM)); % avg clustering coef for empirical network
D = distance_bin(AM); [L,efficiency] = charpath(D); %L = char path length
gamma = C / mean(C_ensemble); %Gamma > 1 suggests greater clustering than random
lambda = L / mean(L_ensemble); %Lambda ~ 1 suggests a comparable avg path length to randomized network
sigma = gamma / lambda; %sigma > 1 indicates small-worldness

%Telesford et al., 2011 technique:
%uses both random and latticized networks to capture both sides of the
%spectrum
omega = (mean(L_ensemble) / L) - (C / mean(C_latt_ensemble));
%omega = 0.0275, 
%omega index ranges between -1 and 1. Values close to 0 are indicative of 
%small-worldness. Positive values suggest more random characteristics and 
%negative values indicate a lattice-like structure