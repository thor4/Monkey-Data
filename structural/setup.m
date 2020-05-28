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

%total number of edges (connections) in network
edges = sum(sum(AM)); sum(deg)/2, sum(id), sum(od) %all yield 399

% %identify joint degree distribution, J: joint degree distribution matrix 
% %(shifted by one) J_od: number of vertices with od>id. J_id: number of 
% %vertices with id>od. J_bl: number of vertices with id=od.
% [J,J_od,J_id,J_bl] = jdegree(AM);
% 
% %identify density (fraction of present connections to possible connections)
kden = density_dir(AM); % 45.86%
% 
% %identify clutering coeff (the fraction of triangles around a node (equiv. 
% %the fraction of node's neighbors that are neighbors of each other)
% C = clustering_coef_bd(AM);
% 

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

%% visualize association matrix

CM = CIJ; %define connectivity matrix

figure(1), clf
% AMv = imagesc(AM);
CMv = imagesc((1:size(CM,2))-0.5, (1:size(CM,1))-0.5, CM);
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
export_fig toeplitz.png % no background
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
load('fitted_models.mat') %saved fits
% yi_poly = fpnresult{1}(xi_deg); %Polynomial best fit of cCDF in-deg fpn
yi_pow = fpnresult{2}(xi_deg); %Power best fit of cCDF in-deg fpn
yi_exp = fpnresult{3}(xi_deg); %Exponential best fit of cCDF in-deg fpn
yi_gauss = fpnresult{4}(xi_deg); %Gaussian best fit of cCDF in-deg fpn
% yo_poly = fpnresult{1}(xo_deg); %Polynomial best fit of cCDF out-deg fpn
yo_pow = fpnresult{2}(xo_deg); %Power best fit of cCDF out-deg fpn
yo_exp = fpnresult{3}(xo_deg); %Exponential best fit of cCDF out-deg fpn
yo_gauss = fpnresult{4}(xo_deg); %Gaussian best fit of cCDF out-deg fpn


%% Make id + od dist cCDF subplots
%load data and prep for plotting
load('tail_fits.mat') %load the tail fits (exp & pow) for id+od dist

load('fitted_models.mat') %load the full fit models and calc y's
% yi_poly = fpnresult{1}(xi_deg); %Polynomial best fit of cCDF in-deg fpn
yi_pow = fpnresult{2}(id_emp(:,1)); %Power best fit of cCDF in-deg fpn
yi_exp = fpnresult{3}(id_emp(:,1)); %Exponential best fit of cCDF in-deg fpn
yi_gauss = fpnresult{4}(id_emp(:,1)); %Gaussian best fit of cCDF in-deg fpn
% yo_poly = fpnresult{1}(xo_deg); %Polynomial best fit of cCDF out-deg fpn
yo_pow = fpnresult{2}(od_emp(:,1)); %Power best fit of cCDF out-deg fpn
yo_exp = fpnresult{3}(od_emp(:,1)); %Exponential best fit of cCDF out-deg fpn
yo_gauss = fpnresult{4}(od_emp(:,1)); %Gaussian best fit of cCDF out-deg fpn

% get subtightplot function from:
% https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot

%setup the subtightplot function parameters to minimize spacing
make_it_tight = true; %used to turn on/off subplot functionality
%set ([vert horiz](axes gap),[lower uppper](margins),[left right](margins))
subplot = @(m,n,p) subtightplot (m, n, p, [0.035 0.01], [0.075 0.05], [0.05 0.01]);
if ~make_it_tight,  clear subplot;  end

xlimsi=[1.5,25]; xlimso=[2,25]; ylims=[0.03,1.25];

idodFig = figure(2); clf %setup parent fig
idTailfit = subplot(2,2,1); %id dist cCDF tail fit
loglog(id_emp(:,1),id_emp(:,2),'ko','MarkerSize',15,'MarkerFaceColor',[1 1 1]); hold on;
loglog(id_pow_tfit(:,1),id_pow_tfit(:,2),'r:','LineWidth',4); %power law fit
loglog(id_exp_tfit(:,1),id_exp_tfit(:,2),'b:','LineWidth',4); hold off; %exp fit
set(gca,'YLim',ylims,'XLim',xlimsi,'FontSize',20,'xticklabel',{[]}); 
text(xlimsi(1),ylims(2)+0.25,'A','FontSize',24,'FontWeight','bold')
ylabel('P(degree \geq x)','FontSize',20); 
title('In-degree Distribution cCDF','FontSize',20); %update for id/od/deg accordingly

idTailfit = subplot(2,2,2); %od dist cCDF tail fit
loglog(od_emp(:,1),od_emp(:,2),'ko','MarkerSize',15,'MarkerFaceColor',[1 1 1]); hold on;
loglog(od_pow_tfit(:,1),od_pow_tfit(:,2),'r:','LineWidth',4); %power law fit
loglog(od_exp_tfit(:,1),od_exp_tfit(:,2),'b:','LineWidth',4); hold off; %exp fit
set(gca,'YLim',ylims,'XLim',xlimso,'FontSize',20,'xticklabel',{[]},...
    'yticklabel',{[]}); 
text(xlimso(1),ylims(2)+0.25,'B','FontSize',24,'FontWeight','bold')
title('Out-degree Distribution cCDF','FontSize',20); %update for id/od/deg accordingly

idFullfit = subplot(2,2,3); %id dist cCDF full fit
loglog(id_emp(:,1),id_emp(:,2),'ko','MarkerSize',15,'MarkerFaceColor',[1 1 1]); hold on;
loglog(id_emp(:,1),yi_exp,'b:','LineWidth',4); %exp fit
loglog(id_emp(:,1),yi_pow,'r:','LineWidth',4); %power law fit
loglog(id_emp(:,1),yi_gauss,'m:','LineWidth',4); hold off; %gaussian fit
xticks([min(id_emp(:,1)) 10 max(id_emp(:,1))])
set(gca,'XLim',xlimsi,'YLim',ylims,'FontSize',20);
text(xlimsi(1),ylims(2)+0.25,'C','FontSize',24,'FontWeight','bold')
ylabel('P(degree \geq x)','FontSize',20); xlabel('x','FontSize',20)

odFullfit = subplot(2,2,4); %od dist cCDF full fit
loglog(od_emp(:,1),od_emp(:,2),'ko','MarkerSize',15,'MarkerFaceColor',[1 1 1]); hold on;
loglog(od_emp(:,1),yo_exp,'b:','LineWidth',4); %exp fit
loglog(od_emp(:,1),yo_pow,'r:','LineWidth',4); %power law fit
loglog(od_emp(:,1),yo_gauss,'m:','LineWidth',4); hold off; %gaussian fit
xticks([min(od_emp(:,1)) 10 max(od_emp(:,1))])
set(gca,'XLim',xlimso,'YLim',ylims,'FontSize',20,'yticklabel',{[]}); 
text(xlimso(1),ylims(2)+0.25,'D','FontSize',24,'FontWeight','bold')
xlabel('x','FontSize',20)

export_fig conn_subplots.eps -transparent % no background
export_fig conn_subplots.png -transparent % no background


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


%% test for power law distribution
% files located here: http://tuvalu.santafe.edu/~aaronc/powerlaws/
[alpha, xmin, L] = plfit(id','finite');
% The output 'alpha' is the maximum likelihood estimate of the scaling
% exponent, 'xmin' is the estimate of the lower bound of the power-law
% behavior, and L is the log-likelihood of the data x>=xmin under the
% fitted power law.
[p, gof] = plpva(id',xmin,'reps',1000); %p=0.0158, reject power law hypo for id deg, not drawn from a power-law dist

[alpha, xmin, L] = plfit(od','finite');
[p, gof] = plpva(od',xmin,'reps',5000); %p=0.1716, >0.1 so power law is plausible hypo for data

% now go to python for curve fitting comparisons
% find in-degree is better approximated by an exponential than a power law

% find out-degree

%visualize od dist along with fitted power-law dist on log-log axes
% h = plplot(od,xmin,alpha); 

%% visualize id+od dist along with fitted power-law dist on log-log axes
% 99% of this taken directly from plplot
% 1.update x with id/od & make sure you have min from prev section
% 2.comment/uncomment xi/xo & yi/yo accordingly
% 3.if not FPN AM, comment out exp fit fig elements: xi/xo & yi/yo, if
% doing this, don't need fig anyway
% reshape input vector
x = reshape(od,numel(od),1); %state whether looking at in/out/total deg
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
%         xi= c(9:end,1); %start from xmin (12)
%         yi= [1 0.790343 0.624642 0.493681 0.390177 0.308374 0.243721 ...
%             0.192623 0.120321]; %exp fit from python from xmin (12) to 21
        xo= c(9:end,1); %start from xmin (10)
        yo= [1 0.821421 0.674732 0.554239 0.455264 0.373963 0.307181 ...
            0.252325 0.207265 0.170252 0.139848 0.114874 0.0775095]; 

        figure;
        h(1) = loglog(c(:,1),c(:,2),'ko','MarkerSize',8,'MarkerFaceColor',[1 1 1]); hold on;
        h(2) = loglog(cf(:,1),cf(:,2),'r--','LineWidth',2); %power law fit
%         h(3) = loglog(xi,yi,'b:','LineWidth',2); hold off; %exp fit
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
%         xi= c(9:end,1); %start from xmin (12)
%         yi= [1 0.790343 0.624642 0.493681 0.390177 0.308374 0.243721 ...
%             0.192623 0.120321]; %exp fit from python from xmin (12) to 21
%         exp fit from python from xmin (10) to 23
        xo= c(6:end,1); %start from xmin (10) 
        yo= [1 0.821421 0.674732 0.554239 0.455264 0.373963 0.307181 ...
            0.252325 0.207265 0.170252 0.139848 0.114874 0.0775095]; 


        figure;
        h(1) = loglog(c(:,1),c(:,2),'ko','MarkerSize',15,'MarkerFaceColor',[1 1 1]); hold on;
        h(2) = loglog(cf(:,1),cf(:,2),'r--','LineWidth',4); %power law fit
%         h(3) = loglog(xi,yi,'b:','LineWidth',4); hold off; %exp fit
%         h(3) = loglog(xo,yo,'b:','LineWidth',4); hold off; %exp fit
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
set(gca,'YLim',[0.03,1.25],'YTick',10.^yrt,'FontSize',20);
ylabel('P(degree \geq x)','FontSize',22); xlabel('x','FontSize',22)
title('Out-degree Distribution cCDF','FontSize',20); %update for id/od/deg accordingly
% legend('empirical data','power-law fit','exponential fit')
% export_fig id_ccdf.eps -transparent % no background
% export_fig id_ccdf.png -transparent % no background
id_emp=c; id_pow_tfit=cf; id_exp_tfit=[xi,yi']; %save id tail fits, then redo for od
od_emp=c; od_pow_tfit=cf; od_exp_tfit=[xo,yo']; %save od tail fits
save('tail_fits.mat','id_emp','od_emp','id_pow_tfit','od_pow_tfit',...
    'id_exp_tfit','od_exp_tfit') %save all the tail fits


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

%init variables
CM = AM; %define connectivity matrix
[f,F]=motif3struct_bin(CM); % real data

% %be sure not to include the diagonal in the mean distance calculation
% rowMean = sum(D,2) ./ sum(D~=0,2); %avg shortest path length
% char_path_length = mean(rowMean); %characteristic path length

%random Maslov-Sneppen re-wiring
% ..the number of iterations should exceed at least 100 times the number of 
% connections in the network (Milo et al., 2004). 100*399 = 39,900 (from
% Bullmore book, but Milo et al, 2004 states 100 is more than adequate and
% 10 appears "adequate" according to its fig 1)
% iter = 100*sum(id); % # of iterations
iter = 10*sum(id); % # of iterations
% networks = 100000; % # of surrogate networks 100 for small world, 1000
% for p-val stats, 100000 for saved null_networks-motifs.mat file
networks = 1000; % # of surrogate networks 100 for small world, 1000 for p-val stats
ensemble = zeros(size(CM,1),size(CM,2),networks); %init ensemble
C_ensemble = zeros(networks,1); %init clust coef ensemble
L_ensemble = zeros(networks,1); %init char path length ensemble
f_ensemble = zeros(13,networks); %init network motif freq fingerprint ensemble
F_ensemble = zeros(13,30,networks); %init node motif freq fingerprint ensemble
p_vals_rand_f = zeros(13,1); %init p-value 
p_vals_rand_F = zeros(13,30); %init p-value 
% doc randmio_dir
tic
parfor i=1:networks
    % R: randomized network, eff: number of actual rewirings carried out
    [R,eff] = randmio_dir(CM, iter); 
    ensemble(:,:,i) = R; %build ensemble of surrogate networks
    C_ensemble(i) = mean(clustering_coef_bd(R)); 
    D_rand = distance_bin(R); %shortest path length for each node
    [L_ensemble(i),~] = charpath(D_rand); 
    [f_ensemble(:,i),F_ensemble(:,:,i)]=motif3struct_bin(R); %network motif freq fingerprint
    p_vals_rand_f = p_vals_rand_f + (f_ensemble(:,i) > f);
    p_vals_rand_F = p_vals_rand_F + (F_ensemble(:,:,i) > F);
end
toc
% change from total times freq count in null networks > empirical to 
% fraction, making it a p-val
p_vals_rand_f = p_vals_rand_f ./ networks; 

% 5991.728168 seconds = 100 surrogate networks (office pc) 39,990 iter
% 600.054129 seconds = 100 surrogate networks (office pc) 3,9990 iter
% 4081.617863 seconds = 100 surr networks (koko-CogNeuroLab) 39,990 iter
% 409.284658 seconds = 100 surr networks (koko-CogNeuroLab) 3,9990 iter
% 4085.557941 seconds = 1000 surr networks (koko-CogNeuroLab) 3,9990 iter
% [id_test,od_test,deg_test] = degrees_dir(ensemble(:,:,57)); %test deg dist
% id==id_test
% od==od_test
% all deg dist are maintained



% while (p_vals_latt_f(9) > 38 || p_vals_latt_f(9)==0 || p_vals_latt_f(13) > 41 || p_vals_latt_F(9,13) > 50)
%lattice
%Rlatt,  latticized network in original node ordering
%Rrp, latticized network in node ordering used for latticization
%ind_rp, node ordering used for latticization
%eff, number of actual rewirings carried out
latt_ensemble = zeros(size(CM,1),size(CM,2),networks); %init ensemble
C_latt_ensemble = zeros(networks,1); %init clust coef ensemble
f_latt_ensemble = zeros(13,networks); %init network motif freq fingerprint ensemble
F_latt_ensemble = zeros(13,30,networks); %init node motif freq fingerprint ensemble
p_vals_latt_f = zeros(13,1); %init p-value 
p_vals_latt_F = zeros(13,30); %init p-value 
tic
parfor i=1:networks
    % R: randomized network, eff: number of actual rewirings carried out
%     [Rlatt,Rrp,ind_rp,eff] = latmio_dir_connected(AM, iter); 
    [Rlatt,Rrp,ind_rp,eff] = latmio_dir(CM, iter); %not fully connected
    latt_ensemble(:,:,i) = Rrp; %build ensemble of surrogate networks
    C_latt_ensemble(i) = mean(clustering_coef_bd(Rrp)); 
    D_latt = distance_bin(Rrp); %shortest path length for each node
    [f_latt_ensemble(:,i),F_latt_ensemble(:,:,i)]=motif3struct_bin(Rrp); %network motif freq fingerprint
    p_vals_latt_f = p_vals_latt_f + (f_latt_ensemble(:,i) > f);
    p_vals_latt_F = p_vals_latt_F + (F_latt_ensemble(:,:,i) > F);
end
toc

%25 seconds = 100 surrogate networks
% 1834.272242 seconds = 100 surr networks, 3,990 iter (office pc)
% 1193.954145 seconds = 100 surr networks, 3,990 iter (koko)
% 12008.572344 seconds = 1000 surr networks, 3,990 iter (koko)
% change from total times freq count in null networks > empirical to 
% fraction, making it a p-val, 3, 5, 9 & 13 sig
p_vals_latt_f = p_vals_latt_f ./ networks; 

%% visualize motif frequency spectra
load('null_networks-motif_and_small_world.mat') %fully connected lattice
load('null_networks-motifs.mat') %not fully connected lattice

% 
% f_rand_ci = zeros(size(f_ensemble,1),2); %init rand 95% ci of mu stat describing motif network fingerprints
% f_latt_ci = zeros(size(f_ensemble,1),2); %init latt 95% ci of mu stat describing motif network fingerprints
% 
% for i=1:size(f_ensemble,1)
%     %each class ID's ensemble is normally distributed, so can fit this dist
%     pd = fitdist(f_ensemble(i,:)','Normal'); %fit normal dist to each class ID
%     ci = paramci(pd); f_rand_ci(i,:) = ci(:,1); %compute & save 95% ci of mu statistic
%     pd = fitdist(f_latt_ensemble(i,:)','Normal'); ci = paramci(pd,'Alpha',0.01); 
%     f_latt_ci(i,:) = ci(:,1); %doing the same for latt as above for rand
% end
% 
% err_low = [mean(f_ensemble,2)-f_rand_ci(:,1) ...
%     mean(f_latt_ensemble,2)-f_latt_ci(:,1)]; %length below the mean null fingerprints
% err_high = [f_rand_ci(:,2)-mean(f_ensemble,2) ...
%     f_latt_ci(:,2)-mean(f_latt_ensemble,2)]; %length below the mean null fingerprints

%confidence intervals are so tight around mean it doesn't make sense to plot them

figure(1), clf
x = (1:13);
vals = [f'; mean(f_ensemble,2)'; mean(f_latt_ensemble,2)'];
b = bar(x,vals); b(1).FaceColor = 'flat'; b(2).FaceColor = 'flat'; 
b(3).FaceColor = 'flat'; hold on
b(1).CData = 1/255*[204 255 153]; %empirical motif fingerprint
b(2).CData = 1/255*[128 128 128]; %random mean motif fingerprint
b(3).CData = 1/255*[192 192 192]; %lattice mean motif fingerprint
hold on
xBar=cell2mat(get(b,'XData')).' + [b.XOffset];  % compute bar centers
% apply error bars to only random and lattice plots
% hEB=errorbar(xBar(:,2:3),vals(2:3,:)',err_low,err_high,'k','LineStyle','none'); 
% yl=ylim; ylim([0,yl(2)]); % bring back to 0
% set(hEB,'LineStyle','none');
%draw significance lines over class ID 9 & 13, first set x val's
sigIDs_x = [b(1).XEndPoints(9) b(2).XEndPoints(9) b(1).XEndPoints(9) b(3).XEndPoints(9) ...
    b(1).XEndPoints(13) b(2).XEndPoints(13) b(1).XEndPoints(13) b(3).XEndPoints(13)];
sigIDs_y = [f(9)+20 f(9)+20 f(9)+40 f(9)+40 f(13)+20 f(13)+20 f(13)+40 f(13)+40];%now set corresponding y val's
% line(sigIDs_x,sigIDs_y) %no good, lines connect the lines
plot(sigIDs_x(1:2), sigIDs_y(1:2), '-k', 'LineWidth',2) %sig ID 9 emp/rand
plot(sigIDs_x(3:4), sigIDs_y(3:4), '-k', 'LineWidth',2) %sig ID 9 emp/latt
% plot(sigIDs_x(5:6), sigIDs_y(5:6), '-k', 'LineWidth',2) %sig ID 13 emp/rand
% plot(sigIDs_x(7:8), sigIDs_y(7:8), '-k', 'LineWidth',2) %sig ID 13 emp/latt
set(gca,'TickLength',[0 0],'FontSize',18) %remove ticks
sig_rand = [xBar(9,1)+0.03 sigIDs_y(1)+5; xBar(9,1)+0.10 sigIDs_y(1)+5; ...
    xBar(9,1)+0.17 sigIDs_y(1)+5; xBar(13,1)+0.03 sigIDs_y(5)+5; ...
    xBar(13,1)+0.10 sigIDs_y(5)+5; xBar(13,1)+0.17 sigIDs_y(5)+5;]; %*** p<.001
sig_latt = [xBar(9,2) sigIDs_y(3)+5; xBar(13,2) sigIDs_y(7)+5]; %* p<=.05
plot(sig_rand(1:3,1),sig_rand(1:3,2), '*k') %*** p<.001
plot(sig_latt(1,1),sig_latt(1,2), '*k') %* p<=.05
ylabel('structural motif count','FontSize',21); xlabel('motif ID (M=3)','FontSize',21)
% legend('Real','Random','Lattice'); title('Motif Frequency Spectra','FontSize',20)
h=gca; h.Color = 'none'; % turn off background color
box off; %take out top and right lines

export_fig motif_spectra.eps -transparent % no background
export_fig motif_spectra.png -transparent % no background

%% Motif analysis
%first must make_motif34lib.m which generates the required motif34lib.mat
%which all motif releated functions require
make_motif34lib

M = 3; %motif size

%next be sure to generate null hypothesis distributions through random and
%lattice-based network ensemble generation

[f,F]=motif3struct_bin(AM); % real data
% The motif frequency of occurrence around an individual node is known as
% the motif fingerprint of that node. The total motif frequency of
% occurrence in the entire network is known as the motif fingerprint of the network. 
%F is node motif frequency fingerprint
%f is network motif frequency fingerprint
f_rand_mean = mean(f_ensemble,2); %compute the mean motif fingerprint over all randomly generated networks
f_rand_std = std(f_ensemble,0,2); %compute the std dev over all rand gen networks
f_rand_z = (f - f_rand_mean) > (0.1 * f_rand_mean); % Milo, 2002 method, only 9 & 13 sig
f_rand_z = (f - f_rand_mean) ./ f_rand_std; 

F_rand_z = (F - mean(F_ensemble,3)) ./ std(F_ensemble,0,3); %compute nodal z-score over all rand networks

f_latt_mean = mean(f_latt_ensemble,2); %compute the mean motif fingerprint over all lattice networks
f_latt_std = std(f_latt_ensemble,0,2); %compute the std dev over all lattice networks
f_latt_z = (f - f_latt_mean) > (0.1 * f_latt_mean); % Milo, 2002 method, 3, 5, & 9 sig
f_latt_z = (f - f_latt_mean) ./ f_latt_std; %   sig

F_latt_z = (F - mean(F_latt_ensemble,3)) ./ std(F_latt_ensemble,0,3); %compute nodal z-score over all latt networks

p_vals_latt_f = p_vals_latt_f ./ networks; p_vals_latt_F = p_vals_latt_F ./ networks; 
p_vals_rand_f = p_vals_rand_f ./ networks; p_vals_rand_F = p_vals_rand_F ./ networks; 

%interweave z-scores for easier copy-pasting to excel for table
F_z = zeros(2*size(F_latt_z,1),size(F_latt_z,2)); %embedding matrix
F_z(1:2:end) = F_rand_z; F_z(2:2:end) = F_latt_z;
f_latt_z'

%% combinatorics to find which nodes are a part of each class ID
C = nchoosek(1:length(AM),M); 
% yields 4060x3 matrix which provides all possible 3-node combinations in 
% the network
C_bd = zeros(M,M,length(C)); %init mat for all possible 3-node bin+dir connectivity schemes in emp FPN
C_M9 = zeros(vals(1,9),1); j=1; %init mat & counter (j) for all M9 isomorphs in emp FPN
C_M13 = zeros(vals(1,13),1); k=1; %init mat & counter (k) for all M13 isomorphs in emp FPN
for i=1:length(C)
    % 1. pull out the connectivity for each node combo yields 3x3 matrix 
    % for each row in C, (save as 3x3x4060 matrix C_bd)
    C_bd(:,:,i) = [AM(C(i,:),C(i,:))]; %binary, directed connectivity for i'th row of C
    % 2. compare with sig motif isomorphs from find_motif34 & save in mat
    if find_motif34(C_bd(:,:,i)) == 9 %is it of class ID 9?
        C_M9(j) = i; j=j+1; %store idx of C
    elseif find_motif34(C_bd(:,:,i)) == 13 %or of class ID 13?
        C_M13(k) = i; k=k+1; %store idx of C
    end
end

%13, 14 & 28 for M9
%6, 7, 13 & 14 for M13 (only frontal areas)
nchoosek(4,3) %4 combinations:
% 6, 7, 13
% 6, 7, 14
% 6, 13, 14
% 7, 13, 14

%% Small-world analysis
%run surrogate networks first to generate ensembles

CM = CIJ; %define connectivity matrix

%Humphries et al. (2006& 2008) method of calculating sigma:
%only uses random networks
C = mean(clustering_coef_bd(CM)); % avg clustering coef for empirical network
D = distance_bin(CM); [L,efficiency] = charpath(D); %L = char path length
gamma = C / mean(C_ensemble); %Gamma > 1 suggests greater clustering than random
lambda = L / mean(L_ensemble); %Lambda ~ 1 suggests a comparable avg path length to randomized network
sigma = gamma / lambda; %sigma > 1 indicates small-worldness (1.1429)

%Telesford et al., 2011 technique:
%uses both random and latticized networks to capture both sides of the
%spectrum
omega = (mean(L_ensemble) / L) - (C / mean(C_latt_ensemble));
%omega = 0.0275 (1,000 networks), 0.0297 (100,000 networks)
%omega index ranges between -1 and 1. Values bet [-0.5,0.5] indicate 
%small-worldness. Positive values suggest more random characteristics and 
%negative values indicate a lattice-like structure. 

%% Modeling analysis to see if Gaussian deg dist leads to a particular kind of motif overrepresentation

%FPN has 30 nodes and 399 directed edges
N = 30; %number of nodes in network
K = 399; %number of directed edges in network
s = 1; %std of toeplitz

%test
N = 30; %number of nodes in network
K = 100; %number of directed edges in network
s = 1; %std of toeplitz

%K=100 was 0.024881 seconds, K=200 was untenable

tic
CIJ = maketoeplitzCIJ(N,K,s);
toc

%identify in-degree, out-degree and in+out=total degree per node
[id,od,deg] = degrees_dir(CIJ);

%step 1: run section %% test for power law distribution
%step 2: run section %% visualize id+od dist along with fitted power-law dist on log-log axes
%step 3: run section %% testing best cCDF plot
%step 4: run section %% surrogate networks
%step 5: run section %% Small-world analysis
%this toeplitz matrix is not small-world: sigma=2.697 & omega=-0.7774
%indeed visualization shows lattice network
%however it is indeed best approximated by a Gaussian

%Re-run this analysis and re-create the association matrix image of toeplitz
%Need to re-write e-mail to Will asking his opinion.
%next, try to figure out how to generate a small-world network with a
%gaussian in & out degree distribution.