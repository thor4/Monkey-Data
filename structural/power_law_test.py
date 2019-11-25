# -*- coding: utf-8 -*-
"""
Fitting power law distributions to data using the Powerlaw package from
https://github.com/jeffalstott/powerlaw
"""
####  DEMO  ###
### Generate Data ##
#import numpy as np
## x = (1-rand(10000,1)).^(-1/(2.5-1));
## where rand is uniformly randomly dist (0,1)
#s = 1-np.random.uniform(0,1,1000) # mimic Clauset's matlab demo
#data = np.array(np.power(s, (-1/(3.5-1)))) # power distribution
#
### Generate fit object ##
#import powerlaw
#results = powerlaw.Fit(data)
#print(results.power_law.alpha)
#print(results.power_law.xmin) # optimal min x value the scaling relationship of power law begins
#print(results.power_law.D) # minimal Kolmogorov-Smirnov distance bet data & fit
#print(results.power_law.sigma) # std error
#
### Generate visualizations of data and fits ##
#import matplotlib.pyplot as plt
#plt.figure(0)
#powerlaw.plot_pdf(data, color='b')
##plt.clf() # clears the entire figure
#plt.figure(1)
#fig2 = results.plot_pdf(color='b', linewidth=2) # pdf of original data
#results.power_law.plot_pdf(color='b', linestyle='--', ax=fig2) # ideal power law fit
#results.plot_ccdf(color='r', linewidth=2, ax=fig2) # ccdf of original data
#results.power_law.plot_ccdf(color='r', linestyle='--', ax=fig2) # power law fit
#
#
### Exploring accessible probabilities ##
#x,y = results.cdf() # sorted original data (x) with assoc probabilities ( y, p(x<k) )
#bin_edges, probability = results.pdf() # left & right edges (bin_edges) of pdf
## small spacing to large (log-scale), with assoc prob (probability, p(x=k) )
#y = results.lognormal.cdf(data=[2,3,4]) # gives assoc prob p(x<k) for each datum
#y = results.lognormal.pdf() # gives assoc prob p(x=k) for all data
#
#
### Comparing candidate distributions ##
#R, p = results.distribution_compare('power_law','exponential', normalized_ratio=True)
## R is log-likelihood ratio of 1st-to-2nd, positive if data is more likely in 
## first dist, neg otherwise. R normalized by its std, p-val is calculated with
## normalized R. exp dist is abs min alt candidate since def of "heavy tail" is
## that it's not exp bounded (Asmussen Sr 2003). if sig diff, then proceed to 
## motivated hypo concerning generative mechanisms to determine which dist to 
## test next.
#results.distribution_compare('power_law','truncated_power_law', normalized_ratio=True)
#results.distribution_compare('power_law','lognormal', normalized_ratio=True)
## likely that neither heavy-tailed dist is a sig stronger fit (p>0.5). can 
## only conclude moderate support for power law, without ruling out possibility
## of other dist

#clear variable workspace
%reset
### IN & OUT-DEGREE DISTRIBUTIONS OF MONKEY FPN ###
import numpy as np
import powerlaw
import matplotlib.pyplot as plt
import scipy.io
import bct

#load adjacency matrix
AM_ = scipy.io.loadmat('AM.mat')
AM = AM_['AM'] # pull out AM
AM.shape # 30x30
AM[0:3,0:3] # starts with 0 and doesn't include last digit (3)


## Import data ##
#bct docs: http://htmlpreview.github.io/?https://github.com/aestrivex/bctpy/blob/master/function_reference.html
#figure out how do id + od
#id = [13 17 14 9 6 19 21 19 15 15 15 11 21 18 12 12 19 2 6 4 7 10 8 16 11 17 16 17 15 14]
id, od, deg = bct.degrees_dir(AM) # compute in-deg, out-deg and total deg

#results = powerlaw.Fit(id) # auto finds xmin=15 which is too restrictive
#results.xmin, results.D #xmin=15, D=0.196435 (MATLAB gives 0.3739)
results = powerlaw.Fit(id,xmin=12) # auto finds xmin=15 which is too restrictive, go with MATLAB xmin=12
results.xmin, results.D #xmin=12, D=0.28646 (MATLAB is 0.2097 and smallest)
results = powerlaw.Fit(od) # create fit object for in-deg data
print(results.power_law.alpha) # id:~4.4456 (3.375 MATLAB), xmin=12 (MATLAB), od: ~3.66 (~3.25 MATLAB), xmin=10.0
print(results.power_law.xmin) # optimal min x value the scaling relationship of power law begins
print(results.power_law.D) # minimal Kolmogorov-Smirnov distance bet data & fit
print(results.power_law.sigma) # std error

#any plot from 'results' variable uses data only from x_min up, in contrast to 
#plots from 'powerlaw' class which uses all of original data
results.supported_distributions #list of supported distributions

## Generate visualizations of data and fits ##
plt.figure(0)
fig1 = powerlaw.plot_pdf(id, color='b')
#x1, x2 = fig1.get_xlim(); y1, y2 = fig1.get_ylim() # load limits for use in nxt fig
#plt.clf() # clears the entire figure
plt.figure(1)
#fig2 = results.plot_pdf(color='b', linewidth=2) # pdf of data only from x_min up
#fig2.set_xlim(x1,x2); fig2.set_ylim(y1,y2);
#results.power_law.plot_pdf(color='b', linestyle='--', ax=fig2) # ideal power law fit
fig2 = results.plot_ccdf(color='b', linewidth=2) # ccf of data only from x_min up
results.power_law.plot_ccdf(color='r', linestyle='--', ax=fig2) # power law fit
results.lognormal.plot_ccdf(ax=fig2, color='m', linestyle='-.')
results.exponential.plot_ccdf(ax=fig2, color='g', linestyle=':')
results.truncated_power_law.plot_ccdf(ax=fig2, color='y', linestyle='--')

plt.figure(2)
fig3 = powerlaw.plot_ccdf(id, color='b') #cCDF of original data
results.power_law.plot_ccdf(color='r', linestyle='--', ax=fig3) # power law fit
results.lognormal.plot_ccdf(ax=fig3, color='m', linestyle='-.')
results.exponential.plot_ccdf(ax=fig3, color='g', linestyle=':')
results.truncated_power_law.plot_ccdf(ax=fig3, color='y', linestyle='--')

## Exploring accessible probabilities ##
x,y = results.cdf() # sorted data only from x_min up (x) with assoc probabilities ( y, p(x<k) )
#bin_edges, probability = results.pdf() # left & right edges (bin_edges) of pdf
## small spacing to large (log-scale), with assoc prob (probability, p(x=k) )
#y = results.lognormal.cdf(data=[2,3,4]) # gives assoc prob p(x<k) for each datum
#y = results.lognormal.pdf() # gives assoc prob p(x=k) for all data

#y = results.exponential.ccdf() # sorted data only from x_min up (x) with assoc probabilities ( y, p(x>=k) )
y = results.exponential.ccdf(data=[12,13,14,15,16,17,18,19,21]) # gives assoc prob p(x<k) for each datum in id
y = results.exponential.ccdf(data=[10,11,12,13,14,15,16,17,18,19,20,21,23]) # gives assoc prob p(x<k) for each datum in od

## Comparing candidate distributions ##
R, p = results.distribution_compare('power_law','exponential', normalized_ratio=True)
# R is log-likelihood ratio of 1st-to-2nd, positive if data is more likely in 
# first dist, neg otherwise. R normalized by its std, p-val is calculated with
# normalized R. exp dist is abs min alt candidate since def of "heavy tail" is
# that it's not exp bounded (Asmussen Sr 2003). if sig diff, then proceed to 
# motivated hypo concerning generative mechanisms to determine which dist to 
# test next.
R,p #(R=-3.3855629932457019, p=0.00071032404153694329) id
#(R=-1.6764042975802269, p=0.093659027710146808) od
# exp is better fit than power-law for tail of in-degree distribution, but not od
results.distribution_compare('exponential', 'stretched_exponential')
# stretched exp is better than exp: (-3.3081976438013103, 0.010104414280270313) for id
#(R=-0.85293484402035791, p=0.19152213548654928) od
#neither distribution is a significantly stronger fit (p>.05) for the following:
results.distribution_compare('exponential','lognormal', normalized_ratio=True)
#(R=-1.179336415953951, p=0.23826424303514571)
results.distribution_compare('exponential', 'lognormal_positive')
#(R=-2.5779328984331871, p=0.23826424303514571)
results.distribution_compare('exponential','truncated_power_law', normalized_ratio=True)
# exp is better than trunc power law(2.4156311709062623, 0.01570796358053525)
results.distribution_compare('stretched_exponential','truncated_power_law')
#(R=3.62209319585141, p=0.11960686333194259)
results.distribution_compare('stretched_exponential','lognormal')
#(R=0.7302647453681228, p=0.14569078178385125)
results.distribution_compare('stretched_exponential','lognormal_positive')
#(R=0.7302647453681228, p=0.14569078178385125)
results.distribution_compare('power_law', 'stretched_exponential')
#(R=-5.13045806393586, p=0.05997056435316063) id
#(R=-2.6363814472541445, p=0.28622857062203577) od
results.distribution_compare('power_law','truncated_power_law', normalized_ratio=True)
#(R=-3.6929254209270748, p=0.082409492672877893) id
#(R=-1.9179450095500206, p=0.083440273525472786) od
results.distribution_compare('power_law','lognormal', normalized_ratio=True)
#(R=-1.6299391402494143, p=0.10311436030588744) id
#(R=-0.97789584568915033, p=0.32812584001695999) od
results.distribution_compare('power_law', 'lognormal_positive')
#(R=-4.4001933185677364, p=0.10311436030588744) id
#(R=-1.7979149474119933, p=0.32812584001695999) od
# negative for out-deg, stop here



# likely that neither heavy-tailed dist is a sig stronger fit (p>0.05). can 
# only conclude moderate support for power law for od, without ruling out possibility
# of other dist


##There is no support for a power law fit for either in or out-degree dist