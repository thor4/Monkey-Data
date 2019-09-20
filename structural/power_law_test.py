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
#results.power_law.plot_pdf(color='b', linestyle='--', ax=fig2) # power law fit
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


### IN & OUT-DEGREE DISTRIBUTIONS OF MONKEY FPN ###
import numpy as np
import powerlaw
import matplotlib.pyplot as plt

## Import data ##
#bct docs: http://htmlpreview.github.io/?https://github.com/aestrivex/bctpy/blob/master/function_reference.html
#figure out how to load adjacency matrix in here from .mat file then do id + od
id = [13 17 14 9 6 19 21 19 15 15 15 11 21 18 12 12 19 2 6 4 7 10 8 16 11 17 .. 
      16 17 15 14]