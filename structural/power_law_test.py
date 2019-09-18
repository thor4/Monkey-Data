# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
# x = (1-rand(10000,1)).^(-1/(2.5-1));
# where rand is uniformly randomly dist (0,1)
s = 1-np.random.uniform(0,1,1000) # mimic Clauset's matlab demo
data = np.array(np.power(s, (-1/(3.5-1)))) # power distribution

import powerlaw
results = powerlaw.Fit(data)
print(results.power_law.alpha)
print(results.power_law.xmin) # optimal min x value the scaling relationship of power law begins
print(results.power_law.D) # minimal Kolmogorov-Smirnov distance bet data & fit
print(results.power_law.sigma) # std error

R, p = results.distribution_compare('power_law', 'lognormal')


import matplotlib.pyplot as plt
plt.figure(0)
powerlaw.plot_pdf(data, color='b')
#plt.clf() # clears the entire figure
plt.figure(1)
fig2 = results.plot_pdf(color='b', linewidth=2) # pdf of original data
results.power_law.plot_pdf(color='b', linestyle='--', ax=fig2) # power law fit
results.plot_ccdf(color='r', linewidth=2, ax=fig2) # ccdf of original data
results.power_law.plot_ccdf(color='r', linestyle='--', ax=fig2) # power law fit


x,y = results.cdf() # sorted original data (x) with assoc probabilities ( y, p(x<k) )
bin_edges, probability = results.pdf() # left & right edges (bin_edges) of pdf
# small spacing to large (log-scale), with assoc prob (probability, p(x=k) )
y = results.lognormal.cdf(data=[2,3,4]) # gives assoc prob p(x<k) for each datum
y = results.lognormal.pdf() # gives assoc prob p(x=k) for all data