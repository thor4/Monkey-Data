# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
# x = (1-rand(10000,1)).^(-1/(2.5-1));
# where rand is uniformly randomly dist (0,1)
s = 1-np.random.uniform(0,1,1000) # mimic Clauset's matlab demo
x = np.power(s, (-1/(2.5-1))) # power distribution
