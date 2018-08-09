# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 10:43:40 2018

@author: bryan
"""

import h5py
import numpy as np
from itertools import *

f = h5py.File('mGoodStableRule1PingRej.h5','r')

f.keys()
f.values()
f.items()

f.name
for name in f:
    print(name)
    
def find_chan(name):
    #find chan in directory
    if 'LIP' in name:
        return name

def print_attrs(name, obj):
    print(name)
    for key, val in obj.attrs.iteritems():
        print("    %s: %s" % (key, val))

g = f.visit(get_names)

f.visititems(print_attrs)
    
d = {'a': 1, 'b': 2}
dvi = f.itervalues()
dvi.next()


g = f['monkey/day']

data = f.get('monkey/day')

type(g) #dataset

data.keys()

data.shape
data.size
data.dtype #object

data[:]

monkey1 = data.get[1,0]


