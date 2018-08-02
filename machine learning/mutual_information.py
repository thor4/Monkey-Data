# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 15:00:46 2018

@author: bryan
"""

import scipy.io
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.model_selection import train_test_split
from sklearn.metrics import mutual_info_score
import matplotlib.pyplot as plt

#split 80/20 train/test sets, train_test_split shuffles by default and stratify
#to ensure biased class ratio is maintained
#X_train1, X_test1, y_train1, y_test1 = train_test_split(X1, y1, test_size=0.2, random_state=42, stratify=y1)

# the histogram of the data
n, bins, patches = plt.hist(x, 50, density=True, facecolor='g', alpha=0.75)
c_xy = np.histogram2d(x, y, bins)[0]


plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title('Histogram of IQ')
plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.axis([40, 160, 0, 0.03])
plt.grid(True)
plt.show()

#feature selection
#maximal information coefficient, standardization of training data yields the same scores
def calc_MI(x, y, bins):
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi

def getChan(monkey, day, chan='', response=''):
    monkey_ = monkey['monkey1'] #pull out monkey1 struct
    day = monkey_['day'] #pull out day object
    days = day[0,0] #pull out day struct
    correct_days = days['correct']; incorrect_days = days['incorrect']
    if not response:
        print("Please specify behavioral response")
    elif response == 'correct':
        cor = correct_days[0,day-1] #pull out correct day
        if not chan:
            print("Please specify recording chan")
        else:
            data_ = cor[chan] #pull out chan data object
            data = data_[0,0] #pull out chan data struct
    else:
        incor = incorrect_days[0,day-1] #pull out incorrect day
        if not chan:
            print("Please specify recording chan")
        else:
            data_ = cor[chan] #pull out chan data object
            data = data_[0,0] #pull out chan data struct
    return data
        

def search_photos(fl_tags, license_val=''):
    if not license_val:
        photos = flickr.photos_search(
            tags=fl_tags,
            is_commons="True")
    else:
        photos = flickr.photos_search(
            tags=fl_tags,
            license=license_val)
    return photos

search_photos("foo") # gives you the same output as before
search_photos("foo", "9") # gives you the search with license as 9

m1d1CorChan3dPFC_ = cor_day1['chan3dPFC'] 
m1d1CorChan3dPFC = m1d1CorChan3dPFC_[0,0]; del m1d1CorChan3dPFC_;
    return mi

m1 = scipy.io.loadmat('m1GoodStableRule1PingRej-split_by_Day_BehResp_and_Chan.mat')

m1d1CorChan1PEC = getChan(m1, 1, chan='chan1PEC', response='correct')

day1mi = calc_MI(m1d1CorChan1PEC[:,0],m1d1CorChan3dPFC[:,0],5)

# https://stackoverflow.com/questions/20491028/optimal-way-to-compute-pairwise-mutual-information-using-numpy

scores = selector.scores_
#extract & plot which of the 250ms timepoints are most important
df = pd.DataFrame(scores);
top250 = df.nlargest(250,0);
topidx = np.sort(top250.index.values);
plt.figure()
plt.hist(topidx, bins=[0,100,200,300,400,500,600,700,810], facecolor='green', alpha=0.75)
plt.xlabel('Time (ms)')
plt.title('Histogram of Mutual Information Score for Each Time Point')
plt.grid(True)
plt.axis('tight')