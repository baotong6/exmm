'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-11-03 09:48:24
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2023-11-03 10:53:45
FilePath: /pulsar/XMMcentral/read_evt.py
Description: 

Copyright (c) 2023 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from stingray.pulse.search import epoch_folding_search, z_n_search
from stingray.pulse.search import search_best_peaks
from stingray.stats import fold_detection_level, z2_n_detection_level
import astropy.units as u
import astropy.constants as c
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }
path='/Users/baotong/Desktop/XMMcentral/'
evtfile=fits.open(path+'src_events1.fits')
evtdata=evtfile[1].data
time=evtdata['time'];energy=evtdata['PI']

nharm = 1
df_min=1e-4
frequencies=np.arange(1/10,1/0.1,df_min)
nbin=20
ntrial = (frequencies[-1] - frequencies[0]) / df_min
freq, zstat = z_n_search(time, frequencies, nbin=nbin, nharm=nharm)
z_detlev = z2_n_detection_level(n=1, epsilon=0.001, ntrial=len(freq)*1)
z_detlev2 = z2_n_detection_level(n=1, epsilon=0.01, ntrial=len(freq)*1)
z_detlev3 = z2_n_detection_level(n=1, epsilon=0.1, ntrial=len(freq)*1)
# ---- PLOTTING --------
plt.figure(figsize=(10,8))
plt.plot(freq, (zstat - nharm), label='$Z_2$ statistics')
plt.axhline(z_detlev - nharm, label='conf = 99.9%',color='r',linestyle='--')
plt.axhline(z_detlev2 - nharm, label='conf = 99%',color='g',linestyle='--')
plt.axhline(z_detlev3 - nharm, label='conf = 90%',color='b',linestyle='--')
plt.text(freq[np.argmax(zstat)],zstat[np.argmax(zstat)]-0.6,
         'P={:.4f} s'.format(1/freq[np.argmax(zstat)]),font1)
# plt.plot(freq, efstat - nbin + 1, color='gray', label='EF statistics', alpha=0.5)
print('Period=',1/freq[zstat.argmax()])
# plt.axvline(1/period, linestyle='--',color='r', lw=1, alpha=0.5, label='Correct frequency')
plt.xlim([frequencies[0], frequencies[-1]])
plt.xlabel('Frequency (Hz)',font1)
plt.ylabel(r'$Z^2$ power',font1)
plt.tick_params(labelsize=16)
plt.legend()
plt.semilogx()
plt.savefig(path+'Z2_src_events1.pdf',bbox_inches='tight', pad_inches=0.01)
plt.show()
