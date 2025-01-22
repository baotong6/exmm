'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-12-07 13:37:01
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2023-12-11 13:33:40
FilePath: /pulsar/XMMcentral/read_nustar_lc.py
Description: 

Copyright (c) 2023 by baotong, All Rights Reserved. 
'''

import numpy as np
import os
import matplotlib.pyplot as plt
import glob
import subprocess
import csv
from stingray.events import EventList
import pandas as pd
from astropy.stats import poisson_conf_interval
import math
import stingray as sr
from stingray.lightcurve import Lightcurve
from astropy.io import fits
import rednoise as rednoise
import hawkeye as hawk

path='/Users/baotong/Desktop/XMMcentral/nustar_xmm/'
file_name='nu80801332002A01_sr_3_10keV.lc'
lcfile=fits.open(path+file_name)
GTI=lcfile[2].data
GTI=np.array(GTI)
GTI_array = np.array([list(x) for x in GTI.tolist()])
t0=GTI_array[:,0];t1=GTI_array[:,1]

# for i in range(len(t0)):
#     plt.plot([t0[i],t1[i]],[1,1],'-')
# plt.show()
gap=[];GTI_use=[]
for i in range(len(t0)-1):
    if (t0[i+1]-t1[i])>1000:
        gap.append([t1[i],t0[i+1]])
for j in range(len(gap)):
    if j==0:
        GTI_use.append([t0[0],gap[0][0]])
    else:
        GTI_use.append([gap[j-1][1],gap[j][0]])
GTI_use.append([gap[-1][1],t1[-1]])
print(len(gap),len(GTI_use))
rate = lcfile[1].data['RATE']
time = lcfile[1].data['TIME']
tstart=lcfile[1].header['TSTART']
time=time+tstart
mask = np.zeros_like(time, dtype=bool)
max_LSP=[];lc_all=[]
for i in range(30):
    np.random.seed(i+200)
    filtered_time_intervals = []
    t_sim_all=[]
    for interval in GTI_use:
        interval_mask = (time >= interval[0]) & (time <= interval[1])
        filtered_time_intervals.append(time[interval_mask])
        lc_cho = Lightcurve(time=time[interval_mask], counts=rate[interval_mask] * (time[1] - time[0]))
        N_cts=int(np.sum(lc_cho.counts));bin=lc_cho.time[1]-lc_cho.time[0]
        T=time[interval_mask][-1]-time[interval_mask][0]
        t = np.random.random(np.random.poisson(N_cts)) * T+ lc_cho.time[0]
        t = np.sort(t)
        evt1 = EventList()
        evt1.time=t
        lc_single=evt1.to_lc(dt=time[1]-time[0])
        lc_all.append(lc_single)
        t_sim_all.extend(t)
    for k in range(len(lc_all)):
        if k==0:
            lc_long=lc_all[0]
        else:
            lc_long = lc_long.join(lc_all[k], skip_checks=True) 

    [FP,out_period,max_NormLSP]=hawk.get_LS(lc_long.time, lc_long.counts, freq=np.arange(1 / T, 1 / 55, 1e-6),
                                                outpath=None, outname=None, save=0, show=0)
    with open(path + file_name[:-3] + '.txt', 'a+') as f:
        f.write(f"{max_NormLSP:.10f}\n")
    max_LSP.append(max_NormLSP)
# np.savetxt(path+file_name[:-3]+'.txt',max_LSP, fmt='%.10f', delimiter='\n')

# T = time[-1] - time[0]
# lc_cho = Lightcurve(time=time, counts=rate * (time[1] - time[0]))
# CR_cho = np.mean(lc_cho.counts) / lc_cho.dt
# epoch = np.array([time[0], time[-1], '11111', T])
#     # frac_rms = np.sqrt(np.var(lc_cho.counts) * lc_cho.counts.size / (lc_cho.counts.size - 1) / np.mean(lc_cho.counts) ** 2)
#     # sigma_range=poisson_conf_interval(lc_cho.counts)
#     # sigma=(sigma_range[1:,]-sigma_range[0,:])/2
#     # mse=np.mean(sigma**2)
#     # frac_rms=np.sqrt((np.var(lc_cho.counts) -mse) /np.mean(lc_cho.counts) ** 2)
#     # print('frms=',frac_rms)
#     # (psd_sim,result_mu) = bestpsd(lc_cho, epoch_info=epoch_info_use,maskfreq=maskfreq)
#     # print(np.sum(lc_cho.counts))
# N_cts=int(np.sum(lc_cho.counts));bin=lc_cho.time[1]-lc_cho.time[0]
# t = np.random.random(np.random.poisson(N_cts)) * T + time[0]
# t = np.sort(t)
# evt = EventList()
# print(len(t))
# # 创建一个布尔掩码，标记哪些元素在bad time intervals中
# mask = np.ones_like(t, dtype=bool)

# for interval in gap:
#     print(interval[0],interval[1])
#     mask &= ~((t >= interval[0]) & (t<= interval[1]))

# evt.time = t[mask]
# print(len(evt.time))
# lc_out = evt.to_lc(dt=bin)
# lc_out.plot()
# plt.show()











