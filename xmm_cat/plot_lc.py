'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-10-31 08:49:38
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-02-21 16:56:50
FilePath: /pulsar/XMMcentral/plot_lc.py
Description: 

Copyright (c) 2023 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import poisson_conf_interval
from scipy.stats import norm
from stingray.events import EventList
from astropy.timeseries import LombScargle
from stingray.lightcurve import Lightcurve
import os
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

def plot_LS_sig(path,dataname,outpath=None,outname=None,save=0,show=1):
     a = fits.open(path + dataname)
     rate = a[1].data['RATE']
     time = a[1].data['TIME']
     T = time[-1] - time[0]
     lc = Lightcurve(time=time, counts=rate * (time[1] - time[0]))
     freq = np.arange(1 / T, 1 / 55, 1e-6)
     x = lc.time;y = lc.counts
     LS = LombScargle(x, y, normalization='standard')
     power = LS.power(freq)
     max_NormLSP = np.max(power)
     period_peak = 1. / freq[np.where(power == np.max(power))][0]
     FP = LS.false_alarm_probability(power.max(), minimum_frequency=freq[0], maximum_frequency=freq[-1],
                                     method='baluev')
     FP_99 = LS.false_alarm_level(0.0027, minimum_frequency=freq[0], maximum_frequency=freq[-1], method='baluev')
     sim_LSP=np.loadtxt(path+dataname[:-3]+'.txt')
     sim_LSP_3sig=np.sort(sim_LSP)[-3]
     sim_LSP_2sig = np.sort(sim_LSP)[-50]
     plt.figure(1, (10, 8))
     plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '-',label=r'$3 \sigma~(baluev)$')
     # plt.plot([freq[0], freq[-1]], [max_NormLSP, max_NormLSP], '--',label='{0:.2f} '.format(sigma)+r'$\sigma$')
     plt.plot([freq[0], freq[-1]], [sim_LSP_3sig, sim_LSP_3sig], '-',label=r'$3 \sigma~(simulated)$')
     plt.plot([freq[0], freq[-1]], [sim_LSP_2sig, sim_LSP_2sig], '-', label=r'$2 \sigma~(simulated)$')
     # plt.title('Period={0:.2f}'.format(period_peak), font1)
     plt.text(freq[np.where(power == np.max(power))][0] * 1.3, max_NormLSP * 0.95, 'P={0:.2f}s'.format(period_peak),
              fontsize=18, fontweight='semibold')
     plt.plot(freq, power)
     plt.semilogx()
     # plt.semilogy()
     out_period = 1. / freq[np.where(power == np.max(power))][0]
     plt.legend()
     plt.xlabel('Frequency (Hz)', font1)
     plt.ylabel('Normalized LS Periodogram',font1)
     plt.tick_params(labelsize=16)
     if save:
         plt.savefig(outpath + outname + '_LS.pdf', bbox_inches='tight', pad_inches=0.01)
     if show:
         plt.show()
     else:
         plt.close()
     return [FP, out_period, max_NormLSP]

def plot_lc(path,dataname,outpath=None,outname=None,save=0,show=1):
    a = fits.open(path + dataname)
    rate = a[1].data['RATE']
    time = a[1].data['TIME']
    T = time[-1] - time[0]
    lc = Lightcurve(time=time, counts=rate * (time[1] - time[0]))
    x = lc.time;y = lc.counts
    plt.plot(x,y)
    plt.show()

if __name__=='__main__':
    path = '/Users/baotong/Desktop/XMMcentral/all_lc/'
    file_names = os.listdir(path);list=[]
    for file_name in file_names:
        if file_name.endswith(".lc"):
            list.append(file_name)
    file_name='0886081101_268.67_26.99_903_pnsrc_2_10keV.lc'
    plot_lc(path,file_name,outpath=path,outname=file_name[:-3],show=1,save=1)
    # plot_LS_sig(path,file_name,outpath=path,outname=file_name[:-3],show=1,save=1)





