import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import poisson_conf_interval
from scipy.stats import norm
from scipy import interpolate
from scipy.optimize import curve_fit
import astropy.units as u
import astropy.constants as c
import stingray as sr
from stingray.events import EventList
from astropy.timeseries import LombScargle
from stingray.lightcurve import Lightcurve
from stingray import Lightcurve, Crossspectrum, sampledata,Powerspectrum,AveragedPowerspectrum
from stingray.simulator import simulator, models
import rednoise
from rednoise import useful_functions as func
from rednoise import powerlw2evt as p2evt
from NGC104.timing_comb import load_data
from rednoise import Vaughan as Vaughan
import hawkeye as hawk
import os

def singleobs2simevt(path,dataname,dt=500,chosen_obs=None,ifoutobs=[],randseed=1,maskfreq=0):
    ## from one obs(usually the longest one obs) to simulate singleobs_evt
    ## if you need only just some obs to be simulated, make ifoutobs not []
    a = fits.open(path + dataname)
    rate = a[1].data['RATE']
    time = a[1].data['TIME']
    T = time[-1] - time[0]
    lc_cho = Lightcurve(time=time, counts=rate * (time[1] - time[0]))
    freq = np.arange(1 / T, 1 / 1, 1e-5)
    CR_cho = np.mean(lc_cho.counts) / lc_cho.dt
    epoch = np.array([time[0], time[-1], '11111', T])
    # frac_rms = np.sqrt(np.var(lc_cho.counts) * lc_cho.counts.size / (lc_cho.counts.size - 1) / np.mean(lc_cho.counts) ** 2)
    # sigma_range=poisson_conf_interval(lc_cho.counts)
    # sigma=(sigma_range[1:,]-sigma_range[0,:])/2
    # mse=np.mean(sigma**2)
    # frac_rms=np.sqrt((np.var(lc_cho.counts) -mse) /np.mean(lc_cho.counts) ** 2)
    # print('frms=',frac_rms)
    # (psd_sim,result_mu) = bestpsd(lc_cho, epoch_info=epoch_info_use,maskfreq=maskfreq)
    # print(np.sum(lc_cho.counts))
    N_cts=int(np.sum(lc_cho.counts));bin=lc_cho.time[1]-lc_cho.time[0]
    t = np.random.random(np.random.poisson(N_cts)) * T + time[0]
    t = np.sort(t)
    evt = EventList()
    evt.time=t
    lc_out = evt.to_lc(dt=bin)

    return lc_out

def sim_LS_onesrc(path,dataname):
    max_LSP=[]
    for i in range(10):
        print(i)
        lc_out=singleobs2simevt(path,dataname)
        T = lc_out.time[-1] - lc_out.time[0]
        # freq = np.arange(1 / T, 1 / 1, 1e-5)
        [FP,out_period,max_NormLSP]=hawk.get_LS(lc_out.time, lc_out.counts, freq=np.arange(1 / T, 1 / 55, 1e-6),
                                                outpath=None, outname=None, save=0, show=1)
        # print(max_NormLSP)
        max_LSP.append(max_NormLSP)
    np.savetxt(path+dataname[:-3]+'.txt',max_LSP, fmt='%.20f', delimiter='\n')
    return max_LSP

def plot_LS_sig(path,dataname,outpath=None,outname=None,save=0,show=1):
     a = fits.open(path + dataname)
     rate = a[1].data['RATE']
     time = a[1].data['TIME']
     T = time[-1] - time[0]
     lc = Lightcurve(time=time, counts=rate * (time[1] - time[0]))
     print('dt=',time[1:]-time[0:-1])
     freq = np.arange(1 / T, 1 / 20, 1e-6)
     x = lc.time
     y = lc.counts
     LS = LombScargle(x, y, normalization='standard')
     power = LS.power(freq)
     max_NormLSP = np.max(power)
     period_peak = 1. / freq[np.where(power == np.max(power))][0]
     FP = LS.false_alarm_probability(power.max(), minimum_frequency=freq[0], maximum_frequency=freq[-1],
                                     method='baluev')
     FP_99 = LS.false_alarm_level(0.0027, minimum_frequency=freq[0], maximum_frequency=freq[-1], method='baluev')
     sim_LSP=np.loadtxt(path+dataname[:-3]+'.txt')
     sim_LSP_3sig=np.sort(sim_LSP)[-4]
     sim_LSP_2sig = np.sort(sim_LSP)[-50]
     Np = 2000
     FAP_N = 1 - (1 - FP) ** Np
     sigma=norm.ppf((1+1-FAP_N)/2)
     print(FP)
     print('FAP_N=',FAP_N)
     print('sigma=',sigma)
     print(sim_LSP_3sig)
     plt.figure(1, (10, 8))
     plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '-',label=r'$3 \sigma~(baluev)$')
     # plt.plot([freq[0], freq[-1]], [max_NormLSP, max_NormLSP], '--',label='{0:.2f} '.format(sigma)+r'$\sigma$')
     plt.plot([freq[0], freq[-1]], [sim_LSP_3sig, sim_LSP_3sig], '-',label=r'$3 \sigma~(simulated)$')
     plt.plot([freq[0], freq[-1]], [sim_LSP_2sig, sim_LSP_2sig], '-', label=r'$2 \sigma~(simulated)$')
     # plt.title('Period={0:.2f}'.format(period_peak), font1)
     plt.text(freq[np.where(power == np.max(power))][0] * 1.3, max_NormLSP * 0.95, 'P={0:.2f}s'.format(period_peak),
              fontsize=18, fontweight='semibold')
     plt.text(freq[np.where(power == np.max(power))][0] * 1.3, max_NormLSP * 0.55, 'sigma={0:.3f}'.format(sigma),
              fontsize=18, fontweight='semibold')

     plt.plot(freq, power)
     plt.semilogx()
     # plt.semilogy()
     out_period = 1. / freq[np.where(power == np.max(power))][0]
     # plt.text(freq[0], max_NormLSP, f'{sigma} '+r'$\sigma$', hawk.font1)
     # plt.text(freq[0], sim_LSP_3sig, r'3$\sigma$ (simulated)', hawk.font1)
     plt.legend()
     plt.xlabel('Frequency (Hz)', hawk.font1)
     plt.ylabel('Normalized LS Periodogram', hawk.font1)
     plt.tick_params(labelsize=16)
     if save:
         plt.savefig(outpath + outname + '_LS_trials.pdf', bbox_inches='tight', pad_inches=0.01)
     if show:
         plt.show()
     else:
         plt.close()
     return [FP, out_period, max_NormLSP]

if __name__=='__main__':
    path = '/Users/baotong/Desktop/XMMcentral/nustar_xmm/'
    dataname='nu80801332002A01_sr_3_10keV.lc'
    # file_names = os.listdir(path);list=[]
    # ###打印所有文件名
    # for file_name in file_names:
    #     if file_name.endswith(".lc"):
    #         list.append(file_name)
    #         plot_LS_sig(path,file_name,
    #                     outpath=path,outname=file_name[:-3],show=1,save=0)
    # plot_LS_sig(path,file_name,outpath=path,outname=file_name[:-3],show=1,save=0)
    
    # singleobs2simevt(path,dataname)
    sim_LS_onesrc(path, dataname)
    # def save_array_to_txt(array, filename):
    #     with open(filename, 'w') as f:
    #         np.savetxt(f, array, fmt='%.10f', delimiter='\n')
    # for dataname in list:
    #     max_LSPlist=sim_LS_onesrc(path, dataname)
    #     save_array_to_txt(array= max_LSPlist,filename=path + dataname[:-3] + '.txt')
    #     del max_LSPlist
#     list=['0840910501_268.37_24.77_18295_m1src_2_10keV.lc','0152920101_267.04_30.10_647_m1src_2_10keV.lc','0201200101_266.19_27.23_15815_m1src_0.2_10keV.lc',
# '0801683101_268.92_29.40_312_m1src_0.2_10keV.lc','0840910501_268.37_24.77_18295_m1src_2_10keV.lc',
# '0861171201_262.75_35.14_3868_m2src_2_10keV.lc',
# '0861171201_262.75_35.14_3868_m2src_2_10keV.lc']
    # dataname='0801683601_268.37_29.95_615_pnsrc_0.2_10keV.lc'
    # plot_LS_sig(path,dataname,outpath=path,outname=dataname[:-3],show=1,save=0)




