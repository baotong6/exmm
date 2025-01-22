
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
import subprocess
import csv
import pandas as pd
import pandas as pd
from astropy.stats import poisson_conf_interval
import math
import stingray as sr
from stingray.lightcurve import Lightcurve
from astropy.io import fits
import rednoise as rednoise
import hawkeye as hawk
import corner
import numpy as np
import emcee
import matplotlib.pyplot as plt
from rednoise import useful_functions as func

font1=func.font1
font2=func.font2
# 定义 Power Law 模型
def powerlaw(x, alpha, beta):
    return alpha * x**(-1)+beta

# 定义似然函数（likelihood）
def ln_likelihood(theta, x, y):
    alpha, beta = theta
    y_model = powerlaw(x, alpha, beta)
    residuals = y - y_model
    chi_squared = np.sum(residuals**2)
    return -0.5 * chi_squared

# 定义先验分布（这里假设均匀先验）

def ln_prior(theta):
    alpha, beta = theta
    if 0.0 < alpha < 10.0 and -5.0 < beta < 5.0:
        return 0.0
    return -np.inf

# 定义全局似然函数（包括先验）
def ln_posterior(theta, x, y):
    ln_prior_val = ln_prior(theta)
    if not np.isfinite(ln_prior_val):
        return -np.inf
    return ln_prior_val + ln_likelihood(theta, x, y)

def plot_mcmc(x,y,CR,filename=None,path=None):
    # 初始化 MCMC 链
    ndim = 2  # 参数的数量
    nwalkers = 100  # Walker 的数量
    nsteps = 200  # 步数
    # nburn = 100  # "burn-in" period to let chains stabilize
    # 创建 EnsembleSampler 对象
    sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_posterior, args=(x, y))
    # 运行 MCMC 链
    # 随机初始化 Walker 的起始位置
    initial_guess = [0,2/CR]  # 初始参数猜测
    starting_guesses = np.zeros((nwalkers, ndim))
    starting_guesses[:, 0] = initial_guess[0]+ 1*np.random.random(nwalkers)
    starting_guesses[:, 1] = initial_guess[1] + np.random.random(nwalkers)

    sampler.run_mcmc(starting_guesses, nsteps)
    # 提取样本
    samples = sampler.chain[:, :, :].reshape((-1, ndim))
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    # trace = sampler.chain[:, nburn:, :].reshape(-1, ndim).T
    trace=samples
    # corner.corner(flat_samples,truths=[0.01,2/CR])
    a = corner.corner(flat_samples, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})
    mu1 = corner.quantile(flat_samples[:, 0], q=[0.5])[0]
    mu2 = corner.quantile(flat_samples[:, 1], q=[0.5])[0]
    plt.savefig(path+f'cornerplot_{file_name}.pdf')
    plt.close()

    p0=trace[:,0];p1=trace[:,1]
    xfit = x
    print(len(xfit))
    print(len(p0))
    print(len(p1))
    # yfit = p0[:,None]*xdata**(-1)*(1+(xdata/1e-2)**4)**(-1/4.)+p1[:,None]
    yfit = p0[:, None] * xfit[::10] ** (-1) +p1[:,None]
    print(yfit)
    mu = yfit.mean(0)
    sig = 1 * yfit.std(0)
    print(mu,sig)
    plt.figure(figsize=(8, 6))
    plt.step(x, y,label='PSD')
    # plt.plot(xfit, mu, '-k',label=r'$\rm Model: P(\nu)=N \nu^{-1}(1+(\frac{\nu}{\nu_0})^4)^{-1/4}+C_p$')
    plt.plot(xfit[::10], mu, '-k', label=r'$\rm Model: P(\nu)=N \nu^{-1}+C_p$')
    plt.fill_between(xfit[::10], mu - 1*sig, mu + 1*sig, color = 'lightgray')
    plt.plot(xfit[::10],np.zeros(len(xfit[::10]))+2/CR,'--',color='orange',label='Poisson noise')
    # ax.plot(np.zeros(len(yfit))+2012.0,yfit,'--',color='yellow')
    plt.xlabel('Frequency (Hz)', font2)
    plt.ylabel(r'$\rm Power~([rms/mean]^2 Hz^{-1})$', font2)
    plt.legend()
    plt.tick_params(labelsize=16)
    plt.loglog()
    plt.savefig(path+f'fitplot_{file_name}.pdf')
    # plt.figure(figsize=(8, 6))
    # plt.step(x, y, label="Data")
    # for alpha, beta in samples[np.random.randint(len(samples), size=100)]:
    #     plt.plot(x, powerlaw(x, alpha, beta), color="k", alpha=0.1)
    # plt.loglog()
    # plt.xlabel('Frequency (Hz)', font2)
    # plt.ylabel(r'$\rm Power~([rms/mean]^2 Hz^{-1})$', font2)
    # plt.legend()
    # plt.savefig(path+f'fitplot_{file_name}.pdf')


if __name__=='__main__':
    path='/Users/baotong/Desktop/XMMcentral/nustar_xmm/'
    # 列出当前目录下的所有文件和文件夹
    file_names = os.listdir(path)
    ##打印所有文件名
    file_names=['nu80801332002A01_sr_3_10keV.lc']
    for file_name in file_names:
        if file_name.endswith(".lc"):
            print(file_name)
            a=fits.open(path+file_name)
            rate=a[1].data['RATE']
            time=a[1].data['TIME']
            T=time[-1]-time[0]
            print(T)
            lc=Lightcurve(time=time,counts=rate*(time[1]-time[0]))
            freq=np.arange(1/T,1/20,1e-6)
            # hawk.get_LS(lc.time,lc.counts,freq,show=1)
            CR = np.mean(lc.counts) / lc.dt
            epoch=np.array([time[0],time[-1],'11111',T])
            hawk.get_LS(time,rate,freq=np.arange(1/T,1/20,1e-6),outpath=None,outname=None,save=0,show=1)
            rednoise.plot_psd(lc, show=0)
            psd = rednoise.plot_psd(lc, norm='frac', show=0, ifexpTfilter=T)
            x= np.array(psd.freq);
            y= np.array(psd.power)
            plot_mcmc(x,y,CR,file_name,path)
    # file_name='0886110501_270.42_23.71_901_m2src_0.2_10keV.lc'
    # a = fits.open(path + file_name)
    # rate = a[1].data['RATE']
    # time = a[1].data['TIME']
    # T = time[-1] - time[0]
    # lc = Lightcurve(time=time, counts=rate * (time[1] - time[0]))
    # CR = np.mean(lc.counts) / lc.dt
    # epoch=np.array([time[0],time[-1],'11111',T])
    # # hawk.get_LS(time,rate,freq=np.arange(1/T,1/1,1e-6),outpath=None,outname=None,save=0,show=1)
    # # rednoise.plot_psd(lc, show=1)
    # psd = rednoise.plot_psd(lc, norm='frac', show=0, ifexpTfilter=T)
    # x= np.array(psd.freq);
    # y= np.array(psd.power)
    # plot_mcmc(x,y,CR,file_name,path)

    