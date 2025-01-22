'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-09-24 18:33:45
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-10-01 10:46:39
FilePath: /code/matches_info.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
import os,sys
import pandas as pd
from convert_date_MJD import convert_to_mjd
path='/Users/baotong/data_GalDisc/data/'
erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
mindex=np.loadtxt('/Users/baotong/data_GalDisc/data/'+'e_matches_index.txt')
nomindex=np.loadtxt('/Users/baotong/data_GalDisc/data/''e_nomatches_index.txt')
mindex=mindex.astype('int64')
nomindex=nomindex.astype('int64')
matches=erass1list[mindex];no_matches=erass1list[nomindex]
matches=matches[np.where(matches['EXT_LIKE']==0)]
# print(len(matches))
ra_e_match=matches['RA'];dec_e_match=matches['DEC']
factor_cps2flux=9.800E-13 

def calculate_ratio_with_error(a, b, delta_a, delta_b):
    """
    计算a/b及其误差
    :param a: 数值a
    :param b: 数值b
    :param delta_a: a的误差
    :param delta_b: b的误差
    :return: 比值和误差
    """
    # 计算比值
    ratio = a / b
    
    # 计算相对误差
    relative_error = np.sqrt((delta_a / a)**2 + (delta_b / b)**2)
    
    # 计算比值的误差
    delta_ratio = ratio * relative_error
    
    return (ratio, delta_ratio)

def output_xmmmatchesinfo():
    # RA_xmm	DEC_xmm	ID_xmm	RA_e	DEC_e	ID_e  sep	 ObsID #
    matchesall=pd.read_csv(path+'match_xmm_eRASS1/'+'merge_res_10sec.txt', sep='\t',index_col=False)
    ra1=matchesall.iloc[:,0];dec1=matchesall.iloc[:,1];ra2=matchesall.iloc[:,3];dec2=matchesall.iloc[:,4]
    id2=matchesall.iloc[:,5];sep=matchesall.iloc[:,6];obsid=matchesall.iloc[:,7];idxmm=matchesall.iloc[:,2]
    # unique_elements = np.unique(id2)
    # num_unique = len(unique_elements)

    safeeid=[];unsafeeid=[];midsafeeid=[];safemindist=[];unsafemindist=[];midsafemindist=[];
    ## safe means no multiple "physically" xmm counterparts for one eRASS source
    ## unsafe means multiple counterparts in single xmm obs
    ## midsafe means all xmm counterparts are one to one in xmm obs, but with large seperations
    ra1=np.array(ra1);dec1=np.array(dec1);obsid=np.array(obsid);sep=np.array(sep)
    ra2=np.array(ra2);dec2=np.array(dec2);idxmm=np.array(idxmm)

        # 读取文件名并存储到列表中
    file_names = []
    with open(path+'xmm_src_list/all.txt', 'r') as file:
        file_names = [line.strip() for line in file]
    xmmmatchinfolist=[]
    for i in range(len(ra_e_match)):
        a1=obsid[np.where((np.abs(ra2-ra_e_match[i])<1e-5)&(np.abs(dec2-dec_e_match[i])<1e-5))]
        a2=idxmm[np.where((np.abs(ra2-ra_e_match[i])<1e-5)&(np.abs(dec2-dec_e_match[i])<1e-5))]
        xmmsrcinfo_list=[]
        for j in range(len(a1)):
            for file_name in file_names:
                if str(a1[j]) in file_name:
                    xmmsrclist=fits.open(path+'xmm_src_list/'+file_name)[1].data
                    xmmsrcinfo=xmmsrclist[a2[j]-1]
                    xmmra=xmmsrcinfo['RA'];xmmdec=xmmsrcinfo['DEC']
                    coordxmm=SkyCoord(ra=xmmra*u.degree,dec=xmmdec*u.degree)
                    coorde=SkyCoord(ra=ra_e_match[i]*u.degree,dec=dec_e_match[i]*u.degree)
                    sep=coordxmm.separation(coorde).arcsec
                    rate=xmmsrcinfo['RATE'];rate_err=xmmsrcinfo['RATE_ERR']
            xmmsrcinfo_list.append([xmmra,xmmdec,rate,rate_err,rate*factor_cps2flux,rate_err*factor_cps2flux,sep,
                                    xmmsrcinfo['SCTS'],xmmsrcinfo['SCTS_pn_1'],xmmsrcinfo['SCTS_pn_2'],xmmsrcinfo['SCTS_pn_3'],
                                    xmmsrcinfo['SCTS_pn_4'],xmmsrcinfo['SCTS_pn_5'],xmmsrcinfo['SCTS_m1_1'],xmmsrcinfo['SCTS_m1_2'],
                                    xmmsrcinfo['SCTS_m1_3'],xmmsrcinfo['SCTS_m1_4'],xmmsrcinfo['SCTS_m1_5'],xmmsrcinfo['SCTS_m2_1'],
                                    xmmsrcinfo['SCTS_m2_2'],xmmsrcinfo['SCTS_m2_3'],xmmsrcinfo['SCTS_m2_4'],xmmsrcinfo['SCTS_m2_5']])
        xmmmatchinfolist.append([a1,a2,xmmsrcinfo_list])

    return (matches,xmmmatchinfolist)

def readobsinfo(search_obsid):
    search_obsid=str(search_obsid).zfill(10)
    #这里的search_obsid是十位的字符串
    filename='obsinfo.csv'
    obsinfo=pd.read_csv(path+'xmm_src_list/'+filename,dtype=str)

    # 查找对应的obsid并获取时间列的值
    time_value = obsinfo.loc[obsinfo['OBSERVATION.OBSERVATION_ID'] == search_obsid, 'OBSERVATION.START_UTC'].values[0]
    MJD=convert_to_mjd(time_value)
    return MJD
def plot_fluxcomp(matches,xmmmatchinfolist):
    era=matches['RA'];edec=matches['DEC']
    ename=matches['IAUNAME']
    eflux=matches['ML_FLUX_1']
    efluxerr=matches['ML_FLUX_ERR_1']
    eMJD=matches['MJD']
    for i in range(len(era)):
        xmmmatchinfo=xmmmatchinfolist[i]
        xmmobsid=xmmmatchinfo[0]
        plt.figure(1)
        for j in range(len(xmmobsid)):
            xmmMJD=readobsinfo(xmmobsid[j])
            xmmra=xmmmatchinfo[2][j][0]
            xmmdec=xmmmatchinfo[2][j][1]
            xmmflux=xmmmatchinfo[2][j][4]
            xmmfluxerr=xmmmatchinfo[2][j][5]
            plt.errorbar(xmmMJD,xmmflux,
                 yerr=xmmfluxerr,fmt='.', 
                 capsize=5,elinewidth=2.5, ecolor='black',
                 color='black',linewidth=2.0,label='XMM')
        plt.errorbar(eMJD[i],eflux[i],
                    yerr=efluxerr[i],fmt='.', 
                    capsize=5,elinewidth=1.5, ecolor='red',
                    color='gray',linewidth=1.0,label='eROSITA')
        plt.xlabel('MJD',fontsize=18)
        plt.ylabel(r'$\rm 0.2-2.3~keV flux~in~erg~cm^{-2}~s^{-1}$',fontsize=18)
        plt.tick_params(labelsize=18)
        # 自定义图例，只显示黑色和红色的图例
        handles, labels = plt.gca().get_legend_handles_labels()
        # 仅保留黑色和红色的图例项
        by_label = dict(zip(labels, handles))
        plt.legend([by_label['XMM'], by_label['eROSITA']], 
                ['XMM', 'eROSITA'])
        plt.savefig(path+f'match_eRASS_figs/{ename[i]}_flux.pdf',
                    bbox_inches='tight', pad_inches=0.05)
        plt.close()
    # print(xmmra)
def plot_fluxratio(matches,xmmmatchinfolist):
    era=matches['RA'];edec=matches['DEC']
    ename=matches['IAUNAME']
    eflux=matches['ML_FLUX_1']
    efluxerr=matches['ML_FLUX_ERR_1']
    eMJD=matches['MJD']
    ratiolist=[];ratioerrlist=[];labellist=[]
    xmmfluxlist1=[];xmmfluxerrlist1=[]
    xmmfluxlist2=[];xmmfluxerrlist2=[]
    for i in range(len(eflux)):
        xmmmatchinfo=xmmmatchinfolist[i]
        xmmobsid=xmmmatchinfo[0]
        xmmflux=np.array(xmmmatchinfo[2])[:,4]
        xmmfluxerr=np.array(xmmmatchinfo[2])[:,5]
        xmmfluxlist1.append(np.max(xmmflux));
        xmmfluxerrlist1.append(xmmfluxerr[np.argmax(xmmflux)])
        xmmfluxlist2.append(np.min(xmmflux))
        xmmfluxerrlist2.append(xmmfluxerr[np.argmin(xmmflux)])
    
        if eflux[i]>=np.max(xmmflux):
            (ratio,ratioerr)=calculate_ratio_with_error(eflux[i], np.min(xmmflux),
                                   efluxerr[i], xmmfluxerr[np.argmin(xmmflux)])
            label=1
        elif eflux[i]<=np.min(xmmflux):
            (ratio,ratioerr)=calculate_ratio_with_error(np.max(xmmflux), eflux[i],
                                    xmmfluxerr[np.argmax(xmmflux)],efluxerr[i])
            label=3
        elif np.min(xmmflux)<eflux[i]<np.max(xmmflux):
            (ratio1,ratio1err)=calculate_ratio_with_error(np.max(xmmflux),eflux[i],
                                   xmmfluxerr[np.argmax(xmmflux)],efluxerr[i])
            (ratio2,ratio2err)=calculate_ratio_with_error(eflux[i], np.min(xmmflux),
                                   efluxerr[i], xmmfluxerr[np.argmin(xmmflux)])
            ratio=np.where(ratio1 >= ratio2, ratio1, ratio2)
            ratioerr=np.where(ratio1 >= ratio2, ratio1err, ratio2err)
            label=2
        ratiolist.append(ratio);ratioerrlist.append(ratioerr)
        labellist.append(label)
        
# 计算 c_err 数组
    labellist=np.array(labellist)
    maxS = np.select([labellist == 1, labellist == 2, labellist == 3], 
                     [eflux, xmmfluxlist1, xmmfluxlist1])
    maxSerr=np.select([labellist == 1, labellist == 2, labellist == 3], 
                     [efluxerr, xmmfluxerrlist1, xmmfluxerrlist1])
    colors = np.where(labellist == 1, 'black', np.where(labellist == 2, 'green', 'red'))
    ratiolist=np.array(ratiolist);ratioerrlist=np.array(ratioerrlist)
    # print(labellist)
    # print(maxS.shape)       # 查看 maxS 的形状
    # print(ratiolist.shape)  # 查看 ratiolist 的形状
    # print(maxSerr.shape)    # 查看 maxSerr 的形状
    # print(ratioerrlist.shape) # 查看 ratioerrlist 的形状
    # print(colors.shape)
    # print(colors)
    plt.figure(1,figsize=(10,8))
    for i in range(len(maxS)):
        plt.errorbar(maxS[i], ratiolist[i], xerr=maxSerr[i], yerr=ratioerrlist[i], fmt='o', 
                    capsize=3, elinewidth=0.5, ecolor=colors[i],
                    color=colors[i], linewidth=1.0)
    plt.xlabel(r'$\rm 0.2-2.3~keV~flux~in~erg~cm^{-2}~s^{-1}$',fontsize=18)
    plt.ylabel(r'$\rm Flux~ratio~(S_{max}/S_{min})$',fontsize=18)
    plt.tick_params(labelsize=18)
    plt.loglog()
    black_patch = plt.Line2D([0], [0], marker='.', color='w', 
                           label='eROSITA (emax)', markerfacecolor='black', markersize=10)
    green_patch = plt.Line2D([0], [0], marker='.', color='w', 
                             label='XMM (emid)', markerfacecolor='green', markersize=10)
    red_patch = plt.Line2D([0], [0], marker='.', color='w', 
                              label='XMM (emin)', markerfacecolor='red', markersize=10)

    plt.legend(handles=[black_patch, green_patch, red_patch], title='')
    plt.savefig(path+f'match_eRASS_figs/ratio_flux.pdf',
                bbox_inches='tight', pad_inches=0.05)
    # plt.show()
    transones=np.argsort(ratiolist)[-10:]
    coord1=SkyCoord(ra=matches[transones]['RA']*u.deg,dec=matches[transones]['DEC']*u.deg)
    # print(coord1)
    # for k in transones:
    #     print(xmmmatchinfolist[k][0],xmmmatchinfolist[k][1])
    #     print(np.array(xmmmatchinfolist[k][2])[:,-1])
    print(maxS)
    return maxS
if __name__ == '__main__':
    (matches,xmmmatchinfolist)=output_xmmmatchesinfo()
    # plot_fluxcomp(matches,xmmmatchinfolist)
    plot_fluxratio(matches,xmmmatchinfolist)
    # readobsinfo('0932201001')
