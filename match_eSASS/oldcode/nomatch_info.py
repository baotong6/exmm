'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-07-16 12:20:19
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-09-27 09:29:20
FilePath: /code/nomatch_info.py
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
import csv
# from match_list import match_twofits
# from match_fits2fits import checkout_safeid,mindist_in_cat
import catalog_num
import list_xmmfits

path='/Users/baotong/data_GalDisc/data/'
erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
mindex=np.loadtxt('/Users/baotong/data_GalDisc/data/'+'e_matches_index.txt')
nomindex=np.loadtxt('/Users/baotong/data_GalDisc/data/''e_nomatches_index.txt')
mindex=mindex.astype('int64')
nomindex=nomindex.astype('int64')
matches=erass1list[mindex];no_matches=erass1list[nomindex]
matches=matches[np.where(matches['EXT_LIKE']==0)]
no_matches=no_matches[np.where(no_matches['EXT_LIKE']==0)]

def esrc_obsid():
    (ra_all,dec_all,err_all,obsid_all)=list_xmmfits.num_cat()
    (ra_xmmcut,dec_xmmcut,err_xmmcut,obsid_xmmcut)=list_xmmfits.cutxmm(ra_all,dec_all,err_all,obsid_all)
    plt.close()

    c1 = SkyCoord(ra=no_matches['RA'] * u.degree, dec=no_matches['DEC'] * u.degree)
    c2 = SkyCoord(ra=ra_xmmcut * u.degree, dec=dec_xmmcut * u.degree)
    idx, d2d, d3d = match_coordinates_sky(matchcoord=c1, catalogcoord=c2, nthneighbor=1)
    esrc_obsid=obsid_xmmcut[idx]

    return esrc_obsid

def plot_expvalue():
    expvalue_match=catalog_num.check_xmmfootprint(matches['RA'],matches['DEC'])
    expvalue_nomatch=catalog_num.check_xmmfootprint(no_matches['RA'],no_matches['DEC'])

    print(expvalue_match,expvalue_nomatch)
    plt.figure(1,(10,8))
    plt.hist(expvalue_match,bins=np.logspace(np.log10(np.min(expvalue_match)),
                                    np.log10(np.max(expvalue_match)),20),histtype='step',color='red',lw=2,label='eRASS matches')
    plt.hist(expvalue_nomatch,bins=np.logspace(np.log10(np.min(expvalue_nomatch)),
                                    np.log10(np.max(expvalue_nomatch)),7),histtype='step',color='green',lw=2,label='eRASS nomatches')
    plt.loglog()
    plt.legend()
    plt.xlabel('Expvalue (second)',fontsize=18)
    plt.ylabel('Number of Sources',fontsize=18)
    plt.tick_params(labelsize=18)
    plt.savefig(path+'hist_expvalue.pdf')
    # plt.close()
    plt.show()
    
def plot_flux_nomatch():
    x=np.arange(1,len(no_matches['ML_FLUX_1'])+1,1)
    y=no_matches['ML_FLUX_1']
    yerr=no_matches['ML_FLUX_ERR_1']
    yerr1=no_matches['ML_FLUX_LOWERR_1']
    yerr2=no_matches['ML_FLUX_UPERR_1']
    netcounts=no_matches['ML_CTS_1']
    print('number of nomatches',len(no_matches))
    expvalue_nomatch=catalog_num.check_xmmfootprint(no_matches['RA'],no_matches['DEC'])
    goodno_matches_index=np.where(expvalue_nomatch>1500)
    badno_matches_index=np.where(expvalue_nomatch<=1500)
    # print(yerr1,yerr2)
    y2=np.zeros(len(y))+1e-14
    y3=np.zeros(len(y))+5e-15
    plt.errorbar(x[goodno_matches_index],y[goodno_matches_index],
                 yerr=yerr[goodno_matches_index],fmt='.', 
                 capsize=5,elinewidth=2.5, ecolor='black',color='black',linewidth=2.0)
    plt.errorbar(x[badno_matches_index],y[badno_matches_index],
                 yerr=yerr[badno_matches_index],fmt='.', 
                 capsize=5,elinewidth=1.5, ecolor='gray',color='gray',linewidth=1.0)
    plt.plot(x,y2,'--',label=r'$5 \sigma$ sensitivity in 0.5-2 keV (time=1000s)')
    plt.plot(x,y3,'--',label=r'$5 \sigma$ sensitivity in 0.5-2 keV (time=5000s)')
    plt.xlabel('ID of the source',fontsize=18)
    plt.ylabel(r'$\rm Flux~in~erg~cm^{-2}~s^{-1}$',fontsize=18)
    plt.tick_params(labelsize=18)
    plt.semilogy()
    plt.legend()
    plt.show()
    # print(netcounts)
    # print(no_matches['RA'],no_matches['DEC']

if __name__ == '__main__':
    # esrc_obsid()
    # plot_expvalue()
    plot_flux_nomatch()

