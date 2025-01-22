'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-07-10 13:28:49
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-07-10 14:16:37
FilePath: /code/catalog_num.py
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
from match_list import match_twofits


path='/Users/baotong/data_GalDisc/data/'
erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
ra2=erass1list['RA']
dec2=erass1list['DEC']
pos2err=erass1list['POS_ERR']
detlike0=erass1list['DET_LIKE_0']

def plot_erass1():
    path='/Users/baotong/eSASS/data/eRASS/'
    erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
    ra2=erass1list['RA']
    dec2=erass1list['DEC']
    pos2err=erass1list['POS_ERR']
    plt.hist(pos2err,bins=np.linspace(0,10,100),histtype='step')
    plt.show()
    coord=SkyCoord(ra=ra2*u.deg,dec=dec2*u.deg)
    esass_galcoord=coord.transform_to('galactic')

    index1=np.where((esass_galcoord.l.value>320)|(esass_galcoord.l.value<9))
    index2=np.where((esass_galcoord.b.value>-3)&(esass_galcoord.b.value<2))
    all_index=np.intersect1d(index1,index2)
    print(len(index1[0]),len(index2[0]),len(all_index))

    ra_exmm=ra2[all_index]
    dec_exmm=dec2[all_index]
    pos2err_exmm=pos2err[all_index]
    plt.scatter(ra_exmm,dec_exmm,marker='+',color='black')
    plt.show()
def check_xmmfootprint(ra,dec):
    expmap=fits.open(path+'mosaic_latest/'+'GalDisc_ima_4exp.fits.gz')
    exp_data=expmap[0].data
    exp_data=exp_data.T
    w = WCS(path+'mosaic_latest/'+'GalDisc_ima_4exp.fits.gz')
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    src_x = np.rint(src_x)
    src_y = np.rint(src_y)
    src_x = src_x.astype(np.int64)
    src_y = src_y.astype(np.int64)
    src_x-=1;src_y-=1
    A=exp_data;B=np.column_stack((src_x,src_y))
    C = np.zeros(len(B))
    mask = (B[:, 0] >= 0) & (B[:, 0] < A.shape[0]) & (B[:, 1] >= 0) & (B[:, 1] < A.shape[1])
    C[mask] = A[B[mask][:, 0], B[mask][:, 1]]
    expvalue_out=C
    return expvalue_out

if __name__ == '__main__':
    print(len(ra2))
    ra2=ra2[np.where(detlike0>10)]
    dec2=dec2[np.where(detlike0>10)]
    print(len(ra2))
    expvalue=check_xmmfootprint(ra2,dec2)
    print(len(np.where(expvalue>100)[0]))








