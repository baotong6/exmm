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
import catalog_num 

path='/Users/baotong/data_GalDisc/data/'
erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
ra2=erass1list['RA']
dec2=erass1list['DEC']
pos2err=erass1list['POS_ERR']
detlike0=erass1list['DET_LIKE_0']
ra2=ra2[np.where(detlike0>10)];dec2=dec2[np.where(detlike0>10)]
pos2err=pos2err[np.where(detlike0>10)]

def num_cat():
    ra_all=[]
    dec_all=[]
    err_all=[]
    obsid_all=[]
    file_names = os.listdir(path+'xmm_src_list/')
    for file_name in file_names:
        if file_name.endswith(".fits"):
            xmmsrclist = fits.open(path+f'xmm_src_list/{file_name}')
            ra1=xmmsrclist[1].data['RA']
            dec1=xmmsrclist[1].data['DEC']
            radec1_err=xmmsrclist[1].data['RADEC_ERR']
            obsid=[file_name[-33:-23] for i in range(len(ra1))]
            # sep_1sig=np.sqrt(np.mean(radec1_err)**2+np.mean(pos2err)**2)
            ra_all.extend(list(ra1));dec_all.extend(list(dec1));err_all.extend(list(radec1_err))
            obsid_all.extend(obsid)
    ra_all=np.array(ra_all);dec_all=np.array(dec_all)
    err_all=np.array(err_all);obsid_all=np.array(obsid_all)
    
    return (ra_all,dec_all,err_all,obsid_all)

def dist_all(ra_all,dec_all,err_all,coordmatch):
    coordinates = SkyCoord(ra=ra_all*u.deg, dec=dec_all*u.deg)
    # 初始化一个数组来存储每个点源到最近邻点源的距离
    min_distances = np.zeros(len(coordinates)) * u.deg
    for i, coord in enumerate(coordinates):
        # print(i)
        dist = coord.separation(coordmatch)
        # 排除与自身的距离
        dist = dist[np.where(dist != 0 * u.deg)]
        # 找到最小距离
        min_dist = np.min(dist)
        # 存储最小距离
        min_distances[i] = min_dist
    return min_distances

def cutesass_rough(ra_all,dec_all,err_all):
## useless , goodbye ##
    ## input from xmm coord
    xmmcoord=SkyCoord(ra=ra_all*u.deg,dec=dec_all*u.deg).transform_to('galactic')
    l=xmmcoord.l.value;b=xmmcoord.b.value

    print(l[np.where(l>200)].min(),l[np.where(l<200)].max())
    print(np.max(xmmcoord.b.value),np.min(xmmcoord.b.value))

    erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
    ra2=erass1list['RA']
    dec2=erass1list['DEC']
    pos2err=erass1list['POS_ERR']
    coord=SkyCoord(ra=ra2*u.deg,dec=dec2*u.deg)
    esass_galcoord=coord.transform_to('galactic')

    index1=np.where((esass_galcoord.l.value>349)|(esass_galcoord.l.value<7.5))
    index2=np.where((esass_galcoord.b.value>-2.5)&(esass_galcoord.b.value<2))
    all_index=np.intersect1d(index1,index2)
    print(len(index1[0]),len(index2[0]),len(all_index))

    ra_exmm=ra2[all_index]
    dec_exmm=dec2[all_index]
    pos2err_exmm=pos2err[all_index]
    return (ra_exmm,dec_exmm,pos2err_exmm)

def cutesass_byexp(ra_all,dec_all):

    expvalue=catalog_num.check_xmmfootprint(ra2,dec2)
    ra2_c=ra2[np.where(expvalue>100)]
    dec2_c=dec2[np.where(expvalue>100)]
    pos2err_c=pos2err[np.where(expvalue>100)]

    return (ra2_c,dec2_c,pos2err_c)

def cutxmm(ra_all,dec_all,err_all,obsid_all):
    xmmcoord=SkyCoord(ra=ra_all*u.deg,dec=dec_all*u.deg).transform_to('galactic')
    l=xmmcoord.l.value;b=xmmcoord.b.value

    print(l[np.where(l>200)].min(),l[np.where(l<200)].max())
    print(np.max(xmmcoord.b.value),np.min(xmmcoord.b.value))
    cutindex=np.where(l>179)[0]
    ra_out=ra_all[cutindex];dec_out=dec_all[cutindex];err_out=err_all[cutindex]
    obsid_out=obsid_all[cutindex]
    return (ra_out,dec_out,err_out,obsid_out)



def plot_mindist():
    min_dist=np.loadtxt(path+'min_dist_xmm2erosita.txt')
    plt.hist(min_dist,bins=np.logspace(np.log10(np.min(min_dist)),
                                       np.log10(np.max(min_dist)),50),histtype='step')
    plt.loglog()
    plt.xlim(0.1,1000)
    plt.savefig(path+'min_dist_xmm2erosita.pdf')


if __name__ == '__main__':
    (ra_all,dec_all,err_all,obsid_all)=num_cat()
    (ra_xmmcut,dec_xmmcut,err_xmmcut,obsid_out)=cutxmm(ra_all,dec_all,err_all,obsid_all)
    (ra_exmm,dec_exmm,pos2err_exmm)=cutesass_byexp(ra_all,dec_all)
    print(len(ra_exmm),len(ra_xmmcut))
    plt.scatter(ra_all,dec_all,marker='o',s=10,edgecolor='gray',
                facecolor='none',label='XMM sources')
    plt.scatter(ra_xmmcut,dec_xmmcut,marker='o',s=10,edgecolor='green',
                facecolor='none',label='XMM sources in western Galactic hemisphere')
    plt.scatter(ra_exmm,dec_exmm,marker='+',s=15,color='black',
                label='eROSITA sources in XMM-obs footprints')
    plt.legend()
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.savefig(path+'xmmerosita_src.pdf')
    coord_erosita=SkyCoord(ra=ra_exmm*u.deg,dec=dec_exmm*u.deg)
    min_distances=dist_all(ra_xmmcut,dec_xmmcut,err_xmmcut,coord_erosita)
    min_dist=min_distances.value*3600
    np.savetxt(path+'min_dist_xmm2erosita.txt',min_dist)
    plot_mindist()
    
