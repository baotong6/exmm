'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-07-10 14:35:43
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-12-04 11:02:27
FilePath: /code/match_eSASS/oldcode/cal_FAP.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
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


def mindist_in_cat(coordinates,coordmatch=None):
    if coordmatch==None:
        coordmatch=coordinates
    min_distances = np.zeros(len(coordinates)) * u.deg
    for i, coord in enumerate(coordinates):
        dist = coord.separation(coordmatch)
        # 排除与自身的距离
        dist = dist[np.where(dist != 0 * u.deg)]
        # 找到最小距离
        min_dist = np.min(dist)
        # 存储最小距离
        min_distances[i] = min_dist
    return min_distances


# erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
# ra2=erass1list['RA']
# dec2=erass1list['DEC']
# pos2err=erass1list['POS_ERR']
# detlike0=erass1list['DET_LIKE_0']
# indexlist=np.arange(0,len(ra2)+1,1)
# indexlist=indexlist[np.where(detlike0>10)]

# ra2=ra2[indexlist];dec2=dec2[indexlist]
# pos2err=pos2err[indexlist]
# detlike0=detlike0[indexlist]

# expvalue=catalog_num.check_xmmfootprint(ra2,dec2)
# indexlist2=np.where(expvalue>100)[0]
# ra2=ra2[indexlist2]
# dec2=dec2[indexlist2]
# pos2err=pos2err[indexlist2]
# detlike0=detlike0[indexlist2]

# # matches=pd.read_csv(path+'match_xmm_eRASS1/'+'merge_res_10sec.txt', sep='\t',index_col=False)
# # ra1=matches.iloc[:,0];dec1=matches.iloc[:,1];ra2_em=matches.iloc[:,3];dec2_em=matches.iloc[:,4]
# # id2=matches.iloc[:,5];sep=matches.iloc[:,6];obsid=matches.iloc[:,7]
# # unique_elements = np.unique(id2)
# # num_unique = len(unique_elements)
# # id_erass=np.arange(1,len(ra2)+1,1)
# # id_nomatch=np.setdiff1d(id_erass,unique_elements)-1

# (ra_all,dec_all,err_all,obsid_all)=list_xmmfits.num_cat()
# (ra_xmmcut,dec_xmmcut,err_xmmcut,obsid_xmmcut)=list_xmmfits.cutxmm(ra_all,dec_all,err_all,obsid_all)
# plt.close()
path='/Users/baotong/data_GalDisc/data/'
erass_usedsrc=fits.open(path+'match_e_xmm/eRASS1_filtered.fits')[1].data
xmmfits=fits.open(path+'xmmdr14s/GalDisc_4xmmdr14s_new_cleaned.fits')[1].data
ncontrib=xmmfits['N_CONTRIB'];
realsrc=xmmfits[np.where((ncontrib>0))]
ra_xmmcut=realsrc['RA'];dec_xmmcut=realsrc['DEC']
# 判断点是否在凸多边形内部的函数
def point_in_hull(point, vertices):
    n = len(vertices)
    inside = False

    p1x, p1y = vertices[0]
    for i in range(1, n+1):
        p2x, p2y = vertices[i % n]
        if point[1] > min(p1y, p2y):
            if point[1] <= max(p1y, p2y):
                if point[0] <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (point[1] - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or point[0] <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside

def random():
    catalog1_coords=np.column_stack((erass_usedsrc['RA'],erass_usedsrc['DEC']))
    # 使用凸多边形算法找到所有坐标的边界
    hull = ConvexHull(catalog1_coords)
    # 获取凸多边形的顶点
    vertices = hull.points[hull.vertices]

    min_ra, max_ra = np.min(catalog1_coords[:, 0]), np.max(catalog1_coords[:, 0])
    min_dec, max_dec = np.min(catalog1_coords[:, 1]), np.max(catalog1_coords[:, 1])
    num_random_points = len(catalog1_coords)
    # 生成随机坐标，确保在凸多边形内
    for i in range(1,100):
        random_coordinates = []
        while len(random_coordinates) < num_random_points:
            ra = np.random.uniform(min_ra, max_ra)
            dec = np.random.uniform(min_dec, max_dec)
            # 检查坐标是否在凸多边形内
            if point_in_hull([ra, dec], vertices):
                random_coordinates.append([ra, dec])
        # 将生成的随机坐标转换为NumPy数组
        random_coordinates = np.array(random_coordinates)
        coords_match=SkyCoord(ra=random_coordinates[:,0]*u.deg,dec=random_coordinates[:,1]*u.deg)
        coords_xmm=SkyCoord(ra=ra_xmmcut*u.deg,dec=dec_xmmcut*u.deg)
        mindist_match=mindist_in_cat(coordinates=coords_match,coordmatch=coords_xmm).value*3600
        np.savetxt(path+'random_eRASSXMM/'+f'random_match{i}.txt',mindist_match)
    # plt.plot(vertices[:, 0], vertices[:, 1], 'r--', lw=2, label='Convex Hull')
    # plt.scatter(random_coordinates[:,0],random_coordinates[:,1])
    # plt.scatter(catalog1_coords[:,0],catalog1_coords[:,1])
    # plt.show()

def plot_random_res():
    mindist_all=[]
    for i in range(1,100):
        mindist=np.loadtxt(path+'random_eRASSXMM/'+f'random_match{i}.txt')
        mindist_all.extend(mindist)
    plt.hist(mindist_all,bins=np.logspace(np.log10(np.min(mindist_all)),
                                       np.log10(np.max(mindist_all)),50),
             histtype='step',cumulative=1,density=1,lw=2,color='red',label='eRASS-XMM min_sep')
    plt.semilogx()
    plt.semilogy()
    plt.plot([16,16],[0,1],'--',color='cyan')
    plt.xlabel("closest separation (arcsec)",fontdict=font1)
    plt.ylabel('FAP',fontdict=font1)
    plt.tick_params(labelsize=18)
    mindist_all=np.array(mindist_all)
    print(len(mindist_all),len(np.where(mindist_all<26)[0]))
    plt.show()
    
if __name__ == '__main__':
    # random()
    plot_random_res()


