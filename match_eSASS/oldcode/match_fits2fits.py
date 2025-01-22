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
import catalog_num

path='/Users/baotong/data_GalDisc/data/'
erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
ra2=erass1list['RA']
dec2=erass1list['DEC']
pos2err=erass1list['POS_ERR']
detlike0=erass1list['DET_LIKE_0']
ra2=ra2[np.where(detlike0>10)];dec2=dec2[np.where(detlike0>10)]
pos2err=pos2err[np.where(detlike0>10)]
detlike0=detlike0[np.where(detlike0>10)]

expvalue=catalog_num.check_xmmfootprint(ra2,dec2)
ra2=ra2[np.where(expvalue>100)]
dec2=dec2[np.where(expvalue>100)]
pos2err=pos2err[np.where(expvalue>100)]
detlike0=detlike0[np.where(expvalue>100)]
## Now the sample is the cutted eRASS1 sample

## using pos_err to cross-match ##
def match_2cat(ra_cat,dec_cat):
    file_names = os.listdir(path+'xmm_src_list/')
    for file_name in file_names:
        if file_name.endswith(".fits"):
            xmmsrclist = fits.open(path+f'xmm_src_list/{file_name}')
            ra1=xmmsrclist[1].data['RA']
            dec1=xmmsrclist[1].data['DEC']
            radec1_err=xmmsrclist[1].data['RADEC_ERR']
            sep_1sig=np.sqrt(np.mean(radec1_err)**2+np.mean(pos2err)**2)
            match_twofits(ra1,dec1,ra_cat,dec_cat,separation=10,
                        outpath=path+'match_xmm_eRASS1/res_obs/',outname=f'res_{file_name[:-23]}_10sec.txt')
            ##千万注意，这里的res中的id1和id2是index+1，不能直接用作索引, sorry for shit mountain.
def merge_res():
    merged_data = []
    for filename in os.listdir(path+'match_xmm_eRASS1/res_obs/'):
        if filename.endswith('_10sec.txt') and filename.startswith('res_'):
            file_path = os.path.join(path+'match_xmm_eRASS1/res_obs/', filename)
            with open(file_path, 'r') as file:
                lines = file.readlines()
                if lines:
                    file_info = os.path.splitext(filename)[0]
                    for line in lines:
                        data_row = line.strip().split()
                        data_row = [float(value) for value in data_row] 
                        data_row.append(filename[-20:-10])
                        merged_data.append(data_row)
    merged_array = np.array(merged_data)
    sorted_array = merged_array[merged_array[:,0].argsort()]

    df = pd.DataFrame(sorted_array,columns=['RA_xmm','DEC_xmm','ID_xmm',
                      'RA_e','DEC_e','ID_e','sep','ObsID'])
    # 将第三列和第六列转换为整数
    df.iloc[:, 2] = df.iloc[:, 2].astype(float).astype(int)
    df.iloc[:, 5] = df.iloc[:, 5].astype(float).astype(int)
    output_file = path+'match_xmm_eRASS1/'+'merge_res_10sec.txt'
    df.to_csv(output_file, sep='\t', index=False, header=True)
    # np.savetxt(output_file, sorted_array, fmt='%10s %10s %10s %10s %10s %10s %10s %10s')

def plot_match():
    matches=pd.read_csv(path+'match_xmm_eRASS1/'+'merge_res_10sec.txt', sep='\t')
    ra1=matches.iloc[:,0];dec1=matches.iloc[:,1];ra2=matches.iloc[:,3];dec2=matches.iloc[:,4]
    id2=matches.iloc[:,5];sep=matches.iloc[:,6];obsid=matches.iloc[:,7]
    plt.figure(1)
    plt.scatter(ra1,dec1,marker='o',s=30,edgecolor='green',facecolor='none')
    plt.scatter(ra2,dec2,marker='+',s=50,color='black')
    plt.savefig(path+'match_xmm_eRASS1/'+'matches_all_10sec.pdf')
    plt.close()
    plt.figure(2)
    plt.hist(sep,bins=np.logspace(np.log10(np.min(sep)),
                                  np.log10(np.max(sep)),20),histtype='step')
    plt.semilogx()
    plt.savefig(path+'match_xmm_eRASS1/'+'sep_histall_10sec.pdf')

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

def checkout_safeid():
    matches=pd.read_csv(path+'match_xmm_eRASS1/'+'merge_res_10sec.txt', sep='\t',index_col=False)
    ra1=matches.iloc[:,0];dec1=matches.iloc[:,1];ra2=matches.iloc[:,3];dec2=matches.iloc[:,4]
    id2=matches.iloc[:,5];sep=matches.iloc[:,6];obsid=matches.iloc[:,7]
    unique_elements = np.unique(id2)
    num_unique = len(unique_elements)

    safeeid=[];unsafeeid=[];midsafeeid=[];safemindist=[];unsafemindist=[];midsafemindist=[];
    ## safe means no multiple "physically" xmm counterparts for one eRASS source
    ## unsafe means multiple counterparts in single xmm obs
    ## midsafe means all xmm counterparts are one to one in xmm obs, but with large seperations
    ra1=np.array(ra1);dec1=np.array(dec1);obsid=np.array(obsid);sep=np.array(sep)
    for eid in unique_elements:
        index=np.where(id2==eid)[0]
        if len(index)>1:
            coords=SkyCoord(ra=ra1[index]*u.deg,dec=dec1[index]*u.deg)
            # print(coords)
            coords_xmmobsid=obsid[index]
            mindist=mindist_in_cat(coords).value*3600
            if np.max(mindist)<5 and len(coords_xmmobsid) == len(set(coords_xmmobsid)):
                safeeid.append(eid)
                safemindist.append(np.max(mindist))
            elif len(coords_xmmobsid) != len(set(coords_xmmobsid)):
                unsafeeid.append(eid)
                unsafemindist.append(np.max(mindist))
            else:
                midsafeeid.append(eid)
                midsafemindist.append(np.max(mindist))

        else:safeeid.append(eid)
    print(len(safeeid),len(unsafeeid),len(midsafemindist),num_unique)
    plt.hist(np.concatenate((safemindist,unsafemindist,midsafemindist)),bins=30,histtype='step')
    plt.savefig(path+'mindist_hist_amongxmmsrc.pdf')
    plt.close()
    print(unsafeeid)
    all_indices=[]
    for value in safeeid:
        indices = [index for index, a_value in enumerate(id2) if a_value == value]
        all_indices.extend(indices)
    safeeid_sep=sep[all_indices]

    plt.hist(safeeid_sep,bins=np.logspace(np.log10(np.min(safeeid_sep)),
                                          np.log10(np.max(safeeid_sep)),40),histtype='step')
    plt.semilogx()                                
    plt.savefig(path+'safeeid_matchsep.pdf')
    print(unsafeeid,midsafeeid)
    return (safeeid,unsafeeid,midsafeeid)


            
if __name__ == '__main__':
    # match_2cat(ra_cat=ra2,dec_cat=dec2)
    # merge_res()
    # plot_match()
    (safeeid,unsafeeid,midsafeeid)=checkout_safeid()

    