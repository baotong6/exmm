import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from scipy import interpolate
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord,match_coordinates_sky
import pandas as pd

def read_and_format_file_to_dataframe(file_path):
    """
    Read the file line by line, format each element, and store the result in a Pandas DataFrame.
    Args:
    file_path (str): The file path.
    Returns:
    pandas.DataFrame: DataFrame containing the formatted content of each line.
    """
    data_list = []
    with open(file_path, 'r') as file:
        for line in file:
            data = line.strip().split(',')

            formatted_data = []
            for item in data:
                try:
                    formatted_data.append(float(item))
                except ValueError:
                    formatted_data.append(item)
            data_list.append(formatted_data)
    df = pd.DataFrame(data_list)
    return df

def load_data():
    path='/home/baotong/eRASS1/data/'
    formatted_df = read_and_format_file_to_dataframe(path+'xmm_src_list/all_source_list.txt')
    ra1=np.array(formatted_df[1][1:])
    dec1=np.array(formatted_df[2][1:])
    
    erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
    ra2=erass1list['RA']
    dec2=erass1list['DEC']

    return (ra1,dec1,ra2,dec2)

def match_twofits(ra1,dec1,ra2,dec2,separation=1,outpath='',outname='match_res.txt'):
    c1 = SkyCoord(ra=ra1 * u.degree, dec=dec1 * u.degree)
    c2 = SkyCoord(ra=ra2 * u.degree, dec=dec2 * u.degree)
    idx, d2d, d3d = match_coordinates_sky(matchcoord=c1, catalogcoord=c2, nthneighbor=1)
    max_sep =separation* u.arcsec
    sep_constraint = d2d < max_sep
    c_matches = c1[sep_constraint]
    catalog_matches = c2[idx[sep_constraint]]
    d2d_matches = d2d[sep_constraint]
    judge=np.array(sep_constraint).astype('int')
    match_result=np.column_stack((np.arange(1,len(ra1)+1,1),idx+1,judge))
    # np.savetxt(outpath+f'match_sep{separation}.txt',match_result,fmt='%10d %10d %10d')
    match_index=np.where(judge==1)[0]
    outra1=ra1[match_index];outdec1=dec1[match_index]
    outra2=ra2[idx[match_index]];outdec2=dec2[idx[match_index]]
    sep=d2d[match_index].value*3600
    outresult=np.column_stack((outra1,outdec1,match_index+1,outra2,outdec2,idx[match_index]+1,sep))
    np.savetxt(outpath+outname,outresult,fmt='%10.5f %10.5f %10d %10.5f %10.5f %10d %10.5f')
    # for i in range(len(match_index)):
    #     plt.scatter(ra1[match_index[i]],dec1[match_index[i]],marker='o',s=30,edgecolor='green',facecolor='none')
    #     plt.scatter(ra2[idx[match_index[i]]],dec2[idx[match_index[i]]],marker='+',s=50,color='black')
    # plt.savefig(outpath+'sep.pdf')
    # plt.close()
    # plt.figure(1)
    # plt.hist(sep,bins=50,histtype='step')
    # plt.savefig(outpath+'sep_hist.pdf')
    return idx[sep_constraint]

if __name__ == '__main__':
    (ra1,dec1,ra2,dec2)=load_data()
    match_twofits(ra1,dec1,ra2,dec2,separation=10,outpath='/home/baotong/eRASS1/data/')

