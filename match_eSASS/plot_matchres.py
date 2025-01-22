'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-10-21 16:15:50
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-12-02 10:09:44
FilePath: /code/match_eSASS/plot_matchres.py
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
from make_eRASS1list import filter_sources_by_exposure,filter_paras
def load_3list():
    path='/Users/baotong/data_GalDisc/data/'
    eres=pd.read_csv(path+'match_e_xmm/'+'e_xmmdr14s_match_all.csv')

    edata=fits.open(path+'match_e_xmm/'+'eRASS1_Main.v1.1.fits')[1].data
    xmmfootprint_eid=filter_sources_by_exposure(catalog_fits=path+'match_e_xmm/'+'eRASS1_Main.v1.1.fits',
                               exposure_fits=path+'mosaic_latest/'+'GalDisc_ima_2exp.fits.gz',
                               threshold=1)
    filterindex=filter_paras(catalog_fits=path+'match_e_xmm/'+'eRASS1_Main.v1.1.fits')
    outindex=np.intersect1d(filterindex,xmmfootprint_eid)
    edata_filt=edata[outindex]
    xmmdata = pd.read_csv(path+'mosaic_latest/clean_srclist/'+
                          'filtered_ppsxmmsrclis_corr_ML14_EXT0_goodobs_exp_allreg_ft2sig.csv')

    return eres,edata,xmmdata
def save_df_to_region_file(df, region_output_file):
    """
    Saves the 'ra' and 'dec' values from a DataFrame to a DS9 region file in fk5 format.
    Sources with 'xmm_index' == 0 will be saved as cyan circles, others as green circles.
    Parameters:
    df (pd.DataFrame): DataFrame containing 'e_ra', 'e_dec', and 'xmm_index' columns.
    region_output_file (str): The output file path for the region file.
    """
    # Check if required columns exist
    if 'e_ra' not in df.columns or 'e_dec' not in df.columns or 'xmm_index' not in df.columns:
        raise ValueError("The DataFrame must contain 'e_ra', 'e_dec', and 'xmm_index' columns.")

    with open(region_output_file, 'w') as f:
        f.write('fk5\n')  # Header for fk5 coordinate system
        for _, row in df.iterrows():
            ra = row['e_ra']
            dec = row['e_dec']
            xmm_index = row['xmm_index']
            # Choose color based on xmm_index value
            color = 'cyan' if xmm_index == 0 else 'green'
            # Write the circle with the chosen color
            f.write(f"circle({ra},{dec},{50/3600.0}) # color={color}\n")
    
    print(f"Sources saved to region file {region_output_file}")

    return None

def make_nomatch_reg():
    path='/Users/baotong/data_GalDisc/data/match_e_xmm/'
    df=pd.read_csv(path+'e_xmmdr14s_match_all.csv')
    nomatch_list=np.where(df['separation_arcsec']>15)[0]
    region_output_file=path+'e_xmmdr14s_nomatch.reg'
    with open(region_output_file, 'w') as f:
        f.write('fk5\n')  # Header for fk5 coordinate system
        for _, row in df.iterrows():
            ra = row['e_ra']
            dec = row['e_dec']
            xmm_index = row['xmm_index']
            eid=row['e_id']
            sep=row['separation_arcsec']
            # Choose color based on xmm_index value
            if sep>16.5:
                color='cyan'
            # color = 'cyan' if xmm_index == 0 else 'magenta'
            # Write the circle with the chosen color
                f.write(f"circle({ra},{dec},{60/3600.0}) # color={color}  width=3 font=\"helvetica 24 bold roman\" text={{{int(eid)}}} \n")


def plot_match_res():
    path='/Users/baotong/data_GalDisc/data/match_e_xmm/'
    df=pd.read_csv(path+'e_xmmdr14s_match_all.csv')
    sep=df['separation_arcsec']
    bins=np.logspace(-1,3,25)
    plt.hist(sep,bins=bins,histtype='step',lw=2)
    plt.semilogx()
    plt.show()
if __name__ == '__main__':
    path='/Users/baotong/data_GalDisc/data/match_e_xmm/'
    # (eres,edata,xmmdata)=load_3list()
    # save_df_to_region_file(df=eres, region_output_file=path+'region_ematch.reg')
    plot_match_res()
    # make_nomatch_reg()

    