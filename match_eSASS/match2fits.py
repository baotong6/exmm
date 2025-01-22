'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-10-21 11:21:19
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2025-01-21 12:41:38
FilePath: /code/match_eSASS/match2fits.py
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
def load_2list(erassfits,xmmcsv):    
    hdul=fits.open(erassfits)[1].data
    xmmfootprint_eid=filter_sources_by_exposure(catalog_fits=path+'match_e_xmm/'+'eRASS1_Main.v1.1.fits',
                               exposure_fits=path+'mosaic_latest/'+'GalDisc_ima_2exp.fits.gz',
                               threshold=1)
    filterindex=filter_paras(catalog_fits=path+'match_e_xmm/'+'eRASS1_Main.v1.1.fits')
    outindex=np.intersect1d(filterindex,xmmfootprint_eid)
    c1_data=hdul[outindex]
    c1_ra = c1_data['ra']
    c1_dec = c1_data['dec']
    # Load c2 from the CSV file
    c2_df = pd.read_csv(xmmcsv)
    c2_ra = c2_df['ra'].values
    c2_dec = c2_df['dec'].values

    coords_e = SkyCoord(ra=c1_ra * u.degree, dec=c1_dec * u.degree)
    coords_xmm = SkyCoord(ra=c2_ra * u.degree, dec=c2_dec * u.degree)

    return coords_e,coords_xmm


def load_2fits():    
    # hdul=fits.open(path+'match_e_xmm/'+'eRASS1_Main.v1.1.fits')[1].data
    # xmmfootprint_eid=filter_sources_by_exposure(catalog_fits=path+'match_e_xmm/'+'eRASS1_Main.v1.1.fits',
    #                            exposure_fits=path+'mosaic_latest/'+'GalDisc_ima_2exp.fits.gz',
    #                            threshold=1)
    # filterindex=filter_paras(catalog_fits=path+'match_e_xmm/'+'eRASS1_Main.v1.1.fits')
    # outindex=np.intersect1d(filterindex,xmmfootprint_eid)
    # c1_data=hdul[outindex]
    c1_data=fits.open(path+'match_e_xmm/'+'eRASS1_filtered.fits')[1].data
    c1_ra = c1_data['RA']
    c1_dec = c1_data['DEC']
    c2_data=fits.open('/Users/baotong/data_GalDisc/data/xmmdr14s/GalDisc_4xmmdr14s_new_cleaned.fits')[1].data
    ncontrib=c2_data['N_CONTRIB'];
    # flag=c2_data['STACK_FLAG'];detml=c2_data['EP_DET_ML'];extent=c2_data['EXTENT']
    # realsrc=c2_data[np.where((ncontrib>0)&(flag<=1)&(detml>6)&(extent<=0))]
    realsrc=c2_data[np.where((ncontrib>0))]
    print('number of realsrc',len(realsrc['RA']))
    c2_ra = realsrc['RA']
    c2_dec = realsrc['DEC']
    coords_e = SkyCoord(ra=c1_ra * u.degree, dec=c1_dec * u.degree)
    coords_xmm = SkyCoord(ra=c2_ra * u.degree, dec=c2_dec * u.degree)

    return coords_e,coords_xmm

def plot_nearest_sep(coord1,coord2):
    nearest_distances = []
    # For each source in c2, find the nearest source in c1 and compute the distance
    for i, coord_c1 in enumerate(coord1):
        # Compute the separations to all sources in c1
        separations = coord_c1.separation(coord2)
        # Find the minimum separation (nearest neighbor)
        min_sep = separations.min()
        # Store the nearest distance
        nearest_distances.append(min_sep.arcsecond)  # Distance in arcseconds
    # Convert to numpy array
    nearest_distances = np.array(nearest_distances)
    # Plot the distance distribution
    bins=np.logspace(-1,3,30)
    plt.hist(nearest_distances, bins=bins, edgecolor='black',histtype='step',lw=2)
    plt.plot([16,16],[0,50],'--',color='cyan',lw=2)
    plt.semilogx()
    plt.semilogy()
    plt.xlabel('Nearest neighbor distance (arcseconds)',fontsize=16)
    plt.ylabel('Number of sources',fontsize=16)
    plt.tick_params(labelsize=16)
    plt.title('Distribution of nearest neighbor distances between c2 and c1',fontsize=16)
    plt.show()
def match_exmm(esasscoord, xmmcoord, output_file):
    """
    Find the nearest source in xmm_coords for each source in esass_coords and save the results to a CSV file.
    Parameters:
        esasscoord (SkyCoord): SkyCoord object for the ESASS sources.
        xmmcoord (SkyCoord): SkyCoord object for the XMM sources.
        output_file (str): Path to save the output CSV file.
    """
    # Initialize output data
    output_data = []
    nomatch_count = 0
    # Find the nearest XMM source for each ESASS source
    for idx, esass_coord in enumerate(esasscoord):
        separation = esass_coord.separation(xmmcoord)
        nearest_index = separation.argmin()
        nearest_sep = separation[nearest_index].arcsecond
        # Handle sources with a separation > 10 arcseconds
        if nearest_sep > 15:
            nomatch_count += 1
        # Append result for the current ESASS source
        output_data.append({
            'e_id': idx,  # Assign a unique ID to each ESASS source (1-based index)
            'e_ra': esass_coord.ra.deg,
            'e_dec': esass_coord.dec.deg,
            'xmm_index': nearest_index,  # Assign a unique ID to each XMM source (1-based index)
            'xmm_ra': xmmcoord[nearest_index].ra.deg,
            'xmm_dec': xmmcoord[nearest_index].dec.deg,
            'separation_arcsec': nearest_sep
        })

    # Convert to DataFrame and save as CSV
    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_file, index=False)

    print(f"Number of unmatched sources with separation > 15 arcseconds: {nomatch_count}")

def plot_POSerr():
    c1_data=fits.open(path+'match_e_xmm/'+'eRASS1_filtered.fits')[1].data
    c2_data=fits.open('/Users/baotong/data_GalDisc/data/xmmdr14s/GalDisc_4xmmdr14s_new_cleaned.fits')[1].data
    ncontrib=c2_data['N_CONTRIB'];
    realsrc=c2_data[np.where((ncontrib>0))]
    poserr1=c1_data['POS_ERR']
    poserr2=realsrc['RADEC_ERR']
    plt.hist(poserr1,bins=20,histtype='step',lw=2,color='k',label=r'1$\sigma$ positional error of eRASS1 sources')
    plt.hist(poserr2,bins=30,histtype='step',lw=2,color='green',label=r'1$\sigma$ positional error of XMM sources')
    plt.xlabel('Arcsec',fontsize=16)
    plt.ylabel('Number of Sources',fontsize=16)
    plt.tick_params(labelsize=16)
    plt.semilogy()
    plt.legend()
    # plt.title(r'1$\sigma$ positional error of eRASS1 and XMM sources',fontsize=16)
    plt.show()

if __name__ == '__main__':
    path='/Users/baotong/data_GalDisc/data/'
    erassfits=path+'match_e_xmm/'+'eRASS1_Main.v1.1.fits'
    # xmmcsv=path+'mosaic_latest/clean_srclist/'+'filtered_ppsxmmsrclis_corr_ML14_EXT0_goodobs_exp_allreg_ft2sig.csv'
    (coords_e,coords_xmm)=load_2fits()
    # match_exmm(esasscoord=coords_e,xmmcoord=coords_xmm,
    #            output_file=path+'match_e_xmm/e_xmmdr14s_match_all.csv')
    # plot_nearest_sep(coords_e,coords_xmm)
    plot_POSerr()

