'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-10-21 20:24:04
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-11-28 11:45:37
FilePath: /code/match_eSASS/make_xmmdr14list.py
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
def filter_sources_by_exposure(catalog_fits, expmap_fits, output_fits):
    # Step 1: Load the catalog FITS file (the source list)
    with fits.open(catalog_fits) as hdul:
        catalog_data = hdul[1].data  # Assuming the data is in the second extension
        catalog_header = hdul[1].header  # To preserve the header for the new output FITS

        # Step 2: Load the exposure map and its WCS information
        with fits.open(expmap_fits) as exp_hdul:
            expmap_data = exp_hdul[0].data  # Exposure map data
            expmap_wcs = WCS(exp_hdul[0].header)  # WCS information for the exposure map

        # Step 3: Extract the RA and Dec from the catalog
        ra = catalog_data['RA']
        dec = catalog_data['DEC']

        # Step 4: Convert RA/Dec to pixel coordinates in the exposure map
        coords = SkyCoord(ra, dec, unit='deg', frame='fk5')
        x_pix, y_pix = expmap_wcs.world_to_pixel(coords)

        # Step 5: Filter sources that fall within the boundaries of the exposure map
        in_bounds = (
            (x_pix >= 0) & (x_pix < expmap_data.shape[1]) &
            (y_pix >= 0) & (y_pix < expmap_data.shape[0])
        )

        # Filter the sources that are within the exposure map boundary
        filtered_ra = ra[in_bounds]
        filtered_dec = dec[in_bounds]
        filtered_x_pix = x_pix[in_bounds].astype(int)
        filtered_y_pix = y_pix[in_bounds].astype(int)

        # Step 6: Check the exposure value at each valid source position
        exposure_at_sources = expmap_data[filtered_y_pix, filtered_x_pix]
        exposure_greater_than_zero = exposure_at_sources > 0
        # Step 7: Select sources with exposure > 0
        final_sources = catalog_data[in_bounds][exposure_greater_than_zero]

        # Step 8: Save the filtered sources into a new FITS file
        hdu = fits.BinTableHDU(data=final_sources, header=catalog_header)
        hdu.writeto(output_fits, overwrite=True)
        print(f"Filtered sources saved to {output_fits}")

def filter_and_save_indices(fits_file, selected_indices_file, unselected_indices_file):
    with fits.open(fits_file) as hdul:
        data = hdul[1].data
        header = hdul[1].header
        # Initialize lists to hold indices for selected and unselected real sources
        selected_indices = []
        unselected_indices = []
        # Temporary list to collect indices for the current source (real + associated)
        current_group_indices = []
        # Iterate over each row in the data
        for i, row in enumerate(data):
            print(i)  # For tracking progress
            if row['N_OBS'] > 0:  # Real source detected
                if current_group_indices:  # If there is already a group, process it
                    # Check if the last source was selected
                    real_source = data[current_group_indices[0]]
                    if (real_source['STACK_FLAG'] <= 1 and 
                        real_source['EP_DET_ML'] >= 6 and 
                        real_source['EXTENT'] == 0 and real_source['EXTENT_ML'] > 0):
                        selected_indices.extend(current_group_indices)
                    else:
                        unselected_indices.extend(current_group_indices)
                # Start a new group with the current real source
                current_group_indices = [i]
            else:
                # Add observational index to the current group
                current_group_indices.append(i)
        # Process the last group after the loop finishes
        if current_group_indices:
            real_source = data[current_group_indices[0]]
            if (real_source['STACK_FLAG'] <= 1 and 
                real_source['EP_DET_ML'] >= 6 and 
                real_source['EXTENT'] == 0 and real_source['EXTENT_ML'] > 0):
                selected_indices.extend(current_group_indices)
            else:
                unselected_indices.extend(current_group_indices)

    # Save the indices to text files
    np.savetxt(selected_indices_file, selected_indices, fmt='%d')
    np.savetxt(unselected_indices_file, unselected_indices, fmt='%d')

    print(f"Selected indices saved to {selected_indices_file}")
    print(f"Unselected indices saved to {unselected_indices_file}")

def save_fits_to_region(fits_file, region_output_file):
    # Step 1: Load the FITS file and extract the relevant columns
    with fits.open(fits_file) as hdul:
        catalog_data = hdul[1].data  # Assuming the data is in the second extension
        ra = catalog_data['RA']
        dec = catalog_data['DEC']
        stack_flag = catalog_data['STACK_FLAG']
        n_obs = catalog_data['N_OBS']
        ep_det_ml = catalog_data['EP_DET_ML']
        extent_ml = catalog_data['EXTENT_ML']
        extent = catalog_data['EXTENT']
    # Step 2: Write to the region file based on the specified conditions
    with open(region_output_file, 'w') as f:
        f.write('fk5\n')  # Write the header indicating the coordinate system
        for r, d, flag, nobs, ep_ml, ext_ml, ext in zip(ra, dec, stack_flag, n_obs, ep_det_ml, extent_ml, extent):
            if (nobs >= 1) and (flag <= 1) and (ep_ml > 6) and (ext == 0):
                f.write(f"circle({r},{d},{40/3600.0})\n")  # Circle with 50" radius

    print(f"Region file saved as {region_output_file}")
if __name__ == '__main__':
    path='/Users/baotong/data_GalDisc/data/'
    # filter_sources_by_exposure(catalog_fits=path+'xmmstack_v3.2_4xmmdr14s.fits.gz', 
    #                            expmap_fits=path+'mosaic_latest/'+'GalDisc_ima_2exp.fits.gz', 
    #                            output_fits=path+'GalDisc_4xmmdr14.fits')
    save_fits_to_region(fits_file=path+'xmmdr14s/GalDisc_4xmmdr14s_new_cleaned.fits',
                        region_output_file=path+'xmmdr14s/GalDisc_4xmmdr14s_new_cleaned.reg')
    # filter_and_save_indices(fits_file=path+'xmmdr14s/'+'GalDisc_4xmmdr14s.fits',
    #                     selected_indices_file=path+'xmmdr14s/'+'GalDisc_selected_newextent_indices.txt', 
    #                     unselected_indices_file=path+'xmmdr14s/'+'GalDisc_unselected_newextent_indices.txt')
    # xmmdr14_obs=fits.open(path+'xmmstack_v3.2_4xmmdr14s_obslist.fits.gz')
    # allobsnow=np.loadtxt(path+'allobs.txt',dtype='str')
    # obsused=np.intersect1d(allobsnow,xmmdr14_obs[1].data['OBS_ID'])
    # np.savetxt(path+'DR14_usedobsid.txt',obsused,fmt='%s')
    