'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-10-21 10:24:10
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2025-01-21 13:44:57
FilePath: /code/match_eSASS/make_eRASS1list.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
from astropy.io import fits
from astropy.table import Table
import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

def check_xmmfootprint(ra, dec):
    path = '/Users/baotong/data_GalDisc/data/'
    expmap = fits.open(path + 'mosaic_latest/' + 'GalDisc_ima_2exp.fits.gz')
    exp_data = expmap[0].data.T
    w = WCS(expmap[0].header)
    src_x, src_y = w.all_world2pix(ra, dec, 1)
    src_x = np.rint(src_x).astype(np.int64) - 1
    src_y = np.rint(src_y).astype(np.int64) - 1
    exp_values = np.zeros(len(ra))
    mask = (src_x >= 0) & (src_x < exp_data.shape[0]) & (src_y >= 0) & (src_y < exp_data.shape[1])
    exp_values[mask] = exp_data[src_x[mask], src_y[mask]]
    return exp_values

def filter_sources_by_exposure(catalog_fits, exposure_fits, threshold=0):
    with fits.open(catalog_fits) as cat_hdu:
        catalog_data = cat_hdu[1].data
        ra = catalog_data['ra']
        dec = catalog_data['dec']
    with fits.open(exposure_fits) as exp_hdu:
        exposure_data = exp_hdu[0].data
        exposure_wcs = WCS(exp_hdu[0].header)
    source_coords = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    x_pix, y_pix = exposure_wcs.world_to_pixel(source_coords)
    x_pix = np.round(x_pix).astype(int)
    y_pix = np.round(y_pix).astype(int)
    in_bounds_mask = (x_pix >= 0) & (x_pix < exposure_data.shape[1]) & \
                     (y_pix >= 0) & (y_pix < exposure_data.shape[0])
    valid_sources_mask = np.zeros(len(x_pix), dtype=bool)
    valid_sources_mask[in_bounds_mask] = exposure_data[y_pix[in_bounds_mask], x_pix[in_bounds_mask]] > threshold
    valid_indices = np.where(valid_sources_mask)[0]
    return valid_indices

def filter_paras(catalog_fits):
    srclist = fits.open(catalog_fits)[1].data
    indexlist = np.where((srclist['DET_LIKE_0'] > 10) & (srclist['EXT_LIKE'] == 0))[0]
    return indexlist

def save_filtered_sources(input_fits, output_fits, indices):
    with fits.open(input_fits) as hdul:
        original_header = hdul[1].header
        original_data = hdul[1].data
        filtered_data = original_data[indices]  # 根据筛选的索引提取数据
        hdu = fits.BinTableHDU(data=filtered_data, header=original_header)
        hdu.writeto(output_fits, overwrite=True)
    print(f"已保存筛选后的数据到 {output_fits}")

if __name__ == '__main__':
    path = '/Users/baotong/data_GalDisc/data/'
    catalog_fits = path + 'match_e_xmm/' + 'eRASS1_Main.v1.1.fits'
    exposure_fits = path + 'mosaic_latest/' + 'GalDisc_ima_2exp.fits.gz'
    output_fits = path + 'eRASS1_filtered_bydet10.fits'

    xmmfootprint_eid = filter_sources_by_exposure(catalog_fits=catalog_fits, exposure_fits=exposure_fits, threshold=1)
    filterindex = filter_paras(catalog_fits=catalog_fits)
    outindex = np.intersect1d(filterindex, xmmfootprint_eid)
    print(f"筛选后源数量: {len(outindex)}")
    save_filtered_sources(input_fits=catalog_fits, output_fits=output_fits, indices=outindex)
