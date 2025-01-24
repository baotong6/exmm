'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-11-25 17:03:59
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2025-01-24 12:24:48
FilePath: /code/match_eSASS/filter_xmmdr14list.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
from astropy.io import fits
import numpy as np
def filter_fits_by_indices(input_fits, output_fits, index_list):
    """
    根据索引列表过滤FITS文件的table，保存结果到新文件。
    
    :param input_fits: str, 输入的FITS文件路径
    :param output_fits: str, 输出的FITS文件路径
    :param index_list: list, 包含要保留的行号的索引列表
    """
    # 打开输入FITS文件
    with fits.open(input_fits) as hdul:
        # 获取数据和头文件
        original_header = hdul[1].header
        original_data = hdul[1].data
        
        # 过滤数据，保留索引列表中的行
        filtered_data = original_data[index_list]
        
        # 创建一个新的HDUList，保留原始的header
        new_hdu = fits.BinTableHDU(data=filtered_data, header=original_header)
        
        # 创建完整的HDUList
        new_hdul = fits.HDUList([hdul[0], new_hdu])
        
        # 保存到新的FITS文件
        new_hdul.writeto(output_fits, overwrite=True)
        
        print(f"新的FITS文件已保存到: {output_fits}")

# 示例调用
path='/Users/baotong/data_GalDisc/data/xmmdr14s/'
input_fits = path+"GalDisc_4xmmdr14s.fits"  
output_fits = path+ "GalDisc_4xmmdr14s_newextent_cleaned.fits" 
index_list = np.loadtxt(path+'GalDisc_selected_newextent_indices.txt') 
index_list = np.array(index_list, dtype=int)
filter_fits_by_indices(input_fits, output_fits, index_list)
