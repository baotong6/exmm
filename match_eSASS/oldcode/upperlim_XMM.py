'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-09-21 15:01:22
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-09-22 13:45:29
FilePath: /code/upperlim_XMM.py
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
import math
from scipy.special import gammainc, gamma, gammaincinv

'''
fx is the source’s ﬂux in a given energy band, 
ti,Bi,EEFi are the exposure time, background level, 
and encircled energy fraction of the camera i, 
ECFi is the energy to count conversion factor of the detector i.
'''

def calculate_UL(N, CL, B,t,ECF,EEF):
    # 计算不完全伽马函数的逆
    X=CL * (gamma(N + 1)-gammainc(N + 1, B)*gamma(N + 1)) 
    + gammainc(N + 1, B)*gamma(N + 1)
    print(X)
    # 计算 UL
    UL = gammaincinv(N+1,X/gamma(N + 1)) - B
    CR_UL=UL/(t*EEF)
    fx=UL/(t*ECF*EEF)

    return UL

# 示例参数
N = 10  # 例如，N = 10
CL = 0.68  # 例如，CL = 0.95
B = 5   # 例如，B = 1

UL = calculate_UL(N, CL, B)
print("UL:", UL)

# print(gamma([11],[4]))
