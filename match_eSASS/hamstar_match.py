'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-12-02 13:38:13
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2025-01-22 14:35:47
FilePath: /code/match_eSASS/hamstar_match.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.coordinates import match_coordinates_sky
import pandas as pd
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from make_eRASS1list import filter_sources_by_exposure,filter_paras
def plot_sthflux():
    path = '/Users/baotong/data_GalDisc/data/'
    # Read the match data
    ematchinfo = pd.read_excel(path + 'match_e_xmm/e_xmmdr14s_match_all.xlsx')
    # Open the FITS file
    xmmdr14slist = fits.open(path + 'xmmdr14s/GalDisc_4xmmdr14s_new_cleaned.fits')[1].data
    # Filter sources with N_CONTRIB > 0
    ncontrib = xmmdr14slist['N_CONTRIB']
    realsrc = xmmdr14slist[np.where(ncontrib > 0)]
    print(f"Number of valid sources: {len(realsrc)}")
    
    # Get matching sources with separation < 17
    sep = ematchinfo['separation_arcsec']
    match_src = ematchinfo[sep < 17]
    # Check if match_src is empty
    if match_src.empty:
        print("No sources matched with separation < 17")
        return
    # Get xmm_index values for matched sources
    xmm_indexes = match_src['xmm_index'].values
    hamindex= match_src['hamstar'].values
    # Access the corresponding sources from the FITS file
    exmmsrc = realsrc[xmm_indexes]
    hamsrc=realsrc[xmm_indexes[np.where(hamindex>0)]]

    # print(exmmsrc['N_CONTRIB'],exmmsrc['EP_HR1'])
    plt.errorbar(x=exmmsrc['EP_CTS'], xerr=exmmsrc['EP_CTS_ERR'], 
                 y=exmmsrc['EP_FLUX'], yerr=exmmsrc['EP_FLUX_ERR'], 
                 fmt='.', capsize=2, elinewidth=1, ecolor='gray', 
                 color='gray', linewidth=1)
    plt.loglog()
    plt.show()
    nan_count = np.sum(np.isnan(exmmsrc['EP_CTS']))
    print(f"Number of NaN values in EP_CTS: {nan_count}")    
    bins1=np.logspace(-15,-10,20)
    bins2=np.logspace(-16,-10,100)
    print(len(hamsrc),len(exmmsrc))
    plt.hist(exmmsrc['EP_FLUX'],bins=bins1,histtype='step',lw=2,color='k',label='XMM-eRASS1 sources')
    plt.hist(hamsrc['EP_FLUX'], bins=bins1, histtype='stepfilled', lw=2, label='XMM-eRASS1-hamstar sources', 
             alpha=0.5, color='orange', edgecolor='k',hatch='//')
    plt.hist(realsrc['EP_FLUX'],bins=bins2,histtype='step',lw=1,color='g',label='XMM sources')
    plt.legend()
    plt.xlabel(r'$\rm Flux~density~(erg~cm^{-2}~s^{-1})$',fontsize=18)
    plt.ylabel('Number of Sources',fontsize=18)
    plt.tick_params(labelsize=18)
    plt.loglog()
    plt.show()


def lum_hamstar_matches():
    path = '/Users/baotong/data_GalDisc/data/'
    hamstar = fits.open(path + 'match_e_xmm/hamstar.fit')[1].data
    ematchinfo = pd.read_excel(path + 'match_e_xmm/e_xmmdr14s_match_all.xlsx',sheet_name='all')
    hamindex=ematchinfo['hamstar']
    sep=ematchinfo['separation_arcsec']
    valid_hamindex=hamindex[np.where((hamindex>0)&(sep<17))[0]]
    valid_hamstar=hamstar[valid_hamindex]
    para=valid_hamstar['plx']
    distkpc=1/para

    xmmdr14slist = fits.open(path + 'xmmdr14s/GalDisc_4xmmdr14s_new_cleaned.fits')[1].data
    ncontrib = xmmdr14slist['N_CONTRIB']
    realsrc = xmmdr14slist[np.where(ncontrib > 0)]
    xmm_indexes = ematchinfo['xmm_index'].values
    hamstar_src=realsrc[xmm_indexes[np.where((hamindex>0)&(sep<17))[0]]]
    flux_xmm=hamstar_src['EP_FLUX']
    Lx_xmm=(flux_xmm*u.erg/(u.cm**2*u.s)*4*3.14*(distkpc*u.kpc)**2).to(u.erg/u.s).value
    bins1=np.logspace(1,3.3,15)
    plt.hist(distkpc*1000,bins=bins1,lw=3,color='gray',histtype='step')
    plt.semilogx()
    plt.xlabel('Distance (pc)',fontsize=18)
    plt.ylabel('Number of Sources',fontsize=18)
    plt.tick_params(labelsize=18)
    plt.show()
    bins=np.logspace(27,33,20)
    plt.hist(Lx_xmm,bins=bins,histtype='step',lw=4,color='orange')
    # print(Lx_xmm)
    plt.xlabel('X-ray Luminosity (erg/s)',fontsize=18)
    plt.ylabel('Number of Sources',fontsize=18)
    plt.tick_params(labelsize=18)
    plt.semilogx()
    plt.show()
    return None

def tonycut_fig():
    path = '/Users/baotong/data_GalDisc/data/'
    xmmdr14slist = fits.open(path + 'xmmdr14s/GalDisc_4xmmdr14s_new_cleaned.fits')[1].data
    ncontrib = xmmdr14slist['N_CONTRIB']
    realsrc = xmmdr14slist[np.where(ncontrib > 0)]
    hamstar = fits.open(path + 'match_e_xmm/HamStar_eRASS1_Main_Likely_Identifications_v1.1.fits')[1].data
    ematchinfo = pd.read_excel(path + 'match_e_xmm/e_xmmdr14s_match_all_starinfo.xlsx')

    xmm_indexes = ematchinfo['xmm_index'].values
    exmmsrc = realsrc[xmm_indexes]
    Fx_xmm=exmmsrc['EP_FLUX']
    Fx_xmm_err=exmmsrc['EP_FLUX_ERR']
    
    starsrc=ematchinfo[ematchinfo['hamstarindex']>0]
    
    para=np.array(starsrc['PLX']);Gmag=np.array(starsrc['G']);
    BP_RP=np.array(starsrc['BP_RP']);eFx=np.array(starsrc['Fx']);
    xmmFx=Fx_xmm[ematchinfo['hamstarindex']>0]
    sep=np.array(starsrc['separation_arcsec'])
    Fx_use=np.zeros(len(starsrc))
    for i in range(len(sep)):
        if sep[i]>17:
            Fx_use[i]=eFx[i]
        else:Fx_use[i]=xmmFx[i]
    xmm_index=np.array(starsrc['xmm_index'])
    distkpc=1/para;distpc=distkpc*1000
    msun=-26.7 
    absmsun=4.66
    Lsun=3.846e33*(u.erg/u.s)  ## all in G band

    abs_Gmag = Gmag - 5 * (np.log10(distpc) - 1)
    Fopt=10**(0.4*(msun-Gmag))*Lsun/(4*np.pi*(1*u.au)**2)
    Fopt=Fopt.to(u.erg/u.cm**2/u.s).value
    distkpc=1/para
    print(distpc)
    plt.scatter(BP_RP,abs_Gmag,marker='.',color='k')
    plt.xlabel('BP-RP',fontsize=16)
    plt.ylabel('absGmag',fontsize=16)
    plt.gca().invert_yaxis()
    plt.tick_params(labelsize=16)
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.scatter(BP_RP,Fx_use/Fopt,marker='.',color='black')
    x1 = np.linspace(-0.4, 0.7, 50)
    y1 = 10 ** (x1 - 3.5)
    x2 = np.linspace(0.7,5,50)
    y2 = 10 ** (x2- 3)
    plt.plot(x1, y1,'-', color='blue')
    plt.plot(x2, y2, label='modified cut in (Rodriguez+,2024)', color='blue')
    plt.plot([0.7,0.7],[y1[-1],y2[0]],'--',color='blue')
    plt.plot([1.5,1.5],[1e-6,1e2],'--',color='gray')
    plt.plot([-0.3,-0.3],[1e-6,1e2],'--',color='gray')
    plt.semilogy()
    
    print(Fx_use)
    print(Fopt)
    x=BP_RP
    logy=np.log10(Fx_use/Fopt)
    index_compact=np.where((logy>x-3)&(x<1.5))[0]
    plt.scatter(x[index_compact],10**logy[index_compact],marker='s',color='red',s=30)
    plt.xlabel('Gaia BP-RP',fontsize=16)
    plt.ylabel(r'$\rm F_x/F_{\rm opt}$',fontsize=16)
    plt.tick_params(labelsize=16)
    plt.legend()
    plt.semilogy()
    plt.show()


if __name__=='__main__':
    # plot_sthflux()
    # lum_hamstar_matches()
    tonycut_fig()