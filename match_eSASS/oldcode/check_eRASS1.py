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
from match_fits2fits import checkout_safeid,mindist_in_cat
import catalog_num
import list_xmmfits
path='/Users/baotong/data_GalDisc/data/'
erass1list=fits.open(path+'eRASS1_Main.v1.1.fits')[1].data
ra2=erass1list['RA']
dec2=erass1list['DEC']
pos2err=erass1list['POS_ERR']
detlike0=erass1list['DET_LIKE_0']
indexlist=np.arange(0,len(ra2)+1,1)
indexlist=indexlist[np.where(detlike0>10)]

ra2=ra2[indexlist];dec2=dec2[indexlist]
pos2err=pos2err[indexlist]
detlike0=detlike0[indexlist]

expvalue=catalog_num.check_xmmfootprint(ra2,dec2)
indexlist2=np.where(expvalue>100)[0]
ra2=ra2[indexlist2]
dec2=dec2[indexlist2]
pos2err=pos2err[indexlist2]
detlike0=detlike0[indexlist2]


matches=pd.read_csv(path+'match_xmm_eRASS1/'+'merge_res_10sec.txt', sep='\t',index_col=False)
ra1=matches.iloc[:,0];dec1=matches.iloc[:,1];ra2_em=matches.iloc[:,3];dec2_em=matches.iloc[:,4]
id2=matches.iloc[:,5];sep=matches.iloc[:,6];obsid=matches.iloc[:,7]
unique_elements = np.unique(id2)
num_unique = len(unique_elements)
(ra_all,dec_all,err_all,obsid_all)=list_xmmfits.num_cat()
(ra_xmmcut,dec_xmmcut,err_xmmcut,obsid_xmmcut)=list_xmmfits.cutxmm(ra_all,dec_all,err_all,obsid_all)
plt.close()

id_erass=np.arange(1,len(ra2)+1,1)
id_nomatch=np.setdiff1d(id_erass,unique_elements)-1
# (safeeid,unsafeeid,midsafeeid)=checkout_safeid()

matches=erass1list[indexlist][indexlist2][unique_elements-1]
no_matches=erass1list[indexlist][indexlist2][id_nomatch]
np.savetxt('/Users/baotong/data_GalDisc/data/'+'e_matches_index.txt',
           indexlist[indexlist2][unique_elements-1],fmt='%d')
np.savetxt('/Users/baotong/data_GalDisc/data/'+'e_nomatches_index.txt',
           indexlist[indexlist2][id_nomatch],fmt='%d')

matches=matches[np.where(matches['EXT_LIKE']==0)]
no_matches=no_matches[np.where(no_matches['EXT_LIKE']==0)]

def plot_compare():
    plt.figure(1)
    plt.scatter(ra_xmmcut,dec_xmmcut,marker='o',s=10,edgecolor='gray',
                    facecolor='none',label='XMM sources')
    plt.scatter(ra2_em,dec2_em,marker='o',s=52,edgecolor='green',
                    facecolor='none',label='Matched eRASS1 sources')
    plt.scatter(ra2,dec2,marker='+',s=50,color='black',
                    label='eROSITA sources in XMM-obs footprints')
    plt.legend(loc='upper left')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    for i in range(len(id_nomatch)):
        xlimits=[ra2[id_nomatch[i]]-0.25,ra2[id_nomatch[i]]+0.25]
        ylimits=[dec2[id_nomatch[i]]-0.25,dec2[id_nomatch[i]]+0.25]
        plt.xlim(xlimits)
        plt.ylim(ylimits)
        plt.savefig(path+'nomatch_eRASS/'+f'xmm_erass_nomatch_{int(id_nomatch[i])}.pdf')


def plot_matchhist():
    print(len(matches),len(no_matches))
    plt.hist(matches['DET_LIKE_0'],bins=np.logspace(np.log10(np.min(matches['DET_LIKE_0'])),
                                       np.log10(np.max(matches['DET_LIKE_0'])),20),histtype='step',color='red',label='eRASS matches')
    plt.hist(no_matches['DET_LIKE_0'],bins=np.logspace(np.log10(np.min(no_matches['DET_LIKE_0'])),
                                       np.log10(np.max(no_matches['DET_LIKE_0'])),6),histtype='step',color='green',label='eRASS nomatches')
    plt.loglog()
    plt.xlabel('DET_LIKE_n')
    plt.ylabel('Number of Sources')
    plt.savefig(path+'hist_detlike.pdf')
    plt.close()
    
    match_pos2err=matches['POS_ERR']
    nomatch_pos2err=no_matches['POS_ERR']
    plt.hist(match_pos2err,bins=np.logspace(np.log10(np.min(match_pos2err)),
                                       np.log10(np.max(match_pos2err)),20),histtype='step',color='red',label='eRASS matches')
    plt.hist(nomatch_pos2err,bins=np.logspace(np.log10(np.min(nomatch_pos2err)),
                                       np.log10(np.max(nomatch_pos2err)),6),histtype='step',color='green',label='eRASS nomatches')
    plt.loglog()
    plt.xlabel('POS_err(arcsec)')
    plt.ylabel('Number of Sources')
    plt.savefig(path+'hist_POSerr.pdf')
    plt.close()

    match_extlike=matches['EXT_LIKE']
    nomatch_extlike=no_matches['EXT_LIKE']
    print(len(np.where(match_extlike>0)[0]))
    print(len(np.where(nomatch_extlike>0)[0]))

    match_FLUX1=matches['ML_FLUX_1']
    nomatch_FLUX1=no_matches['ML_FLUX_1']
    plt.hist(match_FLUX1,bins=np.logspace(np.log10(np.min(match_FLUX1)),
                                       np.log10(np.max(match_FLUX1)),20),histtype='step',color='red',label='eRASS matches')
    plt.hist(nomatch_FLUX1,bins=np.logspace(np.log10(np.min(nomatch_FLUX1)),
                                       np.log10(np.max(nomatch_FLUX1)),6),histtype='step',color='green',label='eRASS nomatches')
    plt.loglog()
    plt.xlabel('ML_FLUX1')
    plt.ylabel('Number of Sources')
    plt.savefig(path+'hist_MLFLUX.pdf')
    plt.close()

def plot_mindist():
    ra_match=matches['RA']
    dec_match=matches['DEC']
    ra_nomatch=no_matches['RA']
    dec_nomatch=no_matches['DEC']
    coords_match=SkyCoord(ra=ra_match*u.deg,dec=dec_match*u.deg)
    coords_nomatch=SkyCoord(ra=ra_nomatch*u.deg,dec=dec_nomatch*u.deg)
    coords_xmm=SkyCoord(ra=ra_xmmcut*u.deg,dec=dec_xmmcut*u.deg)
    mindist_match=mindist_in_cat(coordinates=coords_match,coordmatch=coords_xmm).value*3600
    mindist_nomatch=mindist_in_cat(coordinates=coords_nomatch,coordmatch=coords_xmm).value*3600
    print(mindist_match,mindist_nomatch)

    plt.hist(mindist_match,bins=np.linspace(np.min(mindist_match),np.max(mindist_match),25),
             histtype='step',color='red',label='eRASS matches')
    plt.hist(mindist_nomatch,bins=np.linspace(10,100,15),
             histtype='step',color='green',label='eRASS nomatches')
    # plt.loglog()
    plt.xlabel('min_distance between eRASS and XMM')
    plt.ylabel('Number of Sources')
    plt.savefig(path+'hist_mindist.pdf')
    plt.close()

if __name__ == '__main__':
    plot_mindist()
    # plot_matchhist()
    # print(unique_elements)
    # print(id_nomatch)
    # print(len(unique_elements),len(id_nomatch))

    # plt.show()
