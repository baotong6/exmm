'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-12-02 12:37:31
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2025-01-23 17:27:32
FilePath: /code/match_eSASS/plot_3catmatch.py
Description: 

Copyright (c) 2025 by baotong, All Rights Reserved. 
'''
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.gridspec as gridspec
from astropy.io import fits
from scipy.optimize import curve_fit
def rayleigh_distribution(r, sigma):
    """Rayleigh distribution function."""
    return (r / sigma**2) * np.exp(-0.5 * r**2 / sigma**2)

def calculate_mse(observed, predicted):
    return np.mean((observed - predicted)**2)

def calculate_chi_square(observed, predicted):
    return np.sum((observed - predicted)**2 / predicted)

def calculate_r_squared(observed, predicted):
    ss_res = np.sum((observed - predicted)**2)
    ss_tot = np.sum((observed - np.mean(observed))**2)
    return 1 - (ss_res / ss_tot)

def plot_gaia_matchsep(save=0,fit_rayleigh=False):
    path = '/Users/baotong/data_GalDisc/data/match_e_xmm/'
    df = pd.read_excel(path+'matched_gaia_results_optimized_20sec.xlsx', sheet_name='Sheet1')
    df2 = pd.read_excel(path+'e_xmmdr14s_match_all_starinfo.xlsx')

    # Task 1: Identify rows with gaia_index > 0 that also appear multiple times
    filtered_df = df[df['gaia_index'] > 0]
    gaia_index_counts = filtered_df['gaia_index'].value_counts()
    duplicate_gaia_indices = gaia_index_counts[gaia_index_counts > 1].index
    duplicate_rows = df[df['gaia_index'].isin(duplicate_gaia_indices)].index.tolist()
    print("Rows with gaia_index > 0 that are duplicated:")
    print(duplicate_rows)

    # Task 2: For rows with hamstar_index > 0, check corresponding gaia_index matches
    matching_rows = df[(df['CTP_ID'] > 0) & (df['CTP_ID'] == df['gaia_index'])].index.tolist()
    dist_hamstar = df['gaia_separation_arcsec'].iloc[matching_rows]
    print(len(dist_hamstar))

    filtered_df = df[df['CTP_ID'] > 0]
    e_id_groups = filtered_df.groupby('e_id')
    no_match_e_ids = []

    for e_id, group in e_id_groups:
        matching_rows = group[group['CTP_ID'] == group['gaia_index']]
        if matching_rows.empty:
            no_match_e_ids.append(e_id)
    
    print("e_id values with CTP_ID > 0 where no rows have matching CTP_ID and gaia_index:")
    print(no_match_e_ids)

    matching_sources = df[(df['CTP_ID'] == df['gaia_index']) & (df['CTP_ID'] > 0) & (df['separation_arcsec'] < 17)]
    ctp_sep = matching_sources['CTP_SEP']
    gaia_sep = matching_sources['gaia_separation_arcsec']
    eidprint = matching_sources['e_id']

    print('Number of XMM matches:', len(df2[df2['separation_arcsec'] < 17]))
    print('Number of Hamstars:', len(df2[df2['CTP_ID'] > 0]))
    print('Number of Hamstars (with XMM matches):', len(df2[(df2['CTP_ID'] > 0) & (df2['separation_arcsec'] < 17)]))
    print('Number of Hamstars (with Gaia and XMM matches):', len(matching_sources))

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(4, 4, figure=fig)

    # Scatter plot
    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_main.scatter(ctp_sep, gaia_sep, label='Sources', color='blue', alpha=0.5)
    ax_main.plot([min(ctp_sep), max(ctp_sep)], [min(ctp_sep), max(ctp_sep)], color='red', linestyle='--', label='y=x')

    mask = gaia_sep > 4
    ax_main.scatter(ctp_sep[mask], gaia_sep[mask], label='gaia_sep > 4', color='red', s=50,marker='o')

    for i in range(len(ctp_sep[mask])):
        ax_main.text(ctp_sep[mask].iloc[i], gaia_sep[mask].iloc[i], str(eidprint[mask].iloc[i]), fontsize=14, color='black', ha='right')

    ax_main.set_xlabel('eRASS1-Gaia Separation (arcsec)', fontsize=16)
    ax_main.set_ylabel('XMM-Gaia Separation (arcsec)', fontsize=16)
    ax_main.tick_params(labelsize=16)
    ax_main.legend()

    ax_top = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    hist_data_top, bins_top, _=ax_top.hist(ctp_sep, bins=np.linspace(0,20,40), color='blue', lw=3,density=1,histtype='step')

    ax_top.set_ylabel('Number of sources', fontsize=16)
    ax_top.set_title('eRASS1-Gaia Separation', fontsize=16)
    ax_top.tick_params(labelsize=16)

    ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax_main)
    hist_data_right, bins_right, _=ax_right.hist(gaia_sep, bins=np.linspace(0,20,50), color='blue', lw=3,density=1,histtype='step', orientation='horizontal')
    ax_right.set_xlabel('Number of sources', fontsize=16)
    ax_right.set_title('XMM-Gaia Separation', fontsize=16)
    ax_right.tick_params(labelsize=16)

    if fit_rayleigh:
        # Bin centers for fitting
        bin_centers_top = (bins_top[:-1] + bins_top[1:]) / 2
        # Fit the Rayleigh distribution
        popttop, _ = curve_fit(rayleigh_distribution, bin_centers_top, hist_data_top, p0=[4.0])
        sigma_fit = popttop[0]

        # Plot the fitted Rayleigh distribution
        fit_ytop = rayleigh_distribution(bin_centers_top, sigma_fit)
        ax_top.plot(bin_centers_top, fit_ytop, 'r*-', label=f'Rayleigh fit: σ={sigma_fit:.2f}')
        ax_top.legend()

        # Fit and plot Rayleigh distribution for ax_right
        bin_centers_right = (bins_right[:-1] + bins_right[1:]) / 2
        popt_right, _ = curve_fit(rayleigh_distribution, bin_centers_right, hist_data_right, p0=[0.2])
        sigma_fit_right = popt_right[0]
        fit_x_right = rayleigh_distribution(bin_centers_right, sigma_fit_right)
        ax_right.plot(fit_x_right, bin_centers_right, 'r*-', label=f'Rayleigh fit: σ={sigma_fit_right:.2f}')
        ax_right.legend()

        # Calculate fit quality metrics
        r_squared_top = calculate_r_squared(hist_data_top, fit_ytop)
        ax_top.text(0.7, 0.2, f'$R^2$: {r_squared_top:.2f}', transform=ax_top.transAxes, fontsize=15, color='red')
        # mse_top = calculate_mse(hist_data_top, fit_ytop)
        # chi_square_top = calculate_chi_square(hist_data_top, fit_ytop)
        # ax_top.text(0.7, 0.9, f'MSE: {mse_top:.2f}', transform=ax_top.transAxes, fontsize=12, color='red')
        # ax_top.text(0.7, 0.85, f'$\chi^2$: {chi_square_top:.2f}', transform=ax_top.transAxes, fontsize=12, color='red')
        r_squared_right = calculate_r_squared(hist_data_right, fit_x_right)
        ax_right.text(0.2, 0.7, f'$R^2$: {r_squared_right:.2f}', transform=ax_right.transAxes, fontsize=15, color='red')

    plt.tight_layout()
    if save:plt.savefig(path+'useful_figs/' +'XMM_GAIA_eROSITA_sep.pdf',bbox_inches='tight', pad_inches=0.1)
    plt.show()


def tonyfig_old(save=0):
    path = '/Users/baotong/data_GalDisc/data/'
    xmmdr14slist = fits.open(path + 'xmmdr14s/GalDisc_4xmmdr14s_new_cleaned.fits')[1].data
    ncontrib = xmmdr14slist['N_CONTRIB']
    realsrc = xmmdr14slist[np.where(ncontrib > 0)]
    path = '/Users/baotong/data_GalDisc/data/match_e_xmm/'
    df = pd.read_excel(path+'matched_gaia_results_optimized_4sec.xlsx', sheet_name='Sheet1')

    # Convert relevant columns to numeric, handling errors due to strings
    numeric_columns = ['separation_arcsec', 'gaia_parallax', 'gaia_gmag', 'gaia_bpmag', 'gaia_rpmag', 'Fx']
    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Filter sources with CTP_ID < 0
    filtered_df = df[(df['CTP_ID'] < 0) & (df['separation_arcsec'] < 17)]
    # Drop rows with NaN in 'e_id' and 'gaia_separation_arcsec' before groupby
    filtered_df = filtered_df.dropna(subset=['e_id', 'gaia_separation_arcsec'])
    filtered_df = filtered_df.loc[filtered_df.groupby('e_id')['gaia_separation_arcsec'].idxmin()]
    print(len(filtered_df))

    # Filter sources with CTP_ID > 0, keep one source per e_id with smallest separation_arcsec
    matching_sources = df[(df['CTP_ID'] > 0) & (df['separation_arcsec'] < 17)]
    # matching_sources = matching_sources.dropna(subset=['e_id', 'gaia_separation_arcsec'])
    ranking_indices = []
    # Group by 'e_id'
    grouped = matching_sources.groupby('e_id')
    for e_id, group in grouped:
        # Sort the group by 'gaia_separation_arcsec'
        sorted_group = group.sort_values(by='gaia_separation_arcsec').reset_index()
        # Check if there's a match where CTP_ID equals gaia_index
        match = sorted_group[sorted_group['CTP_ID'] == sorted_group['gaia_index']]
        if not match.empty:
            rank = match.index[0]  # Assuming only the first match is needed
            ranking_indices.append(rank)
    # Convert the ranking_indices list to a numpy array
    ranking_indices_array = np.array(ranking_indices)
    print(len(ranking_indices_array))
    print(len(np.where(ranking_indices_array==0)[0]))

    matching_sources = df[(df['CTP_ID'] > 0) & (df['CTP_ID'] == df['gaia_index']) & (df['separation_arcsec'] < 17)]
    print(len(matching_sources))

    # Prepare data A for plotting
    filtered_df['xmmFx']=realsrc['EP_FLUX'][filtered_df['xmm_index']]
    filtered_df['BP_RP'] = filtered_df['gaia_bpmag'] - filtered_df['gaia_rpmag']
    filtered_df['distkpc'] = np.abs(1 / filtered_df['gaia_parallax'])
    filtered_df['distpc'] = filtered_df['distkpc'] * 1000
    filtered_df['abs_Gmag'] = filtered_df['gaia_gmag'] - 5 * (np.log10(filtered_df['distpc']) - 1)

    Lsun = 3.846e33 * u.erg / u.s
    msun=-26.7 
    absmsun = 4.66
    # filtered_df['Lopt'] = 10**(0.4 * (absmsun - filtered_df['gaia_gmag'])) * Lsun
    Fopt_filtered=10**(0.4*(msun-np.array(filtered_df['gaia_gmag'])))*Lsun/(4*np.pi*(1*u.au)**2)
    Fopt_filtered=Fopt_filtered.to(u.erg/u.cm**2/u.s).value



    Fx_filtered = filtered_df['xmmFx']
    BP_RP_filtered = filtered_df['BP_RP']
    abs_Gmag_filtered = filtered_df['abs_Gmag']

    # Prepare data for Class B
    matching_sources['xmmFx']=realsrc['EP_FLUX'][matching_sources['xmm_index']]
    matching_sources['BP_RP'] = matching_sources['gaia_bpmag'] - matching_sources['gaia_rpmag']
    matching_sources['distkpc'] = np.abs(1 / matching_sources['gaia_parallax'])
    matching_sources['distpc'] = matching_sources['distkpc'] * 1000
    matching_sources['abs_Gmag'] = matching_sources['gaia_gmag'] - 5 * (np.log10(matching_sources['distpc']) - 1)

    Fopt_matching=10**(0.4*(msun-np.array(matching_sources['gaia_gmag'])))*Lsun/(4*np.pi*(1*u.au)**2)
    Fopt_matching=Fopt_matching.to(u.erg/u.cm**2/u.s).value


    Fx_matching = matching_sources['xmmFx']
    BP_RP_matching = matching_sources['BP_RP']
    abs_Gmag_matching = matching_sources['abs_Gmag']

    plt.figure(1,(6,6))

    hist1,bins1,_=plt.hist(filtered_df['gaia_separation_arcsec'],bins=np.linspace(0,10,25),histtype='step',color='green',lw=3,density=1,label=f'The remaining {len(Fopt_filtered)} XMM matches')
    hist2,bins2,_=plt.hist(matching_sources['gaia_separation_arcsec'],bins=np.linspace(0,10,31),
    histtype='step',color='k',lw=2,density=1,label=f'{len(Fopt_matching)} "good" Hamstars')

    bin1_centers = (bins1[:-1] + bins1[1:]) / 2
    popt1, _ = curve_fit(rayleigh_distribution, bin1_centers, hist1, p0=[1.0])
    sigma_fit = popt1[0]
    fit_ytop = rayleigh_distribution(bin1_centers, sigma_fit)
    r_squared1 = calculate_r_squared(hist1, fit_ytop)
    plt.plot(bin1_centers, fit_ytop, 'g*-', lw=0.5,label=f'Rayleigh fit: σ={sigma_fit:.2f},$R^2$={r_squared1:.2f}')

    bin2_centers = (bins2[:-1] + bins2[1:]) / 2
    popt2, _ = curve_fit(rayleigh_distribution, bin2_centers, hist2, p0=[1.0])
    sigma_fit = popt2[0]
    fit_ytop = rayleigh_distribution(bin2_centers, sigma_fit)
    r_squared2 = calculate_r_squared(hist2, fit_ytop)
    plt.plot(bin2_centers, fit_ytop,'k*-', lw=0.5,label=f'Rayleigh fit: σ={sigma_fit:.2f},$R^2$={r_squared2:.2f}')

    plt.xlabel('XMM-Gaia Separation', fontsize=16)
    plt.ylabel('Normalized', fontsize=16)
    plt.tick_params(labelsize=16)
    plt.legend()
    if save:plt.savefig(path+'useful_figs/' +'Hist_XMM2GAIA_sep.pdf',bbox_inches='tight', pad_inches=0.1)
    plt.show()

    # Plotting Class A and Class B
    plt.figure(figsize=(6, 6))
    plt.scatter(BP_RP_filtered, abs_Gmag_filtered, marker='.', color='green', label='Class A')
    plt.scatter(BP_RP_matching, abs_Gmag_matching, marker='.', color='gray', label='Class B')
    plt.gca().invert_yaxis()
    plt.xlabel('BP-RP', fontsize=16)
    plt.ylabel('Absolute Gmag', fontsize=16)
    plt.tick_params(labelsize=16)
    plt.legend()
    
    plt.show()

    print('Number of no hamstars:',len(Fopt_filtered))
    print('Number of hamstars:',len(Fopt_matching))
    plt.figure(figsize=(10, 6))
    plt.scatter(BP_RP_filtered, Fx_filtered / Fopt_filtered, marker='.', s=100,color='green', label=f'The remaining {len(Fopt_filtered)} XMM matches')
    plt.scatter(BP_RP_matching, Fx_matching / Fopt_matching, marker='.', s=50,color='gray', label=f'{len(Fopt_matching)} "good" Hamstars')
    # plt.scatter(x, y)

    x1 = np.linspace(-0.4, 0.7, 50)
    y1 = 10 ** (x1 - 3.5)
    x2 = np.linspace(0.7, 5, 50)
    y2 = 10 ** (x2 - 3)

    plt.plot(x1, y1, '-', color='blue')
    plt.plot(x2, y2, label='Modified cut in (Rodriguez+, 2024)', color='blue')
    plt.plot([0.7, 0.7], [y1[-1], y2[0]], '--', color='blue')
    plt.plot([1.5, 1.5], [1e-6, 1e2], '--', color='gray')
    plt.plot([-0.3, -0.3], [1e-6, 1e2], '--', color='gray')

    # for i in range(len(BP_RP_filtered)):
    #     plt.text(BP_RP_filtered.iloc[i], Fx_filtered.iloc[i] / Fopt_filtered[i], 
    #              str(filtered_df['e_id'].iloc[i]), fontsize=14, color='red', ha='right')
    
    plt.xlabel('Gaia BP-RP', fontsize=16)
    plt.ylabel(r'$F_x/F_{\rm opt}$', fontsize=16)
    plt.tick_params(labelsize=16)
    plt.legend()
    plt.semilogy()
    if save:plt.savefig(path+'useful_figs/' +'X-RAY_MainSequence_labels.pdf',bbox_inches='tight', pad_inches=0.1)
    plt.show()
def tonyfig(save=0):
    path= '/Users/baotong/data_GalDisc/data/match_e_xmm/'
    tablename='e_xmmdr14s_merge_spec_starinfo.xlsx'
    df=pd.read_excel(path+tablename,sheet_name='label')
    filtered_df = df[(df['sep_exmm'] < 17)]
    filtered_df=filtered_df.dropna(subset=['G_corrected', 'BP_corrected','RP_corrected'])
    Lsun = 3.846e33 * u.erg / u.s
    msun=-26.7 
    absmsun = 4.66
    Fopt=10**(0.4*(msun-np.array(filtered_df['G_corrected'].values)))*Lsun/(4*np.pi*(1*u.au)**2)
    Fopt=Fopt.to(u.erg/u.cm**2/u.s).value
    Fx = filtered_df['xmmflux'].values
    Fx_err=filtered_df['xmmfluxerr'].values
    BP_RP = filtered_df['BP_corrected'].values-filtered_df['RP_corrected'].values
    CTP_SEP = np.array(filtered_df['CTP_SEP'])
    ToType = np.array(filtered_df['ToType'])
    abs_Gmag = filtered_df['G_corrected'].values-5 * (np.log10(filtered_df['distkpc'].values*1000/10))

    indexstar = np.where(CTP_SEP > 0)[0]
    indexnotstar = np.where(CTP_SEP < 0)[0]
    index_hards = np.where(ToType == 'hardS')[0]
    index_CV = np.where(ToType == 'CV')[0]
    index_LMXB = np.where(ToType == 'LMXB')[0]
    index_AGN = np.where(ToType == 'AGN')[0]
    print(indexstar)
    from dustmaps.marshall import MarshallQuery
    from dustmaps.config import config

    # Set the local data storage path for dustmaps
    config['data_dir'] = '/Users/baotong/dustmaps/'  # Adjust this to your local path
    marshallstar = MarshallQuery()
    hamstar = fits.open(path + 'HamStar_eRASS1_Main_Likely_Identifications_v1.1.fits')[1].data
    mask = (
        ~np.isnan(hamstar['PLX']) &  # PLX should not be NaN
        (hamstar['PLX'] > 0) &       # PLX should be positive
        ~np.isnan(hamstar['G']) &    # G should not be NaN
        ~np.isnan(hamstar['BP_RP'])  # BP_RP should not be NaN
    )
    ra=hamstar[mask]['CTP_RA']
    dec=hamstar[mask]['CTP_DEC']
    distkpc=1/hamstar[mask]['PLX']
    coords=SkyCoord(ra=ra,dec=dec,unit=(u.deg, u.deg), 
                    frame='icrs', distance=distkpc*u.kpc)
    ebv_samples = marshallstar(coords)
    R_G, R_BP, R_RP = 10.1154, 12.8461, 7.5513
    A_G = R_G * np.array(ebv_samples)
    A_BP = R_BP * np.array(ebv_samples)
    A_RP = R_RP * np.array(ebv_samples)
    # Correct magnitudes
    corrected_m_g = hamstar[mask]['G'] - A_G
    corrected_bp_rp=hamstar[mask]['BP_RP']-(A_BP-A_RP)
    absG=corrected_m_g-5*(np.log10(distkpc*1000/10))
    plt.scatter(corrected_bp_rp,absG,marker='.', 
            s=5, color='brown',
            label='All sky Hamstars' )
    plt.scatter(BP_RP[indexnotstar],abs_Gmag[indexnotstar],marker='o', 
                s=120, color='green', edgecolor='black', alpha=0.7,
                label='Not stars' )
                #label=f'The remaining {len(indexnotstar)} XMM matches')
    plt.scatter(BP_RP[indexstar], abs_Gmag[indexstar], 
                marker='o', s=60, color='gray', edgecolor='gray', alpha=0.7,
                label='Hamstars')
                #label=f'{len(indexstar)} "good" Hamstars')
    plt.scatter(BP_RP[index_hards], abs_Gmag[index_hards], 
                marker='^', s=150, color='red', edgecolor='black', alpha=0.9, 
                label=f'Hard X-ray sources (kT>3 keV)')
    plt.scatter(BP_RP[index_LMXB], abs_Gmag[index_LMXB], 
                marker='s', s=150, color='cyan', edgecolor='black', alpha=0.9, 
                label=f'LMXB')
    plt.scatter(BP_RP[index_CV], abs_Gmag[index_CV], 
                marker='v', s=150, color='magenta', edgecolor='black', alpha=0.9,
                label=f'CV')
    plt.scatter(BP_RP[index_AGN], abs_Gmag[index_AGN], 
                marker='o', s=150, color='k', edgecolor='black', alpha=0.9,
                label=f'AGN')

    plt.gca().invert_yaxis()
    plt.xlim(-1,7)
    plt.ylim(16,-10)
    plt.xlabel('BP-RP', fontsize=16)
    plt.ylabel('Absolute Gmag', fontsize=16)
    plt.tick_params(labelsize=16)
    plt.legend()
    if save:plt.savefig(path+'useful_figs/' +'absGmag_BP-RP_legend.pdf',bbox_inches='tight', pad_inches=0.1)
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.scatter(BP_RP[indexnotstar], Fx[indexnotstar]/ Fopt[indexnotstar], 
                marker='.', s=120, color='green', edgecolor='black', alpha=0.7,
                label='Not stars' )
                # label=f'The remaining {len(indexnotstar)} XMM matches')
    plt.scatter(BP_RP[indexstar], Fx[indexstar]/ Fopt[indexstar], 
                marker='.', s=50, color='gray', edgecolor='gray', alpha=0.7, 
                label='Hamstars')
                # label=f'{len(indexstar)} "good" Hamstars')
    plt.scatter(BP_RP[index_hards], Fx[index_hards]/ Fopt[index_hards], 
                marker='^', s=150, color='red', edgecolor='black', alpha=0.9,
                label=f'Hard X-ray sources (kT>3 keV)')
    plt.scatter(BP_RP[index_LMXB], Fx[index_LMXB]/ Fopt[index_LMXB], 
                marker='s', s=150, color='cyan', edgecolor='black', alpha=0.9, 
                label=f'LMXB')
    plt.scatter(BP_RP[index_CV], Fx[index_CV]/ Fopt[index_CV], 
                marker='v', s=150, color='magenta', edgecolor='black', alpha=0.9, 
                label=f'CV')
    plt.scatter(BP_RP[index_AGN], Fx[index_AGN]/ Fopt[index_AGN], 
                marker='o', s=150, color='k', edgecolor='black', alpha=0.9, 
                label=f'AGN')

    x1 = np.linspace(-0.4, 5, 100)
    y1 = 10 ** (x1 - 3.5)
    x2 = np.linspace(0.8, 5, 50)
    y2 = 10 ** (x2 - 3)

    plt.plot(x1, y1, '-', color='k')
    plt.plot(x2, y2, '--',color='blue',label='Modified cut in (Rodriguez+, 2024)')
    plt.plot([0.7, 0.8], [10 ** (0.7 - 3.5), 10 ** (0.8 - 3)], '--', color='blue')
    plt.plot([1.5, 1.5], [1e-6, 1e2], '--', color='gray')
    plt.plot([-0.3, -0.3], [1e-6, 1e2], '--', color='gray')

    # for i in range(len(BP_RP_filtered)):
    #     plt.text(BP_RP_filtered.iloc[i], Fx_filtered.iloc[i] / Fopt_filtered[i], 
    #              str(filtered_df['e_id'].iloc[i]), fontsize=14, color='red', ha='right')
    
    plt.xlabel('Gaia BP-RP', fontsize=16)
    plt.ylabel(r'$F_x/F_{\rm opt}$', fontsize=16)
    plt.tick_params(labelsize=16)
    plt.legend()
    plt.semilogy()
    if save:plt.savefig(path+'useful_figs/' +'X-RAY_MainSequence_legend.pdf',bbox_inches='tight', pad_inches=0.1)
    plt.show()

    return None
def plot_XLF(save=0):
    # Define the file path and table name
    path = '/Users/baotong/data_GalDisc/data/match_e_xmm/'
    tablename = 'matched_gaia_results_optimized_4sec.xlsx'
    df = pd.read_excel(path + tablename, sheet_name='ClosestMatches')
    df = df[df['sep_exmm'] < 17]
    
    # Extract necessary columns
    eid = df['e_id'].values
    fx = df['xmmflux'].values
    fx_err = df['xmmfluxerr'].values
    distkpc = df['distkpc'].values
    Ttype = df['Ttype'].values
    
    # Calculate luminosity Lx and its uncertainty Lxerr
    Lx = (fx * u.erg / (u.cm**2 * u.s) * 4 * np.pi * (distkpc * u.kpc)**2).to(u.erg / u.s).value
    Lxu = ((fx + fx_err) * u.erg / (u.cm**2 * u.s) * 4 * np.pi * (distkpc * u.kpc)**2).to(u.erg / u.s).value
    Lxerr = Lxu - Lx
    
    valid_indices = ~np.isnan(Lx) & ~np.isnan(Lxerr)
    
    # Separate indices based on Ttype
    star_indices = np.where((Ttype == 'star') & valid_indices)[0]
    notstar_indices = np.where((Ttype == 'not') & valid_indices)[0]
    # Plotting for star sources
    plt.errorbar(x=eid[star_indices], y=Lx[star_indices], yerr=Lxerr[star_indices],
            fmt='.', capsize=2, elinewidth=1, ecolor='r',color='r', linewidth=1.5)
    plt.errorbar(x=eid[notstar_indices], y=Lx[notstar_indices], yerr=Lxerr[notstar_indices],
            fmt='.', capsize=2, elinewidth=1, ecolor='green',color='green', linewidth=1.5)

    plt.xlabel('e_id')
    plt.ylabel(r'$L_X$ (erg/s)')
    plt.semilogy()  # Use logarithmic scale for y-axis
    plt.title('X-ray Luminosity Function (XLF) for Stars')
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.tight_layout()
    if save:
        plt.savefig(path+'useful_figs/' +'XMMflux_scatter.pdf',bbox_inches='tight', pad_inches=0.1)
    plt.show()

    Lx_sorted1 = Lx[star_indices]
    Lx_sorted2 = Lx[notstar_indices]
    # Create the weights for the histogram to reflect N(>Lx)
    # weights1 = np.ones_like(Lx_sorted1) / len(Lx_sorted1)
    # weights2 = np.ones_like(Lx_sorted2) / len(Lx_sorted2)
    plt.figure(figsize=(10, 6))
    plt.hist(Lx_sorted1, bins=np.logspace(27,35,18),  
             lw=2,histtype='step',color='red', cumulative=-1,label='Stars')
    
    plt.hist(Lx_sorted2, bins=np.logspace(29,37,18), 
         lw=3, histtype='step', color='green', cumulative=-1,label='Not Stars')

    # Configure histogram plot
    plt.xlabel(r'$L_X$ (erg/s)',fontsize=16)
    plt.ylabel('N (>L)',fontsize=16)
    plt.loglog()  # Use logarithmic scale for x-axis
    plt.title('Luminosity Function',fontsize=16)
    plt.legend()
    # plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.tight_layout()
    plt.tick_params(labelsize=16)
    if save:
        plt.savefig(path+'useful_figs/' +'XLF.pdf',bbox_inches='tight', pad_inches=0.1)
    plt.show()

def plot_histLF(save=0):
    # Define the file path and table name
    path = '/Users/baotong/data_GalDisc/data/match_e_xmm/'
    tablename = 'matched_gaia_results_optimized_4sec.xlsx'

    df = pd.read_excel(path + tablename, sheet_name='ClosestMatches')
    df = df[df['sep_exmm'] < 17]

    fx = df['xmmflux'].values
    fx_err = df['xmmfluxerr'].values
    distkpc = df['distkpc'].values
    Ttype = df['Ttype'].values
    
    # Calculate luminosity Lx
    Lx = (fx * u.erg / (u.cm**2 * u.s) * 4 * np.pi * (distkpc * u.kpc)**2).to(u.erg / u.s).value
    
    valid_indices = ~np.isnan(Lx)
    star_indices = np.where((Ttype == 'star') & valid_indices)[0]
    notstar_indices = np.where((Ttype == 'not') & valid_indices)[0]

    # Bin edges for histogram
    bins = np.linspace(27, 37, 11)  # Adjust these based on your data range
    counts_star, bin_edges_star = np.histogram(np.log10(Lx[star_indices]), bins=bins)
    bin_centers_star = (bin_edges_star[:-1] + bin_edges_star[1:]) / 2
    bin_widths_star = bin_edges_star[1:] - bin_edges_star[:-1]
    errors_star = np.sqrt(counts_star)
    
    counts_notstar, bin_edges_notstar = np.histogram(np.log10(Lx[notstar_indices]), bins=bins)
    bin_centers_notstar = (bin_edges_notstar[:-1] + bin_edges_notstar[1:]) / 2
    bin_widths_notstar = bin_edges_notstar[1:] - bin_edges_notstar[:-1]
    errors_notstar = np.sqrt(counts_notstar)

    plt.figure(figsize=(10, 6))
    
    # Plot histograms for visual comparison
    plt.hist(np.log10(Lx[star_indices]), bins=bins, alpha=0.7, color='gray', label='Stars')
    plt.hist(np.log10(Lx[notstar_indices]), bins=bins, alpha=0.9, color='green', label='Not Stars')
    
    # Plot error bars
    plt.errorbar(bin_centers_star, counts_star, yerr=errors_star, xerr=bin_widths_star / 2,
                 fmt='o', color='gray', capsize=3)
    plt.errorbar(bin_centers_notstar, counts_notstar, yerr=errors_notstar, xerr=bin_widths_notstar / 2,
                 fmt='o', color='green', capsize=3)
    
    # Logarithmic scale for x and y
    plt.yscale('log')
    plt.xlabel(r'$log_{10}(L_X)$ (erg/s)', fontsize=16)
    plt.ylabel('Number of Sources', fontsize=16)
    # plt.title('Luminosity Function', fontsize=16)
    plt.legend()
    plt.tight_layout()
    plt.tick_params(labelsize=16)
    
    if save:
        plt.savefig(path + 'useful_figs/' + 'histLF.pdf', bbox_inches='tight', pad_inches=0.1)
    
    plt.show()

if __name__ == '__main__':
    # plot_gaia_matchsep(save=1,fit_rayleigh=1)
    tonyfig(save=1)
    # plot_XLF(1)
    # plot_histLF(save=1)
