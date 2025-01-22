import numpy as np
import pandas as pd
from astropy.io import fits

def get_hamstar_columns(fits_file):
    """
    Retrieve column names from the hamstar.fit file.
    """
    with fits.open(fits_file) as hdul:
        hamstar_columns = hdul[1].columns.names
    return hamstar_columns

def write_hamstar(maincat_file, ematchinfo_file, fits_file, output_xlsx, selected_columns=None):
    """
    Write matched hamstarindex, CTP_ID, and other custom columns into an xlsx file.
    Columns to write can be customized with selected_columns.
    """
    # Load hamstar, main catalog, and match information data
    hamstar = fits.open(fits_file)[1].data
    maincat = fits.open(maincat_file)[1].data
    ematchinfo = pd.read_csv(ematchinfo_file)
    
    # Initialize DataFrame from maincat data
    df = pd.DataFrame({
        'IAUNAME': maincat['IAUNAME'],
        'e_id': ematchinfo['e_id'],
        'e_ra': ematchinfo['e_ra'],
        'e_dec': ematchinfo['e_dec'],
        'xmm_index': ematchinfo['xmm_index'],
        'xmm_ra': ematchinfo['xmm_ra'],
        'xmm_dec': ematchinfo['xmm_dec'],
        'separation_arcsec': ematchinfo['separation_arcsec']
    })
    
    # Initialize new columns
    df['hamstarindex'] = '-1'
    df['CTP_ID'] = '-1'  # Ensure CTP_ID is initialized as string
    
    # Check if specific columns are selected, otherwise use default columns
    if selected_columns is None:
        selected_columns = ['RAJ2000', 'DEJ2000', 'ERO', 'p-coronal']  # Default columns
    
    # Initialize selected columns with '-1' (as string)
    for col in selected_columns:
        if col not in ['hamstarindex', 'CTP_ID']:  # Exclude already initialized columns
            df[col] = '-1'
    
    # Match sources
    for i, row in df.iterrows():
        match = np.where(hamstar['ERO_IAUNAME'] == row['IAUNAME'])[0]
        if len(match) == 1:
            match_index = match[0]
            df.at[i, 'hamstarindex'] = str(match_index)
            df.at[i, 'CTP_ID'] = str(hamstar['CTP_ID'][match_index])  # Ensure CTP_ID is string
            for col in selected_columns:
                df.at[i, col] = str(hamstar[col][match_index])
    
    # Write to a new xlsx file, ensuring no type conversion issues
    df.to_excel(output_xlsx, index=False)
    print(f"Data has been written to {output_xlsx}")

# Retrieve all column names from hamstar.fit
hamstar_file = '/Users/baotong/data_GalDisc/data/match_e_xmm/HamStar_eRASS1_Main_Likely_Identifications_v1.1.fits'
hamstar_columns = get_hamstar_columns(hamstar_file)
print("Columns in Hamstar file:", hamstar_columns)

# Example usage
maincat_file = '/Users/baotong/data_GalDisc/data/match_e_xmm/eRASS1_filtered.fits'  # Main catalog file
ematchinfo_file = '/Users/baotong/data_GalDisc/data/match_e_xmm/e_xmmdr14s_match_all.csv'  # Match info CSV file
output_xlsx = '/Users/baotong/data_GalDisc/data/match_e_xmm/e_xmmdr14s_match_all_starinfo.xlsx'  # Output xlsx file path
selected_columns = ['p_coronal', 'CTP_SEP', 'Fx', 'G', 'BP_RP', 'PLX', 'SIMBAD_NAME', 'SIMBAD_OTYPE',
                    'CORONAL', 'OB_STAR', 'FLAG_OPT']  # Select columns to write
write_hamstar(maincat_file, ematchinfo_file, hamstar_file, output_xlsx, selected_columns)
