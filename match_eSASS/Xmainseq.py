from dustmaps.marshall import MarshallQuery
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd
from dustmaps.config import config

# Set the local data storage path for dustmaps
config['data_dir'] = '/Users/baotong/dustmaps/'  # Adjust this to your local path

def read_and_filter_sources(file_path, ra_col, dec_col, dist_col, mag_cols):
    """
    Reads an Excel file and returns the DataFrame, keeping all rows but marking invalid ones with NaN in the new columns.

    Parameters:
        file_path (str): Path to the input Excel file.
        ra_col (str): Name of the RA column.
        dec_col (str): Name of the Dec column.
        dist_col (str): Name of the distance column.
        mag_cols (list): List of column names for observed magnitudes (G, BP, RP).

    Returns:
        pd.DataFrame: DataFrame with all rows, adding NaN where data is missing.
    """
    # Read the input Excel file
    df = pd.read_excel(file_path, sheet_name='label')
    
    # Ensure all required columns are present
    required_cols = [ra_col, dec_col, dist_col] + mag_cols
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    
    return df

def correct_magnitudes(input_table, ra_col, dec_col, dist_col, g_mag_col, bp_mag_col, rp_mag_col):
    """
    Corrects Gaia G, BP, and RP magnitudes for extinction using the Bayestar extinction model.
    Keeps all rows from the original table and adds NaN for rows with missing values.

    Parameters:
        input_table (pd.DataFrame): Input table containing source coordinates and magnitudes.
        ra_col (str): Name of the RA column in the table.
        dec_col (str): Name of the Dec column in the table.
        dist_col (str): Name of the distance column (in kpc).
        g_mag_col (str): Name of the G-band observed magnitude column.
        bp_mag_col (str): Name of the BP-band observed magnitude column.
        rp_mag_col (str): Name of the RP-band observed magnitude column.

    Returns:
        pd.DataFrame: A new DataFrame with added columns for corrected G, BP, and RP magnitudes.
    """
    # Extract the required columns
    required_columns = [ra_col, dec_col, dist_col, g_mag_col, bp_mag_col, rp_mag_col]
    for col in required_columns:
        if col not in input_table.columns:
            raise ValueError(f"Missing required column: {col}")
    
    # Get values from input_table
    ra = input_table[ra_col].values
    dec = input_table[dec_col].values
    distances = input_table[dist_col].values
    m_g = input_table[g_mag_col].values
    m_bp = input_table[bp_mag_col].values
    m_rp = input_table[rp_mag_col].values

    # Create a mask for valid values
    mask = ~np.isnan(ra) & ~np.isnan(dec) & ~np.isnan(distances) & ~np.isnan(m_g) & ~np.isnan(m_bp) & ~np.isnan(m_rp)
    
    # Create SkyCoord object only for valid rows
    coords = SkyCoord(ra=ra[mask], dec=dec[mask], unit=(u.deg, u.deg), frame='icrs', distance=distances[mask]*u.kpc)
    marshallstar = MarshallQuery()
    
    try:
        ebv_samples = marshallstar(coords)  # Returns an array of E(B-V) samples
    except Exception as e:
        raise RuntimeError(f"Error querying Bayestar extinction model: {e}")
    
    # Extinction correction coefficients for G, BP, and RP bands
    R_G, R_BP, R_RP = 10.1154, 12.8461, 7.5513
    A_G = R_G * np.array(ebv_samples)
    A_BP = R_BP * np.array(ebv_samples)
    A_RP = R_RP * np.array(ebv_samples)

    # Correct magnitudes
    corrected_m_g = m_g[mask] - A_G
    corrected_m_bp = m_bp[mask] - A_BP
    corrected_m_rp = m_rp[mask] - A_RP

    # Create a copy of the original table to preserve all rows
    output_table = input_table.copy()

    # Apply the corrections only where data was valid
    output_table.loc[mask, 'A_G'] = A_G
    output_table.loc[mask, 'A_BP'] = A_BP
    output_table.loc[mask, 'A_RP'] = A_RP
    output_table.loc[mask, 'G_corrected'] = corrected_m_g
    output_table.loc[mask, 'BP_corrected'] = corrected_m_bp
    output_table.loc[mask, 'RP_corrected'] = corrected_m_rp

    # For rows where data was invalid (NaN), these new columns will stay NaN
    return output_table

def process_and_save(file_path, output_path, ra_col, dec_col, dist_col, g_mag_col, bp_mag_col, rp_mag_col):
    """
    Reads an input Excel file, filters sources, corrects magnitudes for extinction,
    and saves the results to a new Excel file, while keeping the original number of rows.

    Parameters:
        file_path (str): Path to the input Excel file.
        output_path (str): Path to save the output Excel file.
        ra_col (str): Name of the RA column.
        dec_col (str): Name of the Dec column.
        dist_col (str): Name of the distance column.
        g_mag_col (str): Name of the G-band observed magnitude column.
        bp_mag_col (str): Name of the BP-band observed magnitude column.
        rp_mag_col (str): Name of the RP-band observed magnitude column.
    """
    # Read the input file
    input_df = read_and_filter_sources(file_path, ra_col, dec_col, dist_col, 
                                       mag_cols=[g_mag_col, bp_mag_col, rp_mag_col])
    
    # Correct magnitudes for extinction
    corrected_df = correct_magnitudes(input_df, ra_col, dec_col, dist_col, 
                                      g_mag_col, bp_mag_col, rp_mag_col)
    
    # Save the corrected data to a new Excel file
    corrected_df.to_excel(output_path, index=False)
    print(f"Corrected magnitudes saved to {output_path}")

# Example usage:
if __name__ == "__main__":
    input_file = '/Users/baotong/data_GalDisc/data/match_e_xmm/e_xmmdr14s_merge_spec_starinfo.xlsx'  # Replace with your input file path
    output_file = '/Users/baotong/data_GalDisc/data/match_e_xmm/output_e_xmmdr14s_merge_spec_starinfo_correctmag.xlsx'  # Replace with your desired output file path
    # Column names in the input Excel file
    ra_column = 'xmm_ra'
    dec_column = 'xmm_dec'
    distance_column = 'distkpc'
    g_mag_column = 'gaia_gmag'
    bp_mag_column = 'gaia_bpmag'
    rp_mag_column = 'gaia_rpmag'
    
    # Process the file and save results
    process_and_save(input_file, output_file, ra_col=ra_column, dec_col=dec_column, 
                     dist_col=distance_column, g_mag_col=g_mag_column, 
                     bp_mag_col=bp_mag_column, rp_mag_col=rp_mag_column)
