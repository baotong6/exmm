'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2025-01-13 13:56:42
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2025-01-13 15:06:45
FilePath: /code/latex/add_spec2table.py
Description: 

Copyright (c) 2025 by baotong, All Rights Reserved. 
'''
import pandas as pd


def merge_detector_tables(mos1_path, mos2_path, pn_path, output_path, key_column='source_id'):
    """
    Merge information from three CSV tables (mos1, mos2, pn) into a single table.
    
    The merged table will contain information from pn, mos2, or mos1 based on availability,
    and will include a column indicating the data sources for each row.

    Parameters:
        mos1_path (str): Path to the mos1 CSV file.
        mos2_path (str): Path to the mos2 CSV file.
        pn_path (str): Path to the pn CSV file.
        output_path (str): Path to save the merged table.
        key_column (str): Column used for merging (default is 'source_id').
    """
    # Load the tables
    mos1 = pd.read_csv(mos1_path)
    mos2 = pd.read_csv(mos2_path)
    pn = pd.read_csv(pn_path)
    
    # Add a column to each table to identify the source of the data
    mos1['detector'] = 'mos1'
    mos2['detector'] = 'mos2'
    pn['detector'] = 'pn'
    
    # Concatenate the tables, keeping track of the detector source
    all_data = pd.concat([mos1, mos2, pn], ignore_index=True)
    
    # Create a new DataFrame to hold the merged results
    merged_data = all_data.groupby(key_column).apply(merge_rows).reset_index(drop=True)
    
    # Save the merged table
    merged_data.to_csv(output_path, index=False)

def merge_rows(group):
    """
    Merge rows from the same source_id based on the detector priority: pn > mos2 > mos1.
    
    This function creates a new row that includes:
    - Columns from pn if available, otherwise from mos2, otherwise from mos1.
    - A 'detector_sources' column that records which detectors had data.
    """
    # Determine available detectors for this group
    detectors = '+'.join(group['detector'].unique())
    
    # Select the preferred row based on the detector priority
    if 'pn' in group['detector'].values:
        preferred_row = group[group['detector'] == 'pn'].iloc[0]
    elif 'mos2' in group['detector'].values:
        preferred_row = group[group['detector'] == 'mos2'].iloc[0]
    else:
        preferred_row = group[group['detector'] == 'mos1'].iloc[0]
    
    # Add the detector sources information
    preferred_row['detector_sources'] = detectors
    
    return preferred_row

def merge_tables(table1_path, table2_path, table1_format, table2_format, output_path, columns_to_add):
    """
    Merge specific columns from table2 into table1 based on matching 'source_id' and 'xmm_index'.
    
    Parameters:
        table1_path (str): Path to the first table (table1).
        table2_path (str): Path to the second table (table2).
        table1_format (str): Format of the first table ('xlsx' or 'csv').
        table2_format (str): Format of the second table ('xlsx' or 'csv').
        output_path (str): Path to save the merged table.
        columns_to_add (list): List of columns to add from table2 to table1.
    """
    # Load table1 based on its format
    if table1_format == 'xlsx':
        table1 = pd.read_excel(table1_path)
    elif table1_format == 'csv':
        table1 = pd.read_csv(table1_path)
    else:
        raise ValueError("Unsupported format for table1. Use 'xlsx' or 'csv'.")

    # Load table2 based on its format
    if table2_format == 'xlsx':
        table2 = pd.read_excel(table2_path)
    elif table2_format == 'csv':
        table2 = pd.read_csv(table2_path)
    else:
        raise ValueError("Unsupported format for table2. Use 'xlsx' or 'csv'.")

    # Ensure columns_to_add exist in table2
    for col in columns_to_add:
        if col not in table2.columns:
            raise ValueError(f"Column '{col}' not found in table2.")
    
    # Merge table2 into table1 based on the matching source_id (table2) and xmm_index (table1)
    merged_table = table1.merge(table2[['source_id'] + columns_to_add], 
                                left_on='xmm_index', 
                                right_on='source_id', 
                                how='left')
    
    # Drop the redundant 'source_id' column if needed
    merged_table.drop(columns=['source_id'], inplace=True)
    merged_table = merged_table.astype(str)
    # Save the merged table
    if output_path.endswith('.xlsx'):
        merged_table.to_excel(output_path, index=False)
    elif output_path.endswith('.csv'):
        merged_table.to_csv(output_path, index=False)
    else:
        raise ValueError("Unsupported output format. Use 'xlsx' or 'csv'.")
    
if __name__ == "__main__":
    # Example usage
    path='/Users/baotong/data_GalDisc/data/'
    # merge_detector_tables(path+'combspec/fitresult/'+'fit_mos1_results_tbabs_apec_2sig.csv', 
    #                       path+'combspec/fitresult/'+'fit_mos2_results_tbabs_apec_2sig.csv', 
    #                       path+'combspec/fitresult/'+'fit_pn_results_tbabs_apec_2sig.csv', 
    #                       path+'combspec/fitresult/'+'merged_tbabs_apec_2sig.csv')

    merge_tables(path+'match_e_xmm/'+'matched_gaia_results_optimized_4sec.xlsx', 
                path+'combspec/fitresult/'+'merged_tbabs_apec_2sig.csv', 
                'xlsx', 'csv', path+'match_e_xmm/'+'e_xmmdr14s_merge_spec_starinfo.xlsx', 
                columns_to_add=['detector', 'nh','nh_e1','nh_e2','kT','kT_e1',
                                'kT_e2','abund','abund_e1','abund_e2',
                                'chi2','dof','red_chi2','detector_sources'])
