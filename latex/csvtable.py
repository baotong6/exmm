'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-12-19 15:01:05
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2025-01-13 16:44:18
FilePath: /code/latex/csvtable.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
import pandas as pd
import os

def read_file(input_file,sheetname=None):
    """
    Reads the input file, which can be in CSV or Excel format.

    :param input_file: Path to the input file.
    :return: pandas DataFrame.
    """
    file_extension = os.path.splitext(input_file)[1].lower()
    if file_extension == '.csv':
        return pd.read_csv(input_file)
    elif file_extension == '.xlsx':
        if sheetname:
            return pd.read_excel(input_file,sheet_name=sheetname)
        else:return pd.read_excel(input_file)
    else:
        raise ValueError("Unsupported file format. Please provide a .csv or .xlsx file.")
def add_latex_column(df, lower_col, upper_col, base_col, new_col_name, precision=2):
    """
    Adds a new column to the DataFrame with LaTeX formatted values.

    :param df: pandas DataFrame.
    :param lower_col: Column name for the lower bound values.
    :param upper_col: Column name for the upper bound values.
    :param base_col: Column name for the base values.
    :param new_col_name: Name for the new column to be added.
    :param precision: Number of decimal places to round the values.
    :return: DataFrame with the new LaTeX formatted column.
    """
    def format_latex(base, lower, upper):
        base = round(base, precision)
        lower = round(lower, precision)
        upper = round(upper, precision)
        return f"${base}_{{-{lower}}}^{{+{upper}}}$"

    df[new_col_name] = df.apply(lambda row: format_latex(row[base_col], row[lower_col], row[upper_col]), axis=1)
    return df
def csv_to_latex(input_file, sheetname,output_file, selected_rows=None, selected_columns=None, precision=None):
    """
    Converts a CSV or Excel file to LaTeX table code with customizable options.

    :param input_file: Path to the input file, can be CSV or Excel format.
    :param output_file: Path to the output LaTeX file.
    :param selected_rows: Lambda function to filter rows, e.g., lambda x: x['source_id'] > 500.
    :param selected_columns: List of columns to retain, e.g., ['nh', 'kT1', 'abund'].
    :param precision: Dictionary specifying decimal precision for each column, e.g., {'nh': 2, 'kT1': 3}.
    """
    # Read the input file
    df = read_file(input_file,sheetname)
    df = add_latex_column(df, 'kT_e1', 'kT_e2', 'kT', 'kT_latex', precision=2)
    # Filter rows
    if selected_rows is not None:
        df = df[df.apply(selected_rows, axis=1)]

    # Retain only the selected columns
    if selected_columns is not None:
        df = df[selected_columns]

    # Format numerical values based on precision
    if precision is not None:
        for col, prec in precision.items():
            if col in df.columns:
                df[col] = df[col].map(lambda x: f"{x:.{prec}f}")

    # Replace underscores in column names for LaTeX compatibility
    df.columns = [col.replace('_', '\\_') for col in df.columns]
    selected_columns = [col.replace('_', '\\_') for col in selected_columns]

    # Construct LaTeX table
    latex_code = r"""
\documentclass[a4paper,10pt]{article}
\usepackage{booktabs}
\usepackage{longtable}
\begin{document}
\section*{Fit Results}
\begin{longtable}{%s}
\toprule
%s \\
\midrule
\endfirsthead
\toprule
%s \\
\midrule
\endhead
\bottomrule
\endfoot
""" % (
        "c" * len(selected_columns),
        " & ".join(selected_columns),
        " & ".join(selected_columns),
    )

    for _, row in df.iterrows():
        row_data = [str(row[col]) for col in df.columns]
        latex_code += " & ".join(row_data) + r" \\" + "\n"

    latex_code += """
\end{longtable}
\end{document}
"""

    # Write to output file
    with open(output_file, "w") as f:
        f.write(latex_code)

# Example usage
if __name__ == "__main__":
    ## spectra
    # path = '/Users/baotong/data_GalDisc/data/combspec/fitresult/'
    # csv_to_latex(
    #     input_file=path+"fit_mos2_results_tbabs_2apec_2sig.csv",
    #     output_file=path+'texfile/'+"fit_mos2_results_tbabs_2apec_2sig.tex",
    #     selected_rows=lambda x: x['source_id'] > 0,
    #     selected_columns=['source_id', 'nh', 'kT1', 'kT2', 'abund', 'red_chi2'],
    #     precision={'nh': 2, 'kT1': 2, 'kT2': 2, 'abund': 2, 'red_chi2': 2}
    # )
    ## starinfo
    path2 = '/Users/baotong/data_GalDisc/data/match_e_xmm/'
    df = read_file(path2+"e_xmmdr14s_merge_spec_starinfo.xlsx",sheetname='Sheet1')
    csv_to_latex(
        input_file=path2+"e_xmmdr14s_merge_spec_starinfo.xlsx",sheetname='Sheet1',
        output_file=path2+'texfile/'+"output_table.tex",
        selected_rows=lambda x: x['e_id'] >= 0,
        selected_columns=['e_id','xmm_index','xmm_ra','xmm_dec',
                          'sep_exmm','CTP_SEP','sep_xgaia',
                          'distkpc','SIMBAD_OTYPE','kT_latex','red_chi2'],
        precision={'xmm_ra':5,'xmm_dec':5,'sep_exmm':2,'CTP_SEP':2,'sep_xgaia':2,'distkpc':2,'red_chi2':2}
    )
    
