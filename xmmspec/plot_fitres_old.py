import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_comparison_with_error(csv_file, source_groups, parameter, group_labels, colors, red_chi2_threshold=None, output_plot="scatter_comparison.png"):
    """
    Plot a scatter plot with error bars for a given parameter from a CSV file, comparing multiple source groups.
    
    Parameters:
        csv_file (str): Path to the CSV file.
        source_groups (list of lists): List of source ID lists for different groups to compare.
        parameter (str): The parameter to plot (e.g., 'kT', 'nh', 'abund').
        group_labels (list of str): Labels for each source group.
        colors (list of str): Colors to use for each group.
        red_chi2_threshold (float, optional): Threshold to filter data based on red_chi2 values (default: None).
        output_plot (str): Path to save the output plot (default: "scatter_comparison.png").
    
    Returns:
        None
    """
    # Load the CSV file
    data = pd.read_csv(csv_file)

    # Ensure required columns exist
    param_columns = [parameter, f"{parameter}_e1", f"{parameter}_e2", "red_chi2"]
    for col in param_columns:
        if col not in data.columns:
            raise ValueError(f"Column '{col}' is missing in the CSV file.")
    
    # Apply chi2 filtering if a threshold is specified
    if red_chi2_threshold is not None:
        data = data[data['red_chi2'] < red_chi2_threshold]
        if data.empty:
            raise ValueError(f"No data points with red_chi2 < {red_chi2_threshold}.")

    plt.figure(figsize=(10, 6))

    for group, label, color in zip(source_groups, group_labels, colors):
        # Filter data by the current group of source IDs
        group_data = data[data['source_id'].isin(group)]

        # Extract parameter values and errors
        param_values = group_data[parameter].values
        param_errors_low = group_data[f"{parameter}_e1"].values
        param_errors_high = group_data[f"{parameter}_e2"].values
        source_labels = group_data['source_id'].values

        # Calculate symmetric error bars for visualization
        param_errors = np.array([param_errors_low, param_errors_high])

        # Plot scatter with error bars for the group
        plt.errorbar(
            x=source_labels, 
            y=param_values, 
            yerr=param_errors, 
            fmt="o", 
            color=color,
            ecolor=color, 
            capsize=5, 
            label=label
        )

    plt.xlabel("Source ID")
    plt.ylabel(parameter)
    # plt.semilogy()
    plt.title(f"Comparison of {parameter} Across Groups (red_chi2 < {red_chi2_threshold})")
    plt.xticks(rotation=45)
    plt.legend()
    plt.tight_layout()
    # plt.savefig(output_plot)
    plt.show()


def plot_comparison_with_custom_x(csv_file, source_groups, x_parameter, y_parameter, group_labels, colors, red_chi2_threshold=None, output_plot="scatter_comparison_custom_x.png"):
    """
    Plot a scatter plot with error bars for a given parameter, using a custom x-axis parameter from a CSV file.
    
    Parameters:
        csv_file (str): Path to the CSV file.
        source_groups (list of lists): List of source ID lists for different groups to compare.
        x_parameter (str): The parameter to use for the x-axis (e.g., 'red_chi2', 'nh').
        y_parameter (str): The parameter to plot on the y-axis (e.g., 'kT', 'nh', 'abund').
        group_labels (list of str): Labels for each source group.
        colors (list of str): Colors to use for each group.
        red_chi2_threshold (float, optional): Threshold to filter data based on red_chi2 values (default: None).
        output_plot (str): Path to save the output plot (default: "scatter_comparison_custom_x.png").
    
    Returns:
        None
    """
    # Load the CSV file
    data = pd.read_csv(csv_file)

    # Ensure required columns exist
    required_columns = [x_parameter, y_parameter, f"{y_parameter}_e1", f"{y_parameter}_e2", "red_chi2"]
    for col in required_columns:
        if col not in data.columns:
            raise ValueError(f"Column '{col}' is missing in the CSV file.")
    
    # Apply chi2 filtering if a threshold is specified
    if red_chi2_threshold is not None:
        data = data[data['red_chi2'] < red_chi2_threshold]
        if data.empty:
            raise ValueError(f"No data points with red_chi2 < {red_chi2_threshold}.")

    plt.figure(figsize=(10, 6))

    for group, label, color in zip(source_groups, group_labels, colors):
        # Filter data by the current group of source IDs
        group_data = data[data['source_id'].isin(group)]

        # Extract x-axis and y-axis values and errors
        x_values = group_data[x_parameter].values
        y_values = group_data[y_parameter].values
        y_errors_low = group_data[f"{y_parameter}_e1"].values
        y_errors_high = group_data[f"{y_parameter}_e2"].values

        # Calculate symmetric error bars for visualization
        y_errors = np.array([y_errors_low, y_errors_high])

        # Plot scatter with error bars for the group
        plt.errorbar(
            x=x_values, 
            y=y_values, 
            yerr=y_errors, 
            fmt="o", 
            color=color,
            ecolor=color, 
            capsize=5, 
            label=label
        )

    plt.xlabel(x_parameter)
    plt.ylabel(y_parameter)
    plt.title(f"Comparison of {y_parameter} vs {x_parameter} Across Groups (red_chi2 < {red_chi2_threshold})")
    plt.legend()
    plt.tight_layout()
    plt.semilogy()
    # Save and show the plot
    # plt.savefig(output_plot)
    plt.show()


def plot_histogram_distributions(csv_file, source_groups, parameter, group_labels, colors, bins=30, log_scale=False, output_plot="histogram_distribution.png"):
    """
    Plot histograms of a given parameter for multiple source groups.
    
    Parameters:
        csv_file (str): Path to the CSV file.
        source_groups (list of lists): List of source ID lists for different groups to compare.
        parameter (str): The parameter to plot (e.g., 'kT', 'nh', 'abund').
        group_labels (list of str): Labels for each source group.
        colors (list of str): Colors to use for each group.
        bins (int): Number of bins for the histogram (default: 30).
        log_scale (bool): Whether to use a logarithmic scale for the y-axis (default: False).
        output_plot (str): Path to save the output plot (default: "histogram_distribution.png").
    
    Returns:
        None
    """
    # Load the CSV file
    data = pd.read_csv(csv_file)

    # Ensure required column exists
    if parameter not in data.columns:
        raise ValueError(f"Column '{parameter}' is missing in the CSV file.")

    plt.figure(figsize=(10, 6))

    for group, label, color in zip(source_groups, group_labels, colors):
        # Filter data by the current group of source IDs
        group_data = data[data['source_id'].isin(group)]

        # Extract parameter values
        param_values = group_data[parameter].dropna().values

        # Plot histogram
        plt.hist(
            param_values,
            bins=bins,
            histtype='step',
            color=color,
            label=label,
            density=0
        )

    plt.xlabel(parameter)
    plt.ylabel("Density")
    plt.title(f"Histogram of {parameter} for Different Source Groups")
    # plt.xscale("log")
    if log_scale:
        plt.yscale("log")
    plt.legend()
    plt.tight_layout()

    # Save and show the plot
    # plt.savefig(output_plot)
    plt.show()




if __name__ == "__main__":
    # Example usage
    path = '/Users/baotong/data_GalDisc/data/combspec/fitresult/'

    # Load matching data
    ematch = pd.read_excel('/Users/baotong/data_GalDisc/data/match_e_xmm/e_xmmdr14s_match_all.xlsx', sheet_name='all')
    sep = ematch['separation_arcsec'].values
    cand_ematch_index = np.where(sep < 17)[0]
    starindex = np.where(ematch['hamstar'].values > 0)[0]
    notstarindex = np.where(ematch['hamstar'].values < 0)[0]

    # Define source IDs for each group
    starid = ematch['xmm_index'].values[np.intersect1d(starindex, cand_ematch_index)].astype(int)
    notstarid = ematch['xmm_index'].values[np.intersect1d(notstarindex, cand_ematch_index)].astype(int)

    # CSV file and parameter to compare
    csv_file = path + "fit_mos2_results_tbabs_apec_2sig.csv"
    parameter = "kT1"

    # Call the function to plot
    plot_comparison_with_error(
        csv_file=csv_file,
        source_groups=[starid, notstarid],
        parameter=parameter,
        group_labels=["Star-like sources", "Non-star sources"],
        colors=["blue", "orange"],
        red_chi2_threshold=1.5,
        output_plot="kT_comparison_scatter_chi2.png"
    )

    x_parameter = "kT2"  # Replace with the x-axis parameter you want
    y_parameter = "kT1"  # Replace with the y-axis parameter you want

    # Call the function to plot
    plot_comparison_with_custom_x(
        csv_file=csv_file,
        source_groups=[starid, notstarid],
        x_parameter=x_parameter,
        y_parameter=y_parameter,
        group_labels=["Star-like sources", "Non-star sources"],
        colors=["blue", "orange"],
        red_chi2_threshold=1.2,
        output_plot="comparison_scatter_nh_vs_kT.png"
    )

    parameter = "red_chi2"

    # Call the histogram function
    plot_histogram_distributions(
        csv_file=csv_file,
        source_groups=[starid, notstarid],
        parameter=parameter,
        group_labels=["Star-like sources", "Non-star sources"],
        colors=["blue", "orange"],
        bins=np.linspace(0,3,20),
        log_scale=False,
        output_plot="kT1_histogram.png"
    )
