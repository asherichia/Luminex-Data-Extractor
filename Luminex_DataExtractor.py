#!/usr/bin/env python3
# Improved Luminex GUI - With standard efficiency calculation, sample flagging, and hook effect detection

import pandas as pd
import numpy as np
import os
import datetime
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# 4-parameter logistic on log10(concentration): x = log10(conc)
def four_pl(x, A, B, C, D):
    # y = D + (A-D)/(1 + 10^((C - x)*B))
    return D + (A - D) / (1 + 10 ** ((C - x) * B))

# Invert 4PL to compute log10 concentration from MFI
def invert_log_four_pl(y, A, B, C, D):
    # ratio = (A-D)/(y-D) - 1
    ratio = (A - D) / (y - D) - 1
    # log10 concentration = C - (1/B) * log10(ratio)
    with np.errstate(divide='ignore', invalid='ignore'):
        logx = C - (1.0 / B) * np.log10(ratio)
    return logx

def extract_datetime_from_file(df):
    """Extract datetime from row 3 (index 2), column 0"""
    raw_run_info = df.iloc[2, 0]
    try:
        # Try to parse the datetime
        run_dt = pd.to_datetime(raw_run_info)
        return run_dt
    except Exception:
        # Return current datetime if parsing fails
        return datetime.datetime.now()

def preprocess_csv(input_path, output_prefix=None, plate_format="96-well", is_full_plate=True, num_wells=None):
    """
    Extract the relevant data from a Luminex CSV file, removing metadata
    at the top and bottom of the file. Handles different plate formats.
    
    Args:
        input_path: Path to input CSV file
        output_prefix: Prefix for output files
        plate_format: "96-well" or "384-well"
        is_full_plate: Whether the plate is fully utilized
        num_wells: Number of wells used if not a full plate
    """
    try:
        # Read entire file without header
        df = pd.read_csv(input_path, header=None)
        
        # Extract run timestamp
        run_dt = extract_datetime_from_file(df)
        ts = run_dt.strftime('%Y-%m-%d_%H-%M-%S')
        
        # Set row indices based on plate format and utilization
        meta_top_rows = 51  # This is typically fixed
        
        # Determine data_end_row based on plate format and utilization
        if is_full_plate:
            if plate_format == "96-well":
                data_end_row = 148  # 51 + 96 + 1 header row
            elif plate_format == "384-well":
                data_end_row = 436  # 51 + 384 + 1 header row
            else:
                raise ValueError(f"Unknown plate format: {plate_format}")
        else:
            # For partial plates, calculate based on the specified number of wells
            if num_wells is None or num_wells <= 0:
                raise ValueError("Number of wells must be specified and positive for partial plates")
                
            data_end_row = meta_top_rows + num_wells + 1  # +1 for header row
            print(f"Partial plate with {num_wells} wells: data_end_row set to {data_end_row}")
        
        # Split out metadata
        meta_top = df.iloc[:meta_top_rows]  # First N rows are metadata
        meta_bottom = df.iloc[data_end_row+1:]  # Rows after data_end_row are metadata
        unused = pd.concat([meta_top, meta_bottom], ignore_index=True)
        
        # Create output directory and determine where to save unused data
        output_dir = None
        if output_prefix:
            output_dir = f"{output_prefix}_graphs"
            # Ensure the directory exists
            os.makedirs(output_dir, exist_ok=True)
            # Save in the output directory
            unused_filename = os.path.join(output_dir, f"Unused_data_run_{ts}.csv")
        else:
            # Save in the current directory if no prefix is provided
            unused_filename = f"Unused_data_run_{ts}.csv"
        
        # Save unused metadata to file
        unused.to_csv(unused_filename, index=False, header=False)
        print(f"Saved unused data to: {unused_filename}")  # Debug output
        
        # Main data lies between meta_top_rows and data_end_row (inclusive)
        df_main = df.iloc[meta_top_rows:data_end_row+1].reset_index(drop=True)
        
        # First row in that block is header
        df_main.columns = df_main.iloc[0]
        df_main = df_main.drop(0).reset_index(drop=True)
        
        # Make sure all numeric columns are converted to float
        for col in df_main.columns:
            if col not in ['Sample', 'Batch', 'Bead ID', 'Outlier', 'Location']:  # Skip non-numeric columns
                try:
                    df_main[col] = pd.to_numeric(df_main[col], errors='coerce')
                except:
                    pass
        
        return df_main, unused_filename, ts, output_dir, plate_format
        
    except Exception as e:
        print(f"Error in preprocess_csv: {str(e)}")  # Debug output
        raise Exception(f"Error preprocessing CSV: {str(e)}")

def identify_hook_effect(std_df, analyte, corrected_col):
    """
    Identifies the hook effect in standard curves by finding where MFI starts to decrease
    despite increasing concentration.
    
    Returns:
    - sorted_df: DataFrame with standards sorted by concentration
    - peak_idx: Index of the peak MFI (after which the hook effect occurs)
    - hook_present: Boolean indicating whether a hook effect was detected
    - hook_conc: Concentration at which the hook effect begins
    """
    # Get only rows with valid concentration and MFI values
    mask = (std_df['Concentration'] > 0) & np.isfinite(std_df['Concentration'])
    mask &= (std_df[corrected_col] > 0) & np.isfinite(std_df[corrected_col])
    
    # If not enough valid points, return immediately
    if mask.sum() < 3:
        return std_df.loc[mask].sort_values('Concentration'), -1, False, np.nan
    
    # Sort by concentration (ascending)
    sorted_df = std_df.loc[mask].sort_values('Concentration').reset_index(drop=True)
    
    # Find the peak MFI value
    peak_idx = sorted_df[corrected_col].idxmax()
    peak_mfi = sorted_df.loc[peak_idx, corrected_col]
    peak_conc = sorted_df.loc[peak_idx, 'Concentration']
    
    # Check if the peak is not at the highest concentration (possible hook effect)
    if peak_idx < len(sorted_df) - 1:
        # Check if there's a significant drop after the peak (e.g., >10% decrease)
        last_mfi = sorted_df.loc[len(sorted_df)-1, corrected_col]
        
        # If the drop is more than 5% of the peak, consider it a hook effect
        if (peak_mfi - last_mfi) / peak_mfi > 0.05:  # 5% threshold
            return sorted_df, peak_idx, True, peak_conc
    
    # No hook effect detected, use all standards
    return sorted_df, len(sorted_df) - 1, False, np.nan

def fit_standard_curve(df, analyte, corrected_col, exclude_hook=True):
    """
    Fit 4PL standard curve to the standards data for a given analyte.
    If exclude_hook is True, detects and excludes standards showing hook effect.
    """
    # Identify hook effect in standards (if present)
    sorted_std_df, peak_idx, hook_present, hook_conc = identify_hook_effect(df, analyte, corrected_col)
    
    if hook_present and exclude_hook:
        # Use only standards up to the peak (exclude those showing hook effect)
        fit_df = sorted_std_df.iloc[:peak_idx+1]
        print(f"Hook effect detected for {analyte}. Using {len(fit_df)} out of {len(sorted_std_df)} standards.")
    else:
        # Use all standards
        fit_df = sorted_std_df
    
    # Extract concentration and MFI values for curve fitting
    x = np.log10(fit_df['Concentration'].values)
    y = fit_df[corrected_col].values
    
    if len(x) < 4:  # We need at least 4 points for 4PL fitting
        raise ValueError(f"Not enough valid points ({len(x)}) to fit for {analyte}. Minimum required is 4.")
    
    # Initial guesses
    A0, D0 = np.max(y), np.min(y)
    B0 = 1.0
    C0 = np.median(x)
    p0 = [A0, B0, C0, D0]
    
    # Parameter bounds: A>D, B>0, C within x-range, D<max
    lower = [D0, 0.01, np.min(x), 0]
    upper = [A0 * 1.5, 10, np.max(x), D0 * 1.2]
    
    popt, _ = curve_fit(
        four_pl, x, y, p0=p0,
        bounds=(lower, upper),
        maxfev=20000
    )
    
    # Return the curve parameters, all standard points, and hook information
    return popt, sorted_std_df, hook_present, peak_idx, hook_conc

def calculate_standard_efficiency(std_df, analyte):
    """Calculate % efficiency for each standard point"""
    # Calculate % efficiency as (interpolated concentration / expected concentration) * 100
    std_df[f"{analyte}_Efficiency"] = (std_df[f"{analyte}_Conc"] / std_df['Concentration']) * 100
    
    # Group by expected concentration and calculate mean efficiency for each dilution level
    efficiency_by_dilution = std_df.groupby('Concentration')[f"{analyte}_Efficiency"].mean().reset_index()
    
    return efficiency_by_dilution

def load_cutoff_values(cutoff_file):
    """Load cutoff values from a CSV file with columns 'Analyte' and 'Cutoff'"""
    try:
        cutoffs = pd.read_csv(cutoff_file)
        if not {'Analyte', 'Cutoff'}.issubset(cutoffs.columns):
            raise ValueError("Cutoff file must have 'Analyte' and 'Cutoff' columns")
        
        # Convert to dictionary for faster lookup
        cutoff_dict = dict(zip(cutoffs['Analyte'], cutoffs['Cutoff']))
        return cutoff_dict
    except Exception as e:
        print(f"Error loading cutoff file: {str(e)}")
        return None

def analyze_data(df, analytes, standards_file=None, cutoff_file=None, bg_sample="Background0", 
                output_prefix="", dilution_factor=1.0, detect_hook=True):
    """
    Perform 4PL curve fitting and concentration calculations for selected analytes
    """
    results = {}
    all_efficiency_data = {}  # Store efficiency data for all analytes
    all_hook_data = {}  # Store hook effect data for all analytes
    needs_retest_flags = {}  # Store flags for samples needing retesting
    
    try:
        # Create output directory for graphs and results
        output_dir = f"{output_prefix}_graphs"
        os.makedirs(output_dir, exist_ok=True)
        
        # Create a subfolder for individual analyte files
        analytes_subfolder = os.path.join(output_dir, "analytes")
        os.makedirs(analytes_subfolder, exist_ok=True)
        
        # Load cutoff values if provided
        cutoff_dict = None
        if cutoff_file and os.path.exists(cutoff_file):
            print(f"Loading cutoff values from {cutoff_file}")
            cutoff_dict = load_cutoff_values(cutoff_file)
        
        # Load standards data if provided, otherwise extract from sample names
        if standards_file and os.path.exists(standards_file):
            print(f"Loading standards from {standards_file}")
            std_map = pd.read_csv(standards_file)
            if not {'Sample', 'Concentration'}.issubset(std_map.columns):
                raise ValueError('Standards CSV must have columns: Sample, Concentration')
            df = df.merge(std_map, on='Sample', how='left')
            standards_source = "file"
        else:
            print("No standards file provided. Attempting to infer concentrations from sample names.")
            # Extract standard samples from the dataframe
            std_samples = [col for col in df['Sample'] if 'Standard' in str(col)]
            
            # Create concentration mapping based on standard sample names
            std_conc = {}
            for sample in std_samples:
                try:
                    # Try to extract number from "StandardN" format
                    std_num = int(sample.replace("Standard", ""))
                    # Set concentration based on standard number (modify this logic as needed)
                    std_conc[sample] = 10 ** (5 - 0.5 * (std_num - 1))  # Example mapping
                except:
                    std_conc[sample] = np.nan
            
            # Create and apply standard mapping
            std_map = pd.DataFrame({
                'Sample': list(std_conc.keys()),
                'Concentration': list(std_conc.values())
            })
            df = df.merge(std_map, on='Sample', how='left')
            standards_source = "auto"
            
            # Print the inferred standard concentrations for debugging
            if not std_map.empty:
                print("Inferred standard concentrations:")
                for _, row in std_map.iterrows():
                    print(f"  {row['Sample']}: {row['Concentration']}")
            else:
                print("Warning: No standards found in the data.")
        
        # Process each analyte
        for analyte in analytes:
            print(f"Processing analyte: {analyte}")
            
            # Subtract background
            bg_data = df[df['Sample'] == bg_sample]
            if len(bg_data) == 0:
                print(f"Warning: Background sample '{bg_sample}' not found. Using raw values.")
                df[f"{analyte}_corr"] = df[analyte]
                bg_mean = 0
            else:
                bg_mean = bg_data[analyte].mean()
                df[f"{analyte}_corr"] = df[analyte] - bg_mean
                
            print(f"Background mean for {analyte}: {bg_mean}")
            
            # Get standards only
            std_df = df[df['Concentration'].notna() & (df['Concentration'] > 0)].copy()
            
            if len(std_df) < 5:
                print(f"Warning: Not enough standards for {analyte}. Skipping.")
                continue
                
            try:
                # Fit standard curve with hook effect detection
                popt, all_std_df, hook_present, peak_idx, hook_conc = fit_standard_curve(
                    std_df, analyte, f"{analyte}_corr", exclude_hook=detect_hook
                )
                
                # Store hook effect information
                all_hook_data[analyte] = {
                    'hook_detected': hook_present,
                    'hook_concentration': hook_conc,
                    'standards_excluded': len(all_std_df) - (peak_idx + 1) if hook_present else 0
                }
                
                A, B, C, D = popt
                print(f"{analyte} fit params: A={A:.2f}, B={B:.2f}, C={C:.2f}, D={D:.2f}")
                
                # Determine the proper maximum MFI/concentration to use based on hook detection
                if hook_present and detect_hook:
                    # If hook detected and exclusion enabled, use the peak as the max
                    max_valid_std_conc = all_std_df.loc[peak_idx, 'Concentration']
                    max_std_mfi = all_std_df.loc[peak_idx, f"{analyte}_corr"]
                else:
                    # Otherwise use the highest concentration standard
                    max_valid_std_conc = all_std_df['Concentration'].max()
                    max_std_mfi_samples = all_std_df[all_std_df['Concentration'] == max_valid_std_conc][f"{analyte}_corr"]
                    max_std_mfi = max_std_mfi_samples.mean() if not max_std_mfi_samples.empty else np.nan
                
                # Get the minimum concentration/MFI from standards
                min_std_conc = all_std_df['Concentration'].min()
                min_std_mfi_samples = all_std_df[all_std_df['Concentration'] == min_std_conc][f"{analyte}_corr"]
                min_std_mfi = min_std_mfi_samples.mean() if not min_std_mfi_samples.empty else np.nan
                
                # Calculate concentrations for all samples based on the standard curve
                logx = invert_log_four_pl(df[f"{analyte}_corr"].values, *popt)
                
                # Create a mask for samples that are below the minimum standard MFI or negative
                below_min_mfi_mask = (df[f"{analyte}_corr"] < min_std_mfi) | (df[f"{analyte}_corr"] < 0)
                
                # Apply the 4PL formula but set to NaN for values below min MFI
                df[f"{analyte}_Conc"] = np.where(
                    np.isfinite(logx) & ~below_min_mfi_mask,  # Only valid logx AND not below min MFI
                    10 ** logx,  # Apply calculation
                    np.nan  # Otherwise set to NaN
                )
                
                # Apply dilution factor to get final concentrations ONLY for non-standards and non-background samples
                df[f"{analyte}_Final_Conc"] = df[f"{analyte}_Conc"].copy()  # Start with raw concentration
                
                # Apply dilution factor only to non-standards and non-background samples
                mask = (~df['Sample'].str.contains('Standard', case=False, na=False)) & (df['Sample'] != bg_sample)
                df.loc[mask, f"{analyte}_Final_Conc"] = df.loc[mask, f"{analyte}_Conc"] * dilution_factor
                
                # Calculate standard efficiency for each standard point
                std_df[f"{analyte}_Conc"] = np.where(
                    np.isfinite(invert_log_four_pl(std_df[f"{analyte}_corr"].values, *popt)),
                    10 ** invert_log_four_pl(std_df[f"{analyte}_corr"].values, *popt),
                    np.nan
                )
                
                efficiency_data = calculate_standard_efficiency(std_df, analyte)
                
                # Store the efficiency data
                all_efficiency_data[analyte] = efficiency_data
                
                # Calculate average efficiency for the analyte across all standards
                avg_efficiency = std_df[f"{analyte}_Efficiency"].mean()
                
                # Store the average efficiency in the main dataframe
                df[f"{analyte}_Avg_Efficiency"] = avg_efficiency
                
                # FLAG 1: Check if raw MFI is below the lowest standard MFI or negative after background subtraction
                # First create the flag for all samples
                df[f"{analyte}_Below_Min_MFI"] = ((df[f"{analyte}_corr"] < min_std_mfi) | 
                                                 (df[f"{analyte}_corr"] < 0)).astype(int)
                
                # Then set it to NaN for standards (but NOT for background sample)
                is_std = df['Sample'].str.contains('Standard', case=False, na=False)
                df.loc[is_std, f"{analyte}_Below_Min_MFI"] = np.nan
                
                # Background samples should be marked as below min MFI (value = 1)
                is_bg = (df['Sample'] == bg_sample)
                df.loc[is_bg, f"{analyte}_Below_Min_MFI"] = 1
                
                # FLAG 3: Check if raw MFI is above the reliable range (needs dilution)
                # If hook effect detected, use the hook concentration as the threshold
                if hook_present and detect_hook:
                    # Use the peak MFI from hook detection as the threshold
                    df[f"{analyte}_Above_Max_MFI"] = (df[f"{analyte}_corr"] > max_std_mfi).astype(int)
                    # Add a separate flag for values in the "hook zone" (unreliable interpolation)
                    hook_zone_mask = (df[f"{analyte}_corr"] > max_std_mfi) | \
                                    ((df[f"{analyte}_Conc"] > hook_conc) & np.isfinite(df[f"{analyte}_Conc"]))
                    df[f"{analyte}_In_Hook_Zone"] = hook_zone_mask.astype(int)
                else:
                    # Traditional approach - just check if MFI exceeds the highest standard
                    df[f"{analyte}_Above_Max_MFI"] = (df[f"{analyte}_corr"] > max_std_mfi).astype(int)
                    # No hook zone in this case
                    df[f"{analyte}_In_Hook_Zone"] = 0
                
                # Set to NaN for standards and background
                df.loc[is_std | is_bg, f"{analyte}_Above_Max_MFI"] = np.nan
                df.loc[is_std | is_bg, f"{analyte}_In_Hook_Zone"] = np.nan
                
                # FLAG 2: Check if interpolated value is below the defined cutoff (if provided)
                if cutoff_dict and analyte in cutoff_dict:
                    cutoff_value = cutoff_dict[analyte]
                    # Set flag for all samples 
                    df[f"{analyte}_Below_Cutoff"] = (df[f"{analyte}_Final_Conc"] < cutoff_value).astype(int)
                    
                    # Set to NaN for standards only
                    df.loc[is_std, f"{analyte}_Below_Cutoff"] = np.nan
                    
                    # Background samples should always be below cutoff (set to 1)
                    df.loc[is_bg, f"{analyte}_Below_Cutoff"] = 1
                    
                    # Samples with MFI below min standard should automatically be below cutoff
                    below_min_mfi = (df[f"{analyte}_Below_Min_MFI"] == 1) & ~is_bg
                    df.loc[below_min_mfi, f"{analyte}_Below_Cutoff"] = 1
                else:
                    # If no cutoff provided, set to NA
                    df[f"{analyte}_Below_Cutoff"] = np.nan
                
                # Plot standard curve with both log10 MFI and concentration
                fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 6))
                
                # Left plot: Standard curve with concentration on x-axis
                # Plot standards differently based on hook effect
                if hook_present and detect_hook:
                    # Standards before the hook (included in fit)
                    included_df = all_std_df.iloc[:peak_idx+1]
                    ax1.scatter(included_df['Concentration'], included_df[f"{analyte}_corr"], 
                               color='blue', label='Standards (Used in Fit)')
                    
                    # Standards after the hook (excluded from fit)
                    excluded_df = all_std_df.iloc[peak_idx+1:]
                    if not excluded_df.empty:
                        ax1.scatter(excluded_df['Concentration'], excluded_df[f"{analyte}_corr"], 
                                   color='gray', marker='x', label='Standards (Excluded - Hook Effect)')
                else:
                    # All standards included
                    ax1.scatter(all_std_df['Concentration'], all_std_df[f"{analyte}_corr"], 
                               label='Standards')
                
                # Fitted curve (on finer x grid)
                if hook_present and detect_hook:
                    # Only plot curve up to the hook concentration
                    x_fit = np.logspace(min(np.log10(all_std_df['Concentration'])), 
                                      np.log10(hook_conc), 100)
                else:
                    # Plot across the full range
                    x_fit = np.logspace(min(np.log10(all_std_df['Concentration'])), 
                                      max(np.log10(all_std_df['Concentration'])), 100)
                
                y_fit = four_pl(np.log10(x_fit), *popt)
                ax1.plot(x_fit, y_fit, 'r-', label='4PL Fit')
                
                # Add horizontal line for minimum standard MFI
                ax1.axhline(y=min_std_mfi, color='orange', linestyle='--', alpha=0.7, 
                           label=f'Min Std MFI: {min_std_mfi:.1f}')
                
                # Add horizontal line for maximum reliable standard MFI
                ax1.axhline(y=max_std_mfi, color='red', linestyle='--', alpha=0.7, 
                           label=f'Max Reliable MFI: {max_std_mfi:.1f}')
                
                # If hook effect detected, add vertical line at hook concentration
                if hook_present and detect_hook:
                    ax1.axvline(x=hook_conc, color='purple', linestyle='--', alpha=0.7,
                               label=f'Hook Point: {hook_conc:.2f}')
                
                ax1.set_xscale('log')
                ax1.set_xlabel('Concentration')
                ax1.set_ylabel('Corrected MFI')
                
                # Add hook effect information to title if detected
                if hook_present and detect_hook:
                    ax1.set_title(f'Standard Curve for {analyte}\nHook Effect Detected and Corrected')
                else:
                    ax1.set_title(f'Standard Curve for {analyte}')
                
                ax1.legend()
                ax1.grid(True, which="both", ls="--", alpha=0.3)
                
                # Middle plot: Log10 MFI vs concentration
                if hook_present and detect_hook:
                    # Standards before the hook (included in fit)
                    ax2.scatter(included_df['Concentration'], np.log10(included_df[f"{analyte}_corr"]), 
                               color='blue', label='Standards (Used in Fit)')
                    
                    # Standards after the hook (excluded from fit)
                    if not excluded_df.empty:
                        ax2.scatter(excluded_df['Concentration'], np.log10(excluded_df[f"{analyte}_corr"]), 
                                   color='gray', marker='x', label='Standards (Excluded - Hook Effect)')
                else:
                    # All standards included
                    ax2.scatter(all_std_df['Concentration'], np.log10(all_std_df[f"{analyte}_corr"]),
                               label='Standards')
                
                # Fitted curve on log10 scale
                ax2.plot(x_fit, np.log10(y_fit), 'r-', label='4PL Fit')
                
                # If hook effect detected, add vertical line at hook concentration
                if hook_present and detect_hook:
                    ax2.axvline(x=hook_conc, color='purple', linestyle='--', alpha=0.7,
                               label=f'Hook Point: {hook_conc:.2f}')
                
                ax2.set_xscale('log')
                ax2.set_xlabel('Concentration')
                ax2.set_ylabel('Log10(Corrected MFI)')
                ax2.set_title(f'Log10 Standard Curve for {analyte}')
                ax2.legend()
                ax2.grid(True, which="both", ls="--", alpha=0.3)
                
                # Right plot: Standard Efficiency vs Concentration
                # Filter efficiency data if hook effect detected
                if hook_present and detect_hook:
                    # Only show efficiency for standards used in the fit
                    valid_concentrations = included_df['Concentration'].unique()
                    plot_efficiency_data = efficiency_data[efficiency_data['Concentration'].isin(valid_concentrations)]
                else:
                    plot_efficiency_data = efficiency_data
                
                ax3.scatter(plot_efficiency_data['Concentration'], plot_efficiency_data[f"{analyte}_Efficiency"],
                           marker='o', s=80, label='Efficiency')
                ax3.plot(plot_efficiency_data['Concentration'], plot_efficiency_data[f"{analyte}_Efficiency"],
                         'r-', label='Efficiency Trend')
                
                # Add horizontal line at 100% efficiency for reference
                ax3.axhline(y=100, color='k', linestyle='--', alpha=0.5, label='100% Efficiency')
                
                # Add average efficiency line
                avg_efficiency_valid = plot_efficiency_data[f"{analyte}_Efficiency"].mean()
                ax3.axhline(y=avg_efficiency_valid, color='b', linestyle='-.', alpha=0.5, 
                           label=f'Avg Efficiency: {avg_efficiency_valid:.1f}%')
                
                # If hook effect detected, add vertical line at hook concentration
                if hook_present and detect_hook:
                    ax3.axvline(x=hook_conc, color='purple', linestyle='--', alpha=0.7,
                               label=f'Hook Point: {hook_conc:.2f}')
                
                ax3.set_xscale('log')
                ax3.set_xlabel('Expected Concentration')
                ax3.set_ylabel('% Efficiency')
                ax3.set_title(f'Standard Efficiency for {analyte}')
                
                # Add efficiency values as text labels
                for idx, row in plot_efficiency_data.iterrows():
                    ax3.annotate(f"{row[f'{analyte}_Efficiency']:.1f}%", 
                                xy=(row['Concentration'], row[f'{analyte}_Efficiency']),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=9)
                
                ax3.grid(True, which="both", ls="--", alpha=0.3)
                ax3.legend()
                
                plt.tight_layout()
                
                # Save plot in the main directory (not in subfolder)
                graph_path = os.path.join(output_dir, f"{analyte}_standard_curve.png")
                plt.savefig(graph_path, dpi=300)
                plt.close()
                
                # Save efficiency data in the analytes subfolder
                efficiency_data.to_csv(os.path.join(analytes_subfolder, f"{analyte}_efficiency.csv"), index=False)
                
                # Store results for this analyte including flags
                # Include Location column if it exists in the dataframe
                if 'Location' in df.columns:
                    results[analyte] = df[['Sample', 'Location', analyte, f"{analyte}_corr", f"{analyte}_Conc", 
                                         f"{analyte}_Final_Conc", f"{analyte}_Avg_Efficiency",
                                         f"{analyte}_Below_Min_MFI", f"{analyte}_Below_Cutoff", 
                                         f"{analyte}_Above_Max_MFI", f"{analyte}_In_Hook_Zone"]].copy()
                else:
                    results[analyte] = df[['Sample', analyte, f"{analyte}_corr", f"{analyte}_Conc", 
                                         f"{analyte}_Final_Conc", f"{analyte}_Avg_Efficiency",
                                         f"{analyte}_Below_Min_MFI", f"{analyte}_Below_Cutoff",
                                         f"{analyte}_Above_Max_MFI", f"{analyte}_In_Hook_Zone"]].copy()
                
                # Save individual analyte results in the analytes subfolder
                results[analyte].to_csv(os.path.join(analytes_subfolder, f"{analyte}_results.csv"), index=False)
                
                # Track samples that need retesting at lower dilution
                if 'Location' in df.columns:
                    needs_retest_flags[analyte] = df[[
                        'Sample', 'Location', f"{analyte}_Above_Max_MFI", f"{analyte}_In_Hook_Zone"
                    ]].copy()
                else:
                    needs_retest_flags[analyte] = df[[
                        'Sample', f"{analyte}_Above_Max_MFI", f"{analyte}_In_Hook_Zone"
                    ]].copy()
                
            except Exception as e:
                print(f"Error fitting curve for {analyte}: {str(e)}")
        
        # Save hook effect information
        if all_hook_data:
            hook_df = pd.DataFrame([
                {
                    'Analyte': analyte,
                    'Hook_Detected': data['hook_detected'],
                    'Hook_Concentration': data['hook_concentration'],
                    'Standards_Excluded': data['standards_excluded']
                }
                for analyte, data in all_hook_data.items()
            ])
            hook_df.to_csv(os.path.join(output_dir, f"{output_prefix}_hook_effect_data.csv"), index=False)
        
        # Save combined results
        if results:
            # Create a combined dataframe that respects the unique Sample+Location combinations
            if 'Location' in df.columns:
                # Get unique Sample+Location combinations
                sample_location_df = df[['Sample', 'Location']].drop_duplicates()
                combined = sample_location_df.copy()
            else:
                # If no Location column, just use unique samples
                samples = df['Sample'].unique()
                combined = pd.DataFrame({'Sample': samples})
            
            # First, add a column to identify standards and background
            combined['Is_Standard'] = combined['Sample'].str.contains('Standard', case=False, na=False)
            combined['Is_Background'] = combined['Sample'] == bg_sample
            
            # Add basic analysis results for each analyte (concentrations and flags)
            for analyte in results:
                if 'Location' in df.columns:
                    # If we have Location, merge on both Sample and Location
                    analyte_cols = [col for col in results[analyte].columns if col.startswith(f"{analyte}_")]
                    temp = results[analyte][['Sample', 'Location'] + analyte_cols].copy()
                    combined = combined.merge(temp, on=['Sample', 'Location'], how='left')
                else:
                    # Otherwise merge just on Sample
                    analyte_cols = [f"{analyte}_Conc", f"{analyte}_Final_Conc", f"{analyte}_Avg_Efficiency",
                                  f"{analyte}_Below_Min_MFI", f"{analyte}_Below_Cutoff", 
                                  f"{analyte}_Above_Max_MFI", f"{analyte}_In_Hook_Zone"]
                    temp = results[analyte][['Sample'] + analyte_cols].copy()
                    combined = combined.merge(temp, on='Sample', how='left')
            
            # Create a normalized efficiency data table (long format)
            efficiency_long = pd.DataFrame()
            
            for analyte in all_efficiency_data:
                # Get efficiency data for this analyte
                efficiency_df = all_efficiency_data[analyte].copy()
                
                # Add analyte column
                efficiency_df['Analyte'] = analyte
                
                # Rename columns for standardization
                efficiency_df = efficiency_df.rename(columns={f"{analyte}_Efficiency": "Efficiency"})
                
                # Append to the long format dataframe
                efficiency_long = pd.concat([efficiency_long, efficiency_df], ignore_index=True)
            
            # Save the efficiency data in long format
            if not efficiency_long.empty:
                efficiency_long.to_csv(os.path.join(output_dir, f"{output_prefix}_efficiency_data.csv"), index=False)
            
            # Create a pivot table of the efficiency data for easier viewing
            if not efficiency_long.empty:
                efficiency_pivot = efficiency_long.pivot_table(
                    index="Concentration", 
                    columns="Analyte", 
                    values="Efficiency"
                ).reset_index()
                
                efficiency_pivot.to_csv(os.path.join(output_dir, f"{output_prefix}_efficiency_pivot.csv"), index=False)
            
            # Create and save a report of samples that need retesting
            samples_to_retest = {}
            has_location = 'Location' in df.columns
            
            for analyte, flag_df in needs_retest_flags.items():
                # Exclude standards and background samples
                is_std = flag_df['Sample'].str.contains('Standard', case=False, na=False)
                is_bg = flag_df['Sample'] == bg_sample
                valid_samples = ~(is_std | is_bg)
                
                # Find samples with Above_Max_MFI flag set to 1
                above_max_filter = valid_samples & (flag_df[f"{analyte}_Above_Max_MFI"] == 1)
                
                # Find samples in the hook zone (if applicable)
                in_hook_filter = valid_samples & False  # Initialize as all False
                if detect_hook:
                    in_hook_filter = valid_samples & (flag_df[f"{analyte}_In_Hook_Zone"] == 1)
                
                # Combined filter for any samples needing retesting
                needs_retest_filter = above_max_filter | in_hook_filter
                
                if has_location:
                    # Include Sample and Location for each flagged sample
                    retest_info = flag_df.loc[needs_retest_filter, ['Sample', 'Location']].copy()
                    if not retest_info.empty:
                        # Group samples by their locations
                        samples_to_retest[analyte] = retest_info
                else:
                    # Just include Sample names
                    retest_samples = flag_df.loc[needs_retest_filter, 'Sample'].tolist()
                    if retest_samples:
                        samples_to_retest[analyte] = retest_samples
            
            if samples_to_retest:
                if has_location:
                    # Create a set of unique Sample+Location combinations from all analytes
                    all_sample_locations = pd.concat([df for df in samples_to_retest.values()])
                    all_sample_locations = all_sample_locations.drop_duplicates()
                    
                    # Create a matrix with Sample and Location as identifiers
                    retest_df = pd.DataFrame(columns=['Sample', 'Location'] + analytes)
                    
                    # Add all unique Sample+Location combinations
                    retest_df[['Sample', 'Location']] = all_sample_locations[['Sample', 'Location']]
                    
                    # Fill in the matrix with "Yes" for samples needing retesting
                    for analyte, sample_df in samples_to_retest.items():
                        for _, row in sample_df.iterrows():
                            # Find the matching row in retest_df
                            mask = (retest_df['Sample'] == row['Sample']) & (retest_df['Location'] == row['Location'])
                            retest_df.loc[mask, analyte] = "Yes"
                else:
                    # Create a set of unique samples from all analytes
                    all_samples = set()
                    for samples in samples_to_retest.values():
                        all_samples.update(samples)
                    
                    # Create a matrix with Sample as index
                    retest_df = pd.DataFrame(index=sorted(all_samples), columns=analytes)
                    
                    # Fill in the matrix
                    for analyte, samples in samples_to_retest.items():
                        for sample in samples:
                            retest_df.loc[sample, analyte] = "Yes"
                    
                    # Reset index to convert Sample from index to column
                    retest_df = retest_df.reset_index().rename(columns={"index": "Sample"})
                
                # Save to CSV
                retest_file = os.path.join(output_dir, f"{output_prefix}_samples_to_retest.csv")
                retest_df.to_csv(retest_file, index=False)
                print(f"Saved samples needing retesting to: {retest_file}")
            
            # Save combined results in the output directory
            combined.to_csv(os.path.join(output_dir, f"{output_prefix}_combined_results.csv"), index=False)
            return combined, results, efficiency_long, all_hook_data, samples_to_retest
        else:
            return None, None, None, None
    
    except Exception as e:
        print(f"Error in analysis: {str(e)}")
        raise

class LuminexGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Luminex Plate Reader Analysis Tool")
        self.geometry("900x850")  # Increased default window size
        
        # Ensure application fully quits when window is closed
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        # Create main frame with padding
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Add a canvas with scrollbar for the main content
        canvas = tk.Canvas(main_frame, borderwidth=0)
        scrollable_frame = ttk.Frame(canvas)
        
        vsb = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)
        
        vsb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        
        # Create window for the scrollable frame
        canvas_frame = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        
        # Configure canvas to adjust the scrollregion when the size of scrollable_frame changes
        def configure_scroll_region(event):
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas.itemconfig(canvas_frame, width=canvas.winfo_width())
        
        scrollable_frame.bind("<Configure>", configure_scroll_region)
        canvas.bind('<Configure>', lambda e: canvas.itemconfig(canvas_frame, width=e.width))
        
        # Add mouse wheel scrolling
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        canvas.bind_all("<MouseWheel>", _on_mousewheel)  # Windows
        canvas.bind_all("<Button-4>", lambda e: canvas.yview_scroll(-1, "units"))  # Linux
        canvas.bind_all("<Button-5>", lambda e: canvas.yview_scroll(1, "units"))  # Linux
        
        # File selection section
        file_frame = ttk.LabelFrame(scrollable_frame, text="1. File Selection", padding="10")
        file_frame.pack(fill=tk.X, pady=5)
        
        self.file_path = tk.StringVar()
        ttk.Label(file_frame, text="Input CSV File:").grid(row=0, column=0, sticky=tk.W, pady=5)
        ttk.Entry(file_frame, textvariable=self.file_path, width=70).grid(row=0, column=1, sticky=tk.EW, padx=5, pady=5)
        ttk.Button(file_frame, text="Browse...", command=self.load_file).grid(row=0, column=2, padx=5, pady=5)
        
        # Configure the column weights to make the entry field expand
        file_frame.columnconfigure(1, weight=1)
        
        # File preprocessing section
        preprocess_frame = ttk.LabelFrame(scrollable_frame, text="2. File Preprocessing", padding="10")
        preprocess_frame.pack(fill=tk.X, pady=5)
        
        # Create a frame for plate format options
        plate_format_frame = ttk.Frame(preprocess_frame)
        plate_format_frame.grid(row=0, column=0, columnspan=3, sticky=tk.W, padx=5, pady=5)
        
        # Plate format selection
        ttk.Label(plate_format_frame, text="Plate Format:").grid(row=0, column=0, sticky=tk.W, padx=(0,5))
        self.plate_format_var = tk.StringVar(value="96-well")
        plate_format_combo = ttk.Combobox(plate_format_frame, textvariable=self.plate_format_var, 
                                         values=["96-well", "384-well"], state="readonly", width=15)
        plate_format_combo.grid(row=0, column=1, sticky=tk.W, padx=5)
        
        # Full plate checkbox and well count entry
        self.is_full_plate_var = tk.BooleanVar(value=True)
        full_plate_check = ttk.Checkbutton(plate_format_frame, text="Full Plate", 
                                         variable=self.is_full_plate_var, 
                                         command=self.toggle_well_count)
        full_plate_check.grid(row=0, column=2, sticky=tk.W, padx=20)
        
        # Well count input (initially disabled)
        ttk.Label(plate_format_frame, text="Number of Wells:").grid(row=0, column=3, sticky=tk.W, padx=5)
        self.num_wells_var = tk.StringVar(value="96")
        self.num_wells_entry = ttk.Entry(plate_format_frame, textvariable=self.num_wells_var, width=5, state="disabled")
        self.num_wells_entry.grid(row=0, column=4, sticky=tk.W, padx=5)
        
        # Update wells when plate format changes
        plate_format_combo.bind("<<ComboboxSelected>>", self.update_default_wells)
        
        # Preprocess button
        ttk.Button(preprocess_frame, text="Preprocess File", 
                  command=self.preprocess_file).grid(row=1, column=0, columnspan=2, padx=20, pady=5, sticky=tk.W)
        
        # Options section
        options_frame = ttk.LabelFrame(scrollable_frame, text="3. Analysis Options", padding="10")
        options_frame.pack(fill=tk.X, pady=5)
        
        # Configure column weights for better layout
        options_frame.columnconfigure(1, weight=1)
        options_frame.columnconfigure(2, weight=1)
        
        self.std_file_path = tk.StringVar()
        ttk.Label(options_frame, text="Standards CSV (optional):").grid(row=0, column=0, sticky=tk.W, pady=5)
        ttk.Entry(options_frame, textvariable=self.std_file_path, width=70).grid(row=0, column=1, columnspan=2, sticky=tk.EW, padx=5, pady=5)
        ttk.Button(options_frame, text="Browse...", command=self.load_std_file).grid(row=0, column=3, padx=5, pady=5)
        
        # Add info about standards file
        standards_info = ttk.Label(options_frame, text="(If not provided, concentrations will be auto-calculated from standard names)")
        standards_info.grid(row=1, column=1, columnspan=2, sticky=tk.W, pady=2)
        
        # Add cutoff file selection
        self.cutoff_file_path = tk.StringVar()
        ttk.Label(options_frame, text="Cutoff CSV (optional):").grid(row=2, column=0, sticky=tk.W, pady=5)
        ttk.Entry(options_frame, textvariable=self.cutoff_file_path, width=70).grid(row=2, column=1, columnspan=2, sticky=tk.EW, padx=5, pady=5)
        ttk.Button(options_frame, text="Browse...", command=self.load_cutoff_file).grid(row=2, column=3, padx=5, pady=5)
        
        self.bg_sample_var = tk.StringVar(value="Background0")
        ttk.Label(options_frame, text="Background Sample:").grid(row=3, column=0, sticky=tk.W, pady=5)
        ttk.Entry(options_frame, textvariable=self.bg_sample_var, width=20).grid(row=3, column=1, sticky=tk.W, padx=5, pady=5)
        
        self.dilution_factor_var = tk.StringVar(value="1")
        ttk.Label(options_frame, text="Dilution Multiplier:").grid(row=4, column=0, sticky=tk.W, pady=5)
        ttk.Entry(options_frame, textvariable=self.dilution_factor_var, width=10).grid(row=4, column=1, sticky=tk.W, padx=5, pady=5)
        
        # Add hook effect detection option
        self.detect_hook_var = tk.BooleanVar(value=True)
        ttk.Label(options_frame, text="Detect Hook Effect:").grid(row=5, column=0, sticky=tk.W, pady=5)
        ttk.Checkbutton(options_frame, variable=self.detect_hook_var).grid(row=5, column=1, sticky=tk.W, padx=5, pady=5)
        
        # Add the explanatory text on its own row to prevent it from being cut off
        explanation_label = ttk.Label(options_frame, text="(Enable to exclude standards showing hook effect)")
        explanation_label.grid(row=6, column=1, columnspan=2, sticky=tk.W, pady=2)
        
        self.output_prefix_var = tk.StringVar(value="LuminexResults")
        ttk.Label(options_frame, text="Output Prefix:").grid(row=7, column=0, sticky=tk.W, pady=5)
        ttk.Entry(options_frame, textvariable=self.output_prefix_var, width=20).grid(row=7, column=1, sticky=tk.W, padx=5, pady=5)
        
        # Analyte selection section
        analyte_frame = ttk.LabelFrame(scrollable_frame, text="4. Analyte Selection", padding="10")
        analyte_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(analyte_frame, text="Select Analyte:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.analyte_var = tk.StringVar()
        self.analyte_combo = ttk.Combobox(analyte_frame, textvariable=self.analyte_var, state="readonly", width=70)
        self.analyte_combo.grid(row=0, column=1, sticky=tk.EW, padx=5, pady=5)
        
        # Configure the column to expand
        analyte_frame.columnconfigure(1, weight=1)
        
        run_button = ttk.Button(analyte_frame, text="Run Analysis", command=self.run_analysis)
        run_button.grid(row=1, column=0, columnspan=2, pady=10)
        
        # Results section
        results_frame = ttk.LabelFrame(scrollable_frame, text="5. Results", padding="10")
        results_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Add scrollbar to results text area
        text_scroll = ttk.Scrollbar(results_frame)
        text_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.results_text = tk.Text(results_frame, wrap=tk.WORD, width=80, height=18, 
                                   yscrollcommand=text_scroll.set)
        self.results_text.pack(fill=tk.BOTH, expand=True)
        text_scroll.config(command=self.results_text.yview)
        
        # Status message - MOVED TO BOTTOM OF SCROLLABLE FRAME
        self.status_var = tk.StringVar()
        self.status_var.set("Ready. Please select a file to begin.")
        status_label = ttk.Label(scrollable_frame, textvariable=self.status_var, relief=tk.SUNKEN, padding=5)
        status_label.pack(fill=tk.X, pady=5, ipady=3)  # ipady adds a small vertical padding inside the label
        
        # Initialize variables
        self.file_loaded = False  # Flag to indicate if a file has been selected
        self.file_preprocessed = False  # Flag to indicate if a file has been preprocessed
        self.input_file_path = None  # Store the actual file path
        self.df_main = None
        self.unused_file = None
        self.timestamp = None
        self.analytes = []
        self.initial_output_dir = None
        self.plate_format = None
    
    def load_std_file(self):
        """Load standards CSV file"""
        path = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if path:
            self.std_file_path.set(path)
            messagebox.showinfo(
                "Standards File Selected", 
                "Standards file loaded. This will be used to determine the concentration values for standard samples.\n\n"
                "File must contain 'Sample' and 'Concentration' columns."
            )

    def load_cutoff_file(self):
        """Load cutoff CSV file"""
        path = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if path:
            self.cutoff_file_path.set(path)
    
    def load_file(self):
        """Open file dialog and load selected CSV file, but don't preprocess yet"""
        path = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not path:
            return
            
        # Just store the file path and update the UI
        self.file_path.set(path)
        self.input_file_path = path
        self.file_loaded = True
        self.file_preprocessed = False  # Reset preprocessing status
        
        # Clear previous results
        self.results_text.delete(1.0, tk.END)
        
        # Generate a timestamp for the output prefix
        current_time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.output_prefix_var.set(f"LuminexResults_{current_time}")
        
        self.status_var.set(f"File selected: {os.path.basename(path)}. Please choose plate format and preprocess.")
        self.results_text.insert(tk.END, f"File selected: {path}\n\n")
        self.results_text.insert(tk.END, "Please select the plate format and click 'Preprocess File' to continue.\n")
    
    def toggle_well_count(self):
        """Enable or disable the well count entry based on the full plate checkbox"""
        if self.is_full_plate_var.get():
            # Full plate selected, disable well count entry
            self.num_wells_entry.config(state="disabled")
            # Set the default well count based on plate format
            self.update_default_wells()
        else:
            # Partial plate selected, enable well count entry
            self.num_wells_entry.config(state="normal")
    
    def update_default_wells(self, event=None):
        """Update the default well count based on the selected plate format"""
        if self.plate_format_var.get() == "96-well":
            self.num_wells_var.set("96")
        else:
            self.num_wells_var.set("384")
    
    def preprocess_file(self):
        """Preprocess the loaded file with the selected plate format"""
        if not self.file_loaded or not self.input_file_path:
            messagebox.showwarning("Warning", "Please load a CSV file first.")
            return
            
        try:
            # Get the selected plate format
            plate_format = self.plate_format_var.get()
            
            # Get plate utilization settings
            is_full_plate = self.is_full_plate_var.get()
            
            # Get the number of wells if not a full plate
            num_wells = None
            if not is_full_plate:
                try:
                    num_wells = int(self.num_wells_var.get())
                    
                    # Validate the number of wells
                    max_wells = 96 if plate_format == "96-well" else 384
                    if num_wells <= 0 or num_wells > max_wells:
                        raise ValueError(f"Number of wells must be between 1 and {max_wells} for {plate_format} format")
                        
                except ValueError as e:
                    messagebox.showerror("Error", f"Invalid number of wells: {str(e)}")
                    return
            
            # Generate output prefix (or use existing one)
            output_prefix = self.output_prefix_var.get()
            
            # Create the output directory
            output_dir = f"{output_prefix}_graphs"
            os.makedirs(output_dir, exist_ok=True)
            
            # Now preprocess the file
            self.status_var.set(f"Preprocessing file (format: {plate_format}, {'full' if is_full_plate else 'partial'}) and saving metadata...")
            self.update_idletasks()
            
            df, unused_fname, ts, initial_dir, plate_format = preprocess_csv(
                self.input_file_path, output_prefix, plate_format, is_full_plate, num_wells
            )
            
            self.df_main = df
            self.unused_file = unused_fname
            self.timestamp = ts
            self.initial_output_dir = initial_dir
            self.plate_format = plate_format
            self.file_preprocessed = True
            
            # Update results text
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, f"File preprocessed successfully.\n")
            self.results_text.insert(tk.END, f"Plate format: {plate_format}\n")
            
            # Add information about plate utilization
            if is_full_plate:
                self.results_text.insert(tk.END, f"Full plate: {'96' if plate_format == '96-well' else '384'} wells\n")
            else:
                self.results_text.insert(tk.END, f"Partial plate: {num_wells} wells\n")
                
            self.results_text.insert(tk.END, f"Metadata saved to: {unused_fname}\n\n")
            
            # Identify analytes: columns between 'Sample' and 'Total Events'
            cols = list(df.columns)
            try:
                i1 = cols.index('Sample')
                i2 = cols.index('Total Events')
                self.analytes = cols[i1+1:i2]
                
                if not self.analytes:
                    raise ValueError("No analytes found between 'Sample' and 'Total Events' columns.")
                    
                # Update analyte dropdown
                self.analyte_combo['values'] = ["All"] + self.analytes
                self.analyte_combo.current(0)
                
                self.results_text.insert(tk.END, f"Found {len(self.analytes)} analytes:\n")
                self.results_text.insert(tk.END, ", ".join(self.analytes) + "\n\n")
                self.results_text.insert(tk.END, "Ready for analysis. Select an analyte and click 'Run Analysis'.\n")
                
                self.status_var.set(f"File preprocessed successfully. Ready for analysis.")
                
            except ValueError as e:
                self.status_var.set(f"Error: {str(e)}")
                messagebox.showerror("Error", f"CSV must contain 'Sample' and 'Total Events' columns.")
                
        except Exception as e:
            self.status_var.set(f"Error preprocessing file: {str(e)}")
            messagebox.showerror("Error", f"Failed to preprocess CSV: {str(e)}")
    
    def run_analysis(self):
        """Run analysis on loaded data with selected analytes"""
        if not self.file_preprocessed or self.df_main is None:
            messagebox.showwarning("Warning", "Please load and preprocess a CSV file first.")
            return
            
        choice = self.analyte_var.get()
        if not choice:
            messagebox.showwarning("Warning", "Please select an analyte.")
            return
            
        # Determine which analytes to process
        to_run = self.analytes if choice == "All" else [choice]
        
        self.status_var.set(f"Running analysis on {len(to_run)} analyte(s)...")
        self.update_idletasks()
        
        try:
            # Get output prefix
            output_prefix = self.output_prefix_var.get()
            
            # Get background sample name
            bg_sample = self.bg_sample_var.get()
            
            # Get hook effect detection preference
            detect_hook = self.detect_hook_var.get()
            
            # Get dilution factor with error checking
            try:
                dilution_factor = float(self.dilution_factor_var.get())
                if dilution_factor <= 0:
                    raise ValueError("Dilution factor must be positive.")
            except ValueError as e:
                messagebox.showerror("Error", f"Invalid dilution factor: {str(e)}")
                return
            
            # Get standards file if specified
            std_file = self.std_file_path.get() if self.std_file_path.get() else None
            
            # Get cutoff file if specified
            cutoff_file = self.cutoff_file_path.get() if self.cutoff_file_path.get() else None
            
            # Run analysis
            results = analyze_data(
                df=self.df_main,
                analytes=to_run,
                standards_file=std_file,
                cutoff_file=cutoff_file,
                bg_sample=bg_sample,
                output_prefix=output_prefix,
                dilution_factor=dilution_factor,
                detect_hook=detect_hook
            )
            
            combined_results, detailed_results, efficiency_data, hook_data, samples_to_retest = results
            
            # Move the unused data file to the new output directory if needed
            if hasattr(self, 'unused_file') and hasattr(self, 'initial_output_dir'):
                try:
                    # Define the new output directory
                    new_output_dir = f"{output_prefix}_graphs"
                    
                    # If the initial and new directories are different
                    if self.initial_output_dir != new_output_dir:
                        # Get the filename of the unused data file
                        unused_filename = os.path.basename(self.unused_file)
                        
                        # Path to the new location
                        new_unused_file = os.path.join(new_output_dir, unused_filename)
                        
                        # Move the file to the new directory
                        import shutil
                        shutil.copy2(self.unused_file, new_unused_file)
                        
                        # Delete the original file
                        os.remove(self.unused_file)
                        
                        # Check if the initial directory is empty
                        if len(os.listdir(self.initial_output_dir)) == 0:
                            # Remove the empty directory
                            os.rmdir(self.initial_output_dir)
                            print(f"Removed empty directory: {self.initial_output_dir}")
                        
                        # Update the unused file path
                        self.unused_file = new_unused_file
                except Exception as e:
                    print(f"Warning: Failed to move unused data file: {str(e)}")
            
            if combined_results is not None:
                self.results_text.insert(tk.END, "\nAnalysis completed successfully.\n")
                output_dir = f"{output_prefix}_graphs"
                analytes_subfolder = os.path.join(output_dir, "analytes")
                self.results_text.insert(tk.END, f"All results and graphs saved in: {output_dir}\n")
                self.results_text.insert(tk.END, f"Individual analyte CSV files saved in: {analytes_subfolder}\n")
                self.results_text.insert(tk.END, f"Standard curve plots (PNG files) saved in main directory for easy viewing\n\n")
                
                # Show plate format, dilution factor, and hook detection status
                self.results_text.insert(tk.END, f"Plate format: {self.plate_format}\n")
                self.results_text.insert(tk.END, f"Dilution factor applied: {dilution_factor} (only to non-standard, non-background samples)\n")
                self.results_text.insert(tk.END, f"Hook effect detection: {'Enabled' if detect_hook else 'Disabled'}\n")
                
                # Add information about standard concentrations source
                if 'standards_source' in locals():
                    if standards_source == "auto":
                        self.results_text.insert(tk.END, "\nNOTE: Standard concentrations were automatically calculated from sample names.\n")
                        self.results_text.insert(tk.END, "      This assumes standards follow a dilution series where:\n")
                        self.results_text.insert(tk.END, "      Standard1 = 100000, Standard2 = 31623, Standard3 = 10000, etc.\n")
                        self.results_text.insert(tk.END, "      For more accurate results, provide a standards.csv file with actual concentrations.\n\n")
                    else:
                        self.results_text.insert(tk.END, f"\nStandard concentrations loaded from file: {std_file}\n\n")
                
                # Show hook effect summary if detected and enabled
                if hook_data and detect_hook:
                    hook_detected = any(data['hook_detected'] for data in hook_data.values())
                    if hook_detected:
                        self.results_text.insert(tk.END, "Hook Effect Summary:\n")
                        for analyte, data in hook_data.items():
                            if data['hook_detected']:
                                self.results_text.insert(tk.END, f"  - {analyte}: Hook effect detected at concentration {data['hook_concentration']:.2f}. " +
                                                       f"{data['standards_excluded']} standards excluded from fit.\n")
                        self.results_text.insert(tk.END, "\n")
                
                # Show sample of results
                num_samples = min(5, len(combined_results))
                self.results_text.insert(tk.END, f"Sample of results (first {num_samples} rows):\n")
                
                # Format and display sample results
                result_sample = combined_results.head(num_samples).to_string()
                self.results_text.insert(tk.END, result_sample + "\n\n")
                
                # Note about the flags
                self.results_text.insert(tk.END, "Added flags to the output:\n")
                self.results_text.insert(tk.END, "1. [Analyte]_Below_Min_MFI: Value is 1 if raw MFI is below the lowest standard or negative after background subtraction\n")
                self.results_text.insert(tk.END, "   - Standards: Flag set to NaN (not applicable)\n") 
                self.results_text.insert(tk.END, "   - Background samples: Flag set to 1 (always considered negative)\n")
                self.results_text.insert(tk.END, "2. [Analyte]_Below_Cutoff: Value is 1 if interpolated concentration is below the defined cutoff\n")
                self.results_text.insert(tk.END, "   - Standards: Flag set to NaN (not applicable)\n")
                self.results_text.insert(tk.END, "   - Background samples: Flag set to 1 (always considered negative)\n")
                self.results_text.insert(tk.END, "   - Samples below Min MFI: Automatically flagged as below cutoff\n")
                self.results_text.insert(tk.END, "3. [Analyte]_Above_Max_MFI: Value is 1 if raw MFI is above the highest reliable standard's MFI\n")
                self.results_text.insert(tk.END, "   - Standards and background: Flag set to NaN (not applicable)\n")
                self.results_text.insert(tk.END, "   - Samples above Max MFI: Should be rerun with more dilution\n")
                
                # New hook zone flag
                if detect_hook:
                    self.results_text.insert(tk.END, "4. [Analyte]_In_Hook_Zone: Value is 1 if the sample might be affected by the hook effect\n")
                    self.results_text.insert(tk.END, "   - Standards and background: Flag set to NaN (not applicable)\n")
                    self.results_text.insert(tk.END, "   - Samples in hook zone: May have unreliable interpolated values\n\n")
                else:
                    self.results_text.insert(tk.END, "\n")
                
                # Note about standard efficiency
                self.results_text.insert(tk.END, "Standard efficiency data has been calculated and saved in long format.\n")
                self.results_text.insert(tk.END, "Efficiency data is now available in two formats:\n")
                self.results_text.insert(tk.END, "1. Long format: Analyte, Concentration, Efficiency columns\n")
                self.results_text.insert(tk.END, "2. Pivot table: Concentration rows, Analyte columns\n")
                self.results_text.insert(tk.END, "Efficiency represents how well the fitted curve predicts the expected standard values.\n\n")
                
                # Display samples that need retesting
                if samples_to_retest:
                    retest_file = os.path.join(output_dir, f"{output_prefix}_samples_to_retest.csv")
                    self.results_text.insert(tk.END, "Samples Requiring Retesting at Lower Dilution:\n")
                    
                    # Calculate total samples (handling both DataFrame and list formats)
                    has_location = 'Location' in self.df_main.columns
                    if has_location:
                        total_samples = sum(len(df) for df in samples_to_retest.values())
                    else:
                        total_samples = sum(len(samples) for samples in samples_to_retest.values())
                        
                    self.results_text.insert(tk.END, f"Found {total_samples} sample{'s' if total_samples != 1 else ''} that need{'s' if total_samples == 1 else ''} retesting for one or more analytes.\n")
                    self.results_text.insert(tk.END, f"A detailed report has been saved to: {retest_file}\n\n")
                    
                    # Display a few examples of samples that need retesting
                    max_display = 3  # Maximum number of analytes to display
                    max_samples = 5  # Maximum number of samples to display per analyte
                    
                    analytes_shown = 0
                    for analyte, samples in samples_to_retest.items():
                        if analytes_shown >= max_display:
                            break
                            
                        if has_location:
                            # For DataFrame format with Location
                            num_samples = len(samples)
                            self.results_text.insert(tk.END, f"  {analyte}: {num_samples} sample{'s' if num_samples != 1 else ''} need{'s' if num_samples == 1 else ''} retesting\n")
                            
                            # Show up to max_samples sample names with locations
                            sample_display = samples.head(max_samples)
                            if not sample_display.empty:
                                for _, row in sample_display.iterrows():
                                    self.results_text.insert(tk.END, f"    Sample: {row['Sample']}, Location: {row['Location']}\n")
                                
                                if num_samples > max_samples:
                                    self.results_text.insert(tk.END, f"    ...and {num_samples - max_samples} more...\n")
                        else:
                            # For list format without Location
                            num_samples = len(samples)
                            self.results_text.insert(tk.END, f"  {analyte}: {num_samples} sample{'s' if num_samples != 1 else ''} need{'s' if num_samples == 1 else ''} retesting\n")
                            
                            # Show up to max_samples sample names
                            sample_display = samples[:max_samples]
                            if sample_display:
                                self.results_text.insert(tk.END, f"    Examples: {', '.join(sample_display)}")
                                if len(samples) > max_samples:
                                    self.results_text.insert(tk.END, f" and {len(samples) - max_samples} more...")
                                self.results_text.insert(tk.END, "\n")
                            
                        analytes_shown += 1
                    
                    if analytes_shown < len(samples_to_retest):
                        self.results_text.insert(tk.END, f"  ...and {len(samples_to_retest) - analytes_shown} more analytes with samples needing retesting.\n")
                    
                    self.results_text.insert(tk.END, "\n")
                
                self.status_var.set("Analysis complete.")
                messagebox.showinfo("Success", f"Analysis complete. Results saved to {output_dir}/")
            else:
                self.results_text.insert(tk.END, "\nAnalysis failed to produce results.\n")
                self.status_var.set("Analysis failed.")
                messagebox.showerror("Error", "Analysis failed to produce results.")
                
        except Exception as e:
            self.results_text.insert(tk.END, f"\nError during analysis: {str(e)}\n")
            self.status_var.set("Analysis error.")
            messagebox.showerror("Error", f"Analysis failed: {str(e)}")

    def on_closing(self):
        """Handle window close event to ensure application fully quits"""
        try:
            # Display confirmation dialog
            if messagebox.askokcancel("Quit", "Do you want to quit the application?"):
                # Clean up mouse wheel binding
                self.unbind_all("<MouseWheel>")
                self.unbind_all("<Button-4>")
                self.unbind_all("<Button-5>")
                
                # Destroy the tkinter root window
                self.destroy()
                # Force application exit
                import sys
                sys.exit(0)
        except Exception as e:
            print(f"Error during application shutdown: {str(e)}")
            # Ensure application exits even if there's an error
            import os
            os._exit(0)

if __name__ == '__main__':
    app = LuminexGUI()
    app.mainloop()