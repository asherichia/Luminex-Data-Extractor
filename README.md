
# Luminex Data Extractor

A Python script for processing and analyzing Luminex assay data (currently only tested on Luminex Intelliflex files). This tool enables users to import Luminex CSV output files (e.g. direct exports from Luminex machines), apply a 4-parameter logistic (4PL) standard curve fit, interpolate unknown sample concentrations, and identify additional samples for further dilution/testing.

## Features

- GUI interface for selecting input files and analytes
- Supports standard 4PL curve fitting for each analyte
- Calculates interpolated concentrations for unknowns
- Flags samples that fall outside of fit bounds or show poor fitting
- Generates plots of standard curves with sample points
- Outputs clean, analysis-ready CSV results

## Requirements

- Python 3.8+
- `pandas`
- `numpy`
- `matplotlib`
- `tkinter`
- `scipy`

You can install the required packages using:

```bash
pip install pandas numpy matplotlib scipy
```

> Note: `tkinter` is included in most Python installations by default. On macOS, use `python3 -m tkinter` to test if it's available.

## Standards File Format

The standards CSV file should contain known concentrations for each analyte, formatted as follows:

| Analyte    | Concentration | 
|------------|---------------|
| Standard1  | 1000          |
| Standard2  | 500           |
| Standard3  | 250           | 

- Each row is a standard concentration for the analyte(s)
- In this case, the Standard1 through to Standard8 is used across all analytes


## Output Files

After analysis, the script generates:

### 1. `results.csv`
A summary CSV with interpolated concentrations and flags:

<To insert>



### 2. PNG plots

Each analyte generates a curve-fit plot saved as `StandardCurve_<Analyte>.png` showing:
- Standard curve with fitted 4PL curve
- Standard curve with fitted 4PL curve (Log10)
- Efficiency plot for the standard curves

## Usage

Run the script via command line:

```bash
python Combined_Analysis.py
```

This will launch a GUI where you can:
- Select your Luminex CSV data file
- Choose the standards file
- Select all analytes or specific analytes
- Generate output plots and CSV files, whilst tidying up files

## License

[MIT License](LICENSE)

---

## Author

Ashley Otter  
[GitHub Profile](https://github.com/asherichia)
