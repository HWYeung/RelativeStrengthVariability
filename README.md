# Relative Strength Variability

Code for the paper: *Relative strength variability measures for brain structural connectomes and their relationship with cognitive functioning*: [https://onlinelibrary.wiley.com/doi/10.1002/hbm.70314](https://onlinelibrary.wiley.com/doi/10.1002/hbm.70314)

## Scripts Overview

- `NodeRelativeStrengthVariability.m` – Computes RSV and hierarchical RSV measures.  
- `GSRegressResid.m` – Extracts incremental contributions in predicting the g-factor in progressive models.  
- `LRAMissing.m` – Performs low-rank approximation to extract g-factor for incomplete data.  
- `Hierarchical3D_v2.m` – Generates the hierarchical complexity measure from *The complex hierarchical topology of EEG functional connectivity*.  
- `local_assortativity.m` – Computes local assortativity.

## Required Data

All required files have been uploaded to the UK Biobank to facilitate replication.

## Data Preparation

1. Run `PreparingForAnalysis.m` to prepare the tabular data.  
2. Run `compute_all_graph_measures.m` to extract all connectomes and compute graph measures. Outputs are stored as 3D arrays.

## Association Analysis

1. Run `graph_measures_association_analysis.m` to assess associations of graph metrics with g-factor, controlling for covariates.

## Plots

1. Run `graph_measure_correlation_plotting.m` to create the plot for the correlations among graph metrics
2. Run `Corrplots_varwithbrainmeasures.m` to create the plot for the correlations between the RSV measures and other measures (i.e. age, sex, mean edge weight, global volumetric measures)
3. Run `plot_incremental_variance.m` to create plots for incremental R-squared contribution to predicting g-factor

## Notes

- Ensure all scripts are run in the correct directory and required MATLAB toolboxes are installed.  
- Outputs of each step (connectomes, graph measures, RSV metrics) are required for subsequent analysis steps.
