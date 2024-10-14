# HiChIP Manuscript Code

This repository contains the analysis and plotting scripts for the HiChIP manuscript, designed to process and analyze data from TCGA enhancers and gene interactions.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Code Overview](#code-overview)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/HiChIP-Manuscript.git
cd HiChIP-Manuscript
```

Install the necessary dependencies for R and Python scripts. You may need `R`, `Python 3.x`, and the following libraries:

- **R packages**: `ggplot2`, `dplyr`, etc.
- **Python packages**: `matplotlib`, `numpy`, etc.

For R:

```R
install.packages(c("ggplot2", "dplyr", ...))
```

For Python:

```bash
pip install matplotlib numpy
```

## Usage

The repository contains multiple R and Python scripts, each performing different tasks related to HiChIP analysis. Below are the primary scripts and their usage:

1. **TCGA Enhancer Analysis**:
   - `1_TCGA_Enhancers_Gene_Interactions.r`: Analysis of enhancer-gene interactions using TCGA data.
   - `2_TCGA_Enhancer_Copy_Number_Regression.r`: Performs regression analysis on enhancer copy number data.
   
   Example:
   ```bash
   Rscript 1_TCGA_Enhancers_Gene_Interactions.r
   ```

2. **RNA Modeling**:
   - `3_TCGA_RNA_Modeling.r`: RNA modeling analysis for TCGA data.
   - `4_TCGA_RNA_Modeling_Analysis.r`: Detailed RNA modeling analysis.
   
   Example:
   ```bash
   Rscript 3_TCGA_RNA_Modeling.r
   ```

3. **Visualization**:
   - `Plot_heatmaps.R`: Plots heatmaps for various interaction data.
   - `plot_neoloop_dist.py`: Python script for plotting neo-loop distributions.
   
   Example:
   ```bash
   python plot_neoloop_dist.py
   ```

For more detailed usage, see each script's comments and documentation.

## Code Overview

- **Enhancer Analysis**:
   - `1_TCGA_Enhancers_Gene_Interactions.r`
   - `2_TCGA_Enhancer_Copy_Number_Regression.r`
   - `TCGA_King_Enhancer_Interaction_Functions.r`
  
- **RNA Modeling**:
   - `3_TCGA_RNA_Modeling.r`
   - `4_TCGA_RNA_Modeling_Analysis.r`
   - `TCGA_King_RNA_modeling_functions.r`
  
- **Visualization**:
   - `Plot_heatmaps.R`
   - `plot_neoloop.py`
   - `Plot_EIS.R`

- **Workflow Scripts**:
   - `HiChIP_decomposition_workflow.R`: Main workflow for HiChIP data decomposition.
   - `Non_coding_workflow_Final.R`: Non-coding region workflow.

## Contributing

If you'd like to contribute, please fork the repository and use a feature branch. Pull requests are welcome.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

This project was developed with data from TCGA and contributions from multiple collaborators in the HiChIP community.

