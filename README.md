# Urban mobilization diffusion

This repository contains the data, scripts, and documentation for analyzing urban social mobilization and diffusion dynamics. It focuses on exploring network structures, coalition formation, and how diffusion processes drive contentious actions, as exemplified by the 2011 Peñalolén urban mobilization case.

## Repository Structure
```
urban-diffusion-dynamics/
│
├── README.md               # Project overview, instructions, and repository structure
├── LICENSE                 # License of the project (e.g., MIT License)
│
├── data/
│   ├── raw/                # Original datasets (e.g., CSV files with network and attribute data)
│   └── processed/          # Cleaned or transformed data used in the analysis
│
├── scripts/
│   ├── data_preparation.R  # Script for reading and preprocessing datasets
│   ├── network_analysis.R  # Script for descriptive network analysis and visualization
│   ├── diffusion_analysis.R# Script dedicated to analyzing diffusion dynamics
│   └── ergm_model.R        # Script for constructing and evaluating ERGM models
│
├── notebooks/
│   └── exploratory_analysis.Rmd  # R Markdown notebook for exploratory and reproducible analysis
│
└── docs/
    ├── repository_structure.md   # Additional documentation of project structure and methodology
    └── overleaf_project_link.txt   # File containing the link and connection details to the Overleaf document
```

## Overleaf Documentation

For a detailed LaTeX document outlining the project background, methodology, and analysis, please refer to the Overleaf project. You can access the document using the following link:

[Access Overleaf Document](https://www.overleaf.com/project/6794012e99f7f84f2a7235b8)

If you encounter any issues accessing the Overleaf document or require further information, please contact the project maintainers.

## Project Overview

- **Background:**  
  This study investigates how multiplex networks and diffusion dynamics shape urban social movements. It examines the role of various tie types (trust, resource, kinship, etc.) and network configurations in facilitating coalition building and triggering the diffusion of social mobilization.

- **Objectives:**
  - **Network Analysis:** Evaluate the network structure of urban social movements.
  - **ERGM Modeling:** Use Exponential Random Graph Models (ERGM) to assess the impact of different tie types on mobilization events.
  - **Diffusion Dynamics:** Analyze how information, influence, or innovation diffuses within urban communities.

## Getting Started

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/yourusername/urban-diffusion-dynamics.git
   cd urban-diffusion-dynamics
   ```

2. Install Dependencies:
Ensure you have R installed.
Install the required R packages. In R, run:

```r
install.packages(c("tidyverse", "igraph", "statnet", "ergm", "rmarkdown"))
```
3. Run the Scripts:
Execute data_preparation.R to load and preprocess the datasets.
Run network_analysis.R to perform descriptive network analysis and generate visualizations.
Run diffusion_analysis.R for a detailed examination of diffusion dynamics.
Finally, run ergm_model.R to fit and evaluate the ERGM models.

4. Explore the Notebook:
Open notebooks/exploratory_analysis.Rmd in RStudio to review the exploratory analysis in a reproducible format.
Contributing
Contributions are welcome! Fork the repository, make your changes, and submit a pull request. Your input will help improve the analysis and documentation.

## License
This project is licensed under the MIT License.

## Contact
For any questions or feedback, please contact [Roberto Cantillan] at [ricantillan@uc.cl].
