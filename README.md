# SNV2
Reconstruction of SNV-SNV Interaction Networks via Computational
Detection of Variant Co-expression in RNA-sequencing Data

This toolkit contains the required scripts to transform sequencing files into SNV^2 input files then run the MatrixEQTL R package to determine significant variation-expression relationships.

## Getting Started:
These instructions will get you a copy of the scripts up and running on your machine for development and testing purposes. See Running the scripts for notes on how to use the project on a live system.

### Prerequisites:
* An R (>3.5) installation with the following packages: corrplot, data.table, dplyr, factoextra, ggplot2, ggpubr, gtools, MatrixEQTL, pacman, pheatmap, psych, RColorBrewer, stringr, tidyr
* Each of the following scripts copied to a working directory on your machine:
  
  ```
  A-build_snvsnv_matrix_CLine.R
  B-build_snvsnv_cov_matrix_CLine.R
  C-run_snvsnv_matrixEQTL_cov_CLine.R
  D-scatter_plots_SNV2_CLine.R
  E-plot_snvsnv_eQTL_correlograms_CLine.R
  ```
You can obtain the full toolkit [here.](https://github.com/pdbousounis/SNV2)

*Output *.csv* files must be obtained from our ReadCounts tool (https://github.com/HorvathLab/NGS/tree/master/readCounts) containing the read counts extracted per SNV for each sample


## Running the scripts


### A-build\_snvsnv_matrix_CLine.R
Transforms the read counts into a variant fraction matrix with information from all provided samples

#### Input arguments:
* A directory containing the .csv files from the output of Readcounts
* A user-defined prefix for the output files

#### Output: 
* *All output files are saved in a new directory ((prefix)_SNV2_Output) created within the .csv file input directory (it is created if not already present)*
* A SNV matrix of VAFRNA values for each SNV by sample
* A matrix of SNV locations

#### Sample command:
```
Rscript A-build_snvsnv_matrix_CLine.R /home/readcounts/ nerve_234
```
&nbsp;

***

&nbsp;

### B-build\_snvsnv_cov_matrix_CLine.R (\*see *NOTE* in next section)
Creates a covariant matrix from sample metadata (GTEx used here) including Age, Gender, Ethnicity, and the top 10 prinicpal components

#### Input arguments:
* (prefix)_SNV2_Output directory from previous step
* A directory containing sample metadata files (here construced for GTEx data)
* Desired output file prefix

#### Output:
In the (prefix)_SNV2_Output directory from the previous step:
* A matrix of covariates and top 10 PCs for each sample
* Scree and PC contribution plots (.png)
  
#### Sample command:
```
Rscript nerve_234_SNV2_Output sample_metadata_directory nerve_234
```
&nbsp;

***

&nbsp;  

### C-run\_snvsnv_matrixEQTL_cov_CLine.R

#### Input arguments:
* <prefix>_SNV2_Output directory from previous step
* Desired output file prefix 

* *NOTE*: Covariates file is optional (requires modification of the script as specified in the in-code documentation, line 32)

#### Output:
In the (prefix)_SNV2_Output directory from the previous step:
* One file with the cis eQTLs with a p-value < 0.00001
* One file with the trans eQTLs with a p-value < 0.00001
* QQ-plot of distant vs local p-value distributions (.png)

#### Sample command:
```
Rscript C-run_snvsnv_matrixEQTL_cov_CLine.R nerve_234_SNV2_Output nerve_234
```
&nbsp;

***

&nbsp;

### D-scatter\_plots_SNV2_CLine.R
Visualize individual SNV-SNV correlations across all patients as scatter plots

#### Input arguments:
* (prefix)_SNV2_Output directory from previous step
* Designate data to plot, one of 'cis' or 'trans'
* Desired prefix for output file

#### Output:
In the (prefix)_SNV2_Output directory from the previous step:
* A single .png file of the top 200 significant (FDR < 0.05) SNV-SNV correlations, sorted by lowest FDR value

#### Sample command:
```
Rscript D-scatter_plots_SNV2_CLine.R nerve_234_SNV2_Output cis nerve_234
```
&nbsp;

***

&nbsp;

### E-plot\_snvsnv_eQTL_correlograms_CLine.R
Create SNV-SNV correlation plots for a designted chromosome
* *NOTE*: Requires a SeattleSeq annotation file generated from the SNV location matrix and located within the (prefix)_SNV2_Output directory

#### Input arguments:
* (prefix)_SNV2_Output directory from previous step
* Chromosome to plot 
* Desired prefix for output file

#### Output:
In the (prefix)_SNV2_Output directory from the previous step:
* A correlogram of VAFRNA values for each significant SNV-SNV interaction  

#### Sample command:
```
Rscript nerve_234_SNV2_Output 21 nerve_234
```
&nbsp;

## Authors:
* **Pavlos Bousounis** (pbousounis@protonmail.com)

## Acknowledgements:
- Liam F Spurr, Justin Sein, and Anelia Horvath for support and assistance in the development of this toolkit
- The MatrixEQTL team for their sample code and R package upon which *C-run\_snvsnv_matrixEQTL_cov_CLine.R* is based
