# QUANTIFA: Quantile Analysis of Temperature using Individual Foraminiferal Analyses
QUANTIFA is a user-friendly exploratory and data analysis tool for analyzing individual foraminiferal data, constraining the uncertainties inherent to IFA statisics, and testing the resolving power of IFA populations in detecting changes in both annual and interannual climate variability.

We hope you find this algorithm useful! Read on for brief instructions on how to download and use QUANTIFA.

## Repository Structure
```
|
├── model
│   ├── QUANTIFA_Pac_v00.m        <--- Pacific Variant
│   └── QUANTIFA_nonPac_v00.m     <--- Atlantic and Indian Variant
│
├── data
│   ├── reanalysis
│   │   ├── TA_ORAS5_download.py  <--- Tropical Atlantic Dataset
│   │   ├── TP_ORAS5_download.py  <--- Tropical Pacific Dataset
│   │   └── TI_ORAS5_download.py  <--- Tropical Indian Dataset
│   │
│   └── example_data
│       └── QUANTIFA_TropPac_ExDataset.xlsx
│
└── README.md
```
## Downloading the Algorithm
QUANTIFA can be downloaded directly from this repository. The algorithm comes in two variants: a Pacific variant (QUANTIFA_Pac_vXX.m) and a non-Pacific variant (QUANTIFA_nonPac_vXX.m). Be sure to download the correct variant based on your region of interest.

## Downloading the ORA-S5 Dataset
QUANTIFA is designed to work with a subset of potential temperature data from the Ocean Reanalysis System 5 data assimilation (ORA-S5). There are three such datasets available to download in this repository (depending on what region you're most interested in): one for the tropical Pacific Ocean, one for the tropical Atlantic Ocean, and another for the tropical Indian Ocean. Each of these data files contains a three-dimensional gridded field of potential temperature datat at a 1° x 1° horizontal resolution for 75 depth levels (0 - 5902 m). Each grid box contains a 61-year time series of monthly mean potential temperature data extending from Jan 1958 – Dec 2018.

All of this data means the data files are pretty large (>2 GB!) and cannot be uploaded to GitHub directly. Instead, I have included a small python script for each dataset that will retrieve the data from the host server and reformat it for use with QUANTIFA. Just follow these steps:

**(1)** Download the python script that corresponds to the regional dataset you are interested in (XX_ORAS5_download.py, where XX is either 'TP' for tropical Pacific, 'TA' for tropical Atlantic, or 'TI' for tropical Indian). Save it to the same location as the QUANTIFA algorithm.

**(2)** Open a terminal window (mac) or a command script (windows). Next, navigate to where the .py file is saved within the terminal window. You can do this by simply typing ```cd``` into the command line, followed by the path that points to your working directory (e.g. ```cd Documents/Python/Scripts```) and hit enter. Then type ```python XX_ORAS5_download.py``` (replacing XX with TP, TA, or TI), and hit enter once again.

**(3)** Grab a coffee! After about half an hour, a .mat file containing the ORA-S5 data subset should be saved to your specified directory.

Please contact me if there is any trouble in retrieving these data files (rglaubke@marine.rutgers.edu). For any specific information regarding the ORA-S5 dataset, please refer to Zuo et al. (2019) (doi:10.5194/os-15-779-2019) or visit the Ocean Synthesis/Reanalysis Directory of the Integrated Climate Data Center: https://icdc.cen.uni-hamburg.de/daten/reanalysis-ocean/easy-init-ocean/ecmwf-oras5.html.

## Implementing the Algorithm
For the most part, QUANTIFA should be a simple plug-and-play style algorithm. The first section of code (following the model description) is an input window where the user can upload data and define input conditions. After defining these inputs, the algorithm should run smoothly. Small descriptions of each input parameter are included in the comments of the script for easy reference. You can also find a description of each input parameter in Table 1 of our published paper (see below for citation information).

This algorithm can be used in one of three ways:

**(1)** The user can run an exploratory analysis that estimates the sensitivity of IFA populations at a given location to changes in annual and interannual climate variability. This could be a useful tool for establishing whether a particular geographic location and foraminiferal species is suitable for reconstructing some climate oscillation of interest. To perfom this analysis, comment out the X and Y variables in the input window.

**(2)** The user can take a single IFA population (e.g. a "modern" population from a sediment core-top) and compare it against modern hydrographic variability from the ORA-S5 dataset. To perform this analysis, comment out the X variable in the input window, and assign the IFA population to the Y variable.

**(3)** The user can compare two IFA distributions against one another to identify any statistically signficiant differences. Additionally, the algorithm aids in interpreting any detected differences by comparing the user's results against modeled results from a broad range of potential paleoclimate scenarios. To perform this exercise, assign the reference population (e.g. a "modern" core-top population) to the X variable, and assign the comparison population (e.g. a glacial population) to the Y variable.

## Output Products
QUANTIFA generates three output products: (1) a contour plot illustrating the relative sensitivity of IFA populations at the specified location to changes in annual and interannual climate variability; (2) a quantile-quantile (Q-Q) plot comparing the two inputted IFA populations (or a single IFA population against the ORA-S5 data); and (3) a heatmap illustrating the proportion of signficiantly-different quantiles that align with each hypothetical paleoclimate scenario, where the darker colors indicate hypothetical paleoclimate scenarios that are most consistent with the user's data.

For examples of these figures and how to interpret them, please see the case studies detailed in the discussion of our paper (see citation information below).

## Example Datasets
The datasets we used in our paper to demonstrate the algorithm's utility are also available to download in this repository (QUANTIFA_TropPac_ExDataset.xlsx). These include IFA populations from the coast of New Caledonia (Schmitt et al., 2019), the Line Islands (White et al., 2018), and the eastern equatorial Pacific (Ford et al., 2015). Feel free to use these data to orient yourself with the algorithm (if you do, be sure to use the tropical Pacific variant of QUANTIFA and the TP_ORAS5 dataset). Please see our paper for proper references to these data.

## Citation Information
If this code was helpful to you in your own research, please cite our paper!

Ryan H. Glaubke, Kaustubh Thirumalai, Matthew W. Schmidt, and Jennifer E. Hertzberg (submitted). Discerning Changes in High-Frequency Climate Variability using Geochemical Populations of Individual Foraminifera. Submitted to Paleoceanography and Paleoclimatology. doi:XXXX

We would love to see all of the cool and interesting ways you choose to use this algorithm!

## Version History
### 00 - Beta Release (August 3, 2020)
All code and documentation posted online for manuscipt review purposes.
