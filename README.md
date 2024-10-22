# QUANTIFA: Quantile Analysis of Temperature using Individual Foraminiferal Analyses
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7775163.svg)](https://doi.org/10.5281/zenodo.7775163)

> **Note**
> The latest version of QUANTIFA (2.0.0) includes big changes to the model's core code, but is currently only available for the Tropical Pacific. The publication associated with these changes has just been published:
> Glaubke, R. H., M. W. Schmidt, J. E. Hertzberg, L. G. Ward, F. Marcantonio, D. Schimmenti, K. Thirumalai (2024). Divergent ENSO Responses to Northern Hemisphere Stadials during the Last Deglaciation. Geophysical Research Letters, 51(12), https://doi.org/10.1029/2023gl107634

QUANTIFA is a user-friendly proxy system model for Individual Foraminiferal Analyses (IFA) that combines routines for modeling the sensitivity of IFA populations to changes in climate variability with tools for processing, plotting, and interpreting IFA-Mg/Ca data.

We hope you find this algorithm useful! Read on for brief instructions on how to download and use QUANTIFA.

## Repository Structure
```
│
├── data
│   ├── example_data
│   │   └── QUANTIFA_TropPac_ExDataset.xlsx
│   │ 
│   └──reanalysis
│       ├── TA_ORAS5_download.py  <--- Tropical Atlantic Reanalysis Dataset
│       ├── TI_ORAS5_download.py  <--- Tropical Indian Reanalysis Dataset
│       └── TP_ORAS5_download.py  <--- Tropical Pacific Reanalysis Dataset
│
├── model
│   ├── QUANTIFA_Pac_v200.m       <--- QUANTIFA: Pacific Variant (version 2.0.0)
│   └── QUANTIFA_nonPac_v101.m    <--- QUANTIFA: Atlantic and Indian Variant (version 1.0.1)
│
└── README.md
```
## Downloading the Algorithm
QUANTIFA can be downloaded directly from this repository. The algorithm comes in two variants: a Pacific variant (QUANTIFA_Pac_v###.m) and a non-Pacific variant (QUANTIFA_nonPac_v###.m). Be sure to download the correct variant based on your region of interest.

## Downloading the ORA-S5 Dataset
QUANTIFA is designed to work with a subset of potential temperature data from the Ocean Reanalysis System 5 data assimilation (ORA-S5). There are three such datasets available to download in this repository (depending on what region you're most interested in): one for the tropical Pacific Ocean, one for the tropical Atlantic Ocean, and another for the tropical Indian Ocean. Each of these data files contains a three-dimensional gridded field of potential temperature data at a 1° x 1° horizontal resolution for 75 depth levels (0 - 5902 m). Each grid box contains a 61-year time series of monthly mean potential temperature data extending from Jan 1958 – Dec 2018.

All of this data means the data files are pretty large (>2 GB!) and cannot be uploaded to GitHub directly. Instead, I have included a small python script for each dataset that will retrieve the data from the host server and reformat it for use with QUANTIFA. Just follow these steps:

**(1)** Download the python script that corresponds to the regional dataset you are interested in (XX_ORAS5_download.py, where XX is either 'TP' for tropical Pacific, 'TA' for tropical Atlantic, or 'TI' for tropical Indian). Save it in the same location as the QUANTIFA algorithm.

**(2)** Open a terminal window (mac) or a command script (windows). Next, navigate to where the .py file is saved within the terminal window. You can do this by simply typing ```cd``` into the command line, followed by the path that points to your working directory (e.g. ```cd Documents/Python/Scripts```) and hit enter. Then type ```python XX_ORAS5_download.py``` (replacing XX with TP, TA, or TI), and hit enter once again.

**(3)** Grab a coffee! After a while, a .mat file containing the ORA-S5 data subset should be saved to your specified directory. (NOTE: this step can take some time if your machine's RAM is on the lower end (8 GB). Please be patient! I promise it's working.)

Please contact me if there is any trouble in retrieving these data files (rglaubke@marine.rutgers.edu). For any specific information regarding the ORA-S5 dataset, please refer to Zuo et al. (2019) (doi:10.5194/os-15-779-2019) or visit the Ocean Synthesis/Reanalysis Directory of the Integrated Climate Data Center: https://icdc.cen.uni-hamburg.de/daten/reanalysis-ocean/easy-init-ocean/ecmwf-oras5.html.

## Implementing the Algorithm
For the most part, QUANTIFA should be a simple plug-and-play-style algorithm. The first section of code (following the model description) is an input window where the user can upload data and define input conditions. After defining these inputs, the algorithm should run smoothly. Small descriptions of each input parameter are included in the comments of the script for easy reference. You can also find a description of each input parameter in Table 1 of our original publication (see below for citation information). You can find updates to some of these parameters in the methods of our upcoming publication (stay tuned).

This algorithm can be used in one of three ways:

**(1)** The user can run an exploratory analysis (without the need of any IFA data) that estimates the sensitivity of IFA populations at a given location/depth to changes in annual and interannual climate variability. This is a useful tool for establishing whether a particular location or foraminiferal species is suitable for reconstructing some climate oscillation of interest. To perfom this analysis, comment out the two IFA data variables ```X``` and ```Y``` in the input window.

**(2)** The user can take a single IFA population (e.g. a "modern" population from a sediment core-top) and compare it against modern hydrographic variability from the ORA-S5 dataset. To perform this analysis, comment out the ```X``` variable in the input window, and assign the IFA population to the ```Y``` variable.

**(3)** The user can perform IFA population comparisons to identify differences between two inputted IFA datasets, and build an interpretative framework (including flase positive tests and data-model consistency analysis) to aid in data interpretation. To perform this exercise, assign the reference population (e.g. a "modern" core-top population) to the ```X``` variable, and assign the comparison population (e.g. a glacial population) to the ```Y``` variable.

## Output Products
QUANTIFA generates four primary output products: (1) a conformity contour plot illustrating the relative sensitivity of IFA populations at the specified location and depth to changes in annual and interannual climate variability; (2) a quantile-quantile (Q-Q) plot comparing the two inputted IFA populations (or a single IFA population against the ORA-S5 data); (3) a matrix containing false positive rates (mean and SD) for each individual quantile; and (3) a data-model consistency map illustrating the proportion of signficiant quantiles that align with a suite of hypothetical paleoclimate scenarios.

For examples of these figures and how to interpret them, please see the case studies detailed in the discussion of our initial paper (see citation information below) or the discussion in an our upcoming publication (stay tuned). 

## Example Datasets
The datasets we used in the model application exercises in our paper are available to download in this repository (QUANTIFA_TropPac_ExDataset.xlsx). These include IFA populations from the coast of New Caledonia (Schmitt et al., 2019), the Line Islands (White et al., 2018), and the eastern equatorial Pacific (Ford et al., 2015). Feel free to use these data to orient yourself with the algorithm (if you do, be sure to use the tropical Pacific variant of QUANTIFA and the TP_ORAS5 dataset). Please see our paper for proper references to these data.

## Disclaimer
*"All models are wrong, but some are useful." - George E. P. Box*

I'm a pretty novice programmer. If you happen to stumble upon any bugs/errors in any of the algorithms stored in this repository, PLEASE bring them to my attention! Please also feel free to forward suggestions for making the code run more efficiently.

## Citation Information
If this algorithm was helpful to you in your own research, please cite us!

**Initial Release:**
Ryan H. Glaubke, Kaustubh Thirumalai, Matthew W. Schmidt, and Jennifer E. Hertzberg (2021). Discerning Changes in High-Frequency Climate Variability using Geochemical Populations of Individual Foraminifera. *Paleoceanography and Paleoclimatology*, *36*(2), e2020PA004065. https://doi.org/10.1029/2020PA004065.

**Latest Version (2.0.0)**
Stay tuned...

We would love to see all of the cool and interesting ways you choose to use this algorithm!

## Version History
### v2.0.0 (March 27, 2023)
- Manipulates ENSO amplitude and frequency using a new "ENSO variability" metric.
- Depth input parameter ('dep') now accepts a range of depth and associated weights—rather than a single depth horizon—to more accurately parameterize a species' depth distribution.
- Updates the seasonal bias input parameter ('seas') to weight the picking algorithm towards a species' season of production, rather than picking from that season exclusively.
- Includes a new barycenter calculation that reports the "center of gravity" coordinates of the data-model consistency maps.
- Fixes a bug where the seasonal climatology calculation was including ENSO event years.
- Features a progress bar to help keep track of those runs!
### v1.0.1 (February 24, 2020)
- fixes a typo in the nnz() function during input parsing.
### v1.0.0 - Full Launch (February 2, 2021)
First full launch of QUANTIFA, coinciding with the release of our manuscript.
### v0.0.1 - Beta Release, Again (November 18, 2020)
New release number generated to create a DOI through Zenodo. No content change between v0.0.0 and v0.0.1.
### v0.0.0 - Beta Release (August 3, 2020)
All code and documentation posted online for manuscipt review purposes.
