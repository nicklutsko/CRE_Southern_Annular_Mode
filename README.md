This package of scripts and data was used in the analysis of “Revisiting Cloud Radiative Heating and the Southern Annular Mode” by Wall et al. The scripts in this
package were used for loading and gridding the CloudSat satellite data and for performing the main statistical analysis that relates atmospheric cloud radiative heating
to the Southern Annular Mode. To highlight the analysis methods as clearly and concisely as possible, we have removed code that plots the results.

The user will need to take two preliminary steps before running the scripts. First, all scripts are written for Matlab, so the user will need to have access to this analysis
program. Second, the scripts assume that the user has downloaded CloudSat data products 2B-FLXHR-LIDAR, 2B-GEOPROF-LIDAR, and ECMWF-AUX, which are freely
available at https://www.cloudsat.cira.colostate.edu/data-products. For each data product, all data files need to be saved in a single directory (i.e. the data files from a
given data product should not be separated into different subdirectories based on the time of observation). The user will need to change the code in the scripts so that the
path to the data files matches their machine. Comments are included in the scripts where these changes need to be made.

The four main scripts in the top-most directory are:
• load_CloudSat_height_coordinates_480m.m: This is the main script that loads the
CloudSat data and creates a gridded dataset of zonal-mean pentad-mean values.
This script should be run first.
• calculate_overcast_heating_rates.m: This script computes the average atmospheric
LW cloud radiative heating rate 𝑅 over each cloud regime (L, M, LH, etc.). This script
should be run second.
• “compute_LW_CRH_climatology_and_SAM_regression.m”: This script computes
the long-term mean of 𝑅 and cloud fraction, and it computes regressions of 𝑅 and
cloud fraction against the SAM index, 𝑠.
• “LW_CRH_decomposition.m”: This script computes the contribution of each cloud
regime to 𝑑𝑅/𝑑𝑠.

Finally, the scripts use 𝑠 data that are saved in the file “ERA5_AM_19790101_20181231.nc” in the “data” directory. These data were computed from ERA5 reanalysis wind fields following the method of Thompson and
Woodworth (2014; https://doi.org/10.1175/JAS-D-13-0185.1). Shim Yook and David Thompson computed and shared these data.
