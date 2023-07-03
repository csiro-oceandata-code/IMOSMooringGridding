# IMOSMooringGridding
Code to collate IMOS Deep Water Mooring data into a gridded product, one for each mooring location.

The 'East Australian Current individual mooring gridded product - daily and depth gridded' produced by CSIRO is created using Matlab code. The dataset and full details is available from:

Sloyan, Bernadette; Cowley, Rebecca; Chapman, Chris (2021): East Australian Current individual mooring gridded product - daily and 10m depth gridded. v10. CSIRO. Data Collection. https://doi.org/10.25919/kr1g-ew82.

In addition, the gridded product can be accessed via the IMOS THREDDS Server:
http://thredds.aodn.org.au/thredds/catalog/IMOS/DWM/DA/CSIRO_gridded_all_variables/catalog.html

## **Product Build Overview:**
A full overview of the procedure is available from the CSIRO DAP link above, in the Documents section. The following is a basic overview of the order of code use. After initial setup step 1, control is from each individual mooring deployment code under the folders 'matlab\_EAC\_depyyyy\_yy', where yyyy\_yy represents the years of each deployment. Each deployment has a .m file in the format 'mooring\_EACdddd.m' where 'dddd' represents the depth of the mooring deployment. Settings such as file locations can be adjusted in each of the individual mooring deployment codes.

### Initial stacking process:

1. Download FV01 files from the AODN for each mooring, all deployments
2. For each mooring deployment:

    a. Create a csv file with essential metadata (one file per deployment), with information for all moorings (make_insinfo_fromdeploymentcsv.m). Can be done from the IMOS deployment database, or could equally be done from the FV01 files by reading the global attributes.
  
    b. Read all the FV01 instrument files and save as a Matlab structure (one file per instrument) (get_imosdata.m). File name is the instrument serial number. 
  
    c. Use individual mooring code ('mooring\_EACdddd.m') to implement the stacking routine (stack_mooring.m). Saves individual \*.mat files for each mooring, each deployment.

### Filtering and gridding process:
The stacked files created as outlined in the previous section are now used to create the final product for each mooring. This final product is one netcdf file per mooring, and each contains the variables UCUR, VCUR, TEMP and PSAL on a daily time grid (covering 2012 to 2021) and regular depth grid.

grid_filter_mooring_EAC.m








