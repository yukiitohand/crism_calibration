# crism_toolbox
MRO CRISM data analysis toolkit

This toolbox serves as an MATLAB interface for handling MRO CRISM database (http://pds-geosciences.wustl.edu/missions/mro/crism.htm). You can download images into your local database organized in a same (or similar) way as the pds server. Resolving image locations and filenames can be simply done. 

## Installation
First you need to customize your setting/crismToolbox_default.json:

```json
{
    "localCRISM_PDSrootDir" : "",
    "localCATrootDir"       : "/../cat/CAT_v7_4/",
    "CAT_ver"               : "7.4",
    "local_fldsys"          : "pds_unified",
    "remote_fldsys"         : "pds_mro",
    "pds_unified_URL"       : "crism_pds_archive/",
    "pds_mro_URL"           : "pds-geosciences.wustl.edu/mro/",
    "LUT_OBSID2YYYY"        : "LUT_OBSID2YYYY_DOY.mat"
}
```
What you should specify is `localCRISM_PDSrootDir` where the local databse is created and `localCATrootDir` where CAT_ENVI/ is stored. The others do not need to be changed unless you want to a specific database structure. If you opt to use a different remote server, please contact me. With this setup, local database will be created at
```
[localCRISM_PDSrootDir]/crism_pds_archive/
```
Some more information for the setup file.
* `localCRISM_PDSrootDir`: root directory path for which the database will be downloaded.
* `lovcalCATrootDir`     : directory path for which CAT_ENVI is stored.
* `CAT_ver`              : version of the CAT
* `local_fldsys`         : database system (`'pds_mro'` and `'pds_unified'` is supported now)
* `remote_fldsys`        : database system (`'pds_mro'` is supported)
* `pds_unified_URL`      : root directory name or path for the folder system `'pds_unified'`
* `pds_mro_URL`          : root directory name or path for the folder system `'pds_mro'`
* `LUT_OBSID2YYYY`       : the name of the mat file for which yyyy_doy look up table is stored. The table comes with the toolbox, so you do not need to change.

Second, rename setting/crismToolbox_default.json to setting/crismToolbox.json. and run 
```
> crism_setup
```

## Basic Operations
Reading the image can be performed as simply as
```
> TRRIFdata = CRISMdata('HRL000040FF_07_IF183L_TRR3','');
```
If you haven't downloaded, you can download to the local database by 
```
> TRRIFdata.download(2);
```

