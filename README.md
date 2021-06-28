# Mission ML Data Ready

Generate and temporally sync SoHO and/or SDO Mission image products to make a standardized machine-learning-ready data set.

## Flowchart 

[[https://github.com/cshneider/soho-ml-data-ready/blob/master/SoHO_SDO_Pipeline_Steps1_2.jpeg| alt=Flowchart]]

## Space Weather Application

[[https://github.com/cshneider/soho-ml-data-ready/blob/master/SoHO_SDO_Pipeline_Step3.jpeg| alt=CNN]]

## Instructions

1. wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
2. bash Anaconda3-2020.11-Linux-x86_64.sh
3. source ~/.bashrc
4. conda -V
5. conda create -n Give\_A\_Name\_for\_Your\_Virtual\_Environment python=3.6 anaconda
6. conda install --file Path\_to\_Mission\_Requirements\_File/Mission\_Requirements.txt

## Mission_Data_Gen.py

First generate the data that has the proper size and time cadence as specified by the user and remove any data with missing pixel values. 

### Usage and examples

Mission\_Data\_Gen experiment parameters:

| Input 			   | Description |
| ---------------------| ---------------------------------------------------- |
| -h, --help           |	Show this help message and exit. |
| --date_start         |	yyyy-mm-dd, 1996-01-01 is earliest start. |
| --date_finish        |	yyyy-mm-dd, 2011-05-01 is recommended latest finish. |
| --image_size_output  |	Image size (e.g., 128x128). |
| --time_window	   | Integer time step in hours |
| --flag 			   | Resize strategy. Choose from either "subsample", "interp", "minpool", or "maxpool". |
| --home_dir           | Home directory, e.g., "/home/user/Documents/", need "/" in the end. |
| --products           | Product types. Enter all the following or a subset thereof, in any order, seperated by commas. For the SoHO mission: "EIT195, MDI\_96m, LASCO\_C2, LASCO\_C3, EIT171, EIT304, EIT284" and for the SDO mission: "HMI_720s, AIA94, AIA131, AIA171, AIA193, AIA211, AIA304, AIA335, AIA1600, AIA1700, AIA4500". |
| --fits_headers       | Include header metadata in individual FITS files? Y/y or N/n. Applies primarily to MDI and HMI data from JSOC. Required argument. Faster download without header metadata. | 
| --lev1_LASCO		   | Whether to use level 1 LASCO C2 and C3? Y/y or N/n. If no, then level 0.5 LASCO will be used. Required argument. |
| --email 		   | User's email. Required by DRMS for the JSOC Client in order to have DRMS obtain the calibrated MDI products and all SDO mission products. |
| --mission 		   | Choose from 'SOHO' or 'SDO'. |

__Note:__ 
************
The framework does not allow for a mixture of mission products (e.g., MDI and HMI can not be entered together nor can AIA and EIT be entered together) and the products entered must be correctly paired with their respective missions (i.e., MDI, EIT, LASCO C2, and LASCO C3 for SoHO and HMI, AIA for SDO).
************

__Note:__ 
************
You must have an email address registered with JSOC before you are allowed to make a request. 
See this to register your email address:
http://jsoc.stanford.edu/ajax/register_email.html
************

Example usages:

```python
1. python nohup Mission_Data_Gen.py --products='EIT195, MDI_96m, LASCO_C3' --date_start='1996-01-01' --date_finish='2011-05-01' --image_size_output=128 --time_window=6
--flag=subsample --home_dir=/home/USER/ --fits_headers=N --lev1_LASCO=N --email=USER@Domain --mission=SOHO > LOG.log
```

```python
2. python nohup Mission_Data_Gen.py --products='MDI_96m' --date_start='1999-04-04' --date_finish='1999-04-06' --image_size_output=128 --time_window=6
--flag=subsample --home_dir=/home/USER/ --fits_headers=N --lev1_LASCO=N --email=USER@Domain --mission=SOHO > LOG.log
```

### Output

Example output for MDI_96m with subsample resize strategy to arrive at a final image size of 128x128:
- /home/USER/1999-04-04-00:00:02_to_1999-04-06-22:24:02\_MDI\_96m\_subsample\_6\_LASCOlev1-N\_SOHO_128.h5 --> all .fits files found, chronologically ordered with start time (first slice of cube) to finish time (last slice of cube). nomenclature contains exact time of start and finish .fits files composing the cube together with the product type of the cube, the downscaling strategy and the final image dimension.
- /home/USER/1999-04-04_to_1999-04-06\_MDI\_96m\_times\_subsample\_6\_LASCOlev1-N\_SOHO_128.csv --> contains initial times of all .fits files; may have duplicate times present.
- /home/USER/MDI_96m --> contains all *.fits files. All .fits files are unique.

In case of this time range being run as two seperate time ranges as: --date_start='1999-04-04' --date_finish='1999-04-05' followed by --date_start='1999-04-05' --date_finish='1999-04-06' one would obtain:
- /home/USER/1999-04-04-00:00:02\_to\_1999-04-05-20:48:02\_MDI\_96m\_subsample\_6\_LASCOlev1-N\_SOHO_128.h5 --> individual HDF5 file for that run
- /home/USER/1999-04-04-00:00:02\_to\_1999-04-06-03:12:02\_MDI\_96m\_subsample\_6\_LASCOlev1-N\_SOHO_128.h5 --> individual HDF5 file for that run
- /home/USER/1999-04-04-00:00:02\_to\_1999-04-06-22:24:02\_MDI\_96m\_subsample\_6\_LASCOlev1-N\_SOHO_128.h5 --> combined HDF5 file from both days when second run has finished as this picks up all the .fits files present in the folder.
- /home/USER/1999-04-04_to_1999-04-05\_MDI\_96m\_times\_subsample\_6\_LASCOlev1-N\_SOHO_128.csv --> individual CSV file for that run
- /home/USER/1999-04-05_to_1999-04-06\_MDI\_96m\_times\_subsample\_6\_LASCOlev1-N\_SOHO_128.csv --> individual CSV file for that run
- /home/USER/MDI_96m --> all .fits files from both runs.

Example name of a .fits file: /home/USER/MDI\_96m/SOHO\_MDI\_96m\_19990406031202_128.fits. The date and time information is combined in the file name (i.e., 1999-04-06 03:12:02 --> 19990406031202 ).
NOTE: CSV files pruduced from a split time range are not merged. CSV files are provided as more of a check for the USER. They are not used by any successive programs as all the time information is contained in the individual nomenclature  of each .fits file.

### Detailed program description

This program focuses on data products obtained from the SoHO (Solar and Heliospheric Observatory) and Solar Dynamics Observatory (SDO) missions:

* MDI (Michelson Doppler Imager) @ 96 minute cadence
* LASCO (Large Angle and Spectrometric Coronagraph): C2 (1.5 - 6 solar radii) and C3 (3.7 to 30 solar radii)
* EIT (Extreme ultraviolet Imaging Telescope) @ 171, 195, 284, 304 Angstroms)

and

* HMI (Helioseismic and Magnetic Imager) @ 12 minute cadence
* AIA (Atmospheric Imaging Assembly) See `Including Solar Dynamics Observatory (SDO) NASA mission data products' section here below on more details.

This program returns up to 7 folders, as specified by the user, which are named after their respective products that are queried with SunPy's Federated Internet Data Obtainer (Fido) and the Data Record Management System (DRMS) Python package. The DRMS is a system that was developed by the Joint Science Operation Center (JSOC).

For calibrated MDI images, DRMS queries Stanford University's JSOC.
These have the JSOC series 'mdi.fd\_M\_96m\_lev182'.

Each folder contains .fits (FITS = flexible Image Transport System) files which have been sieved for:

* appropriate image size
* filtered for image integrity (e.g., an absence of thresholded value of pixels in the original image)
* time stepped according to the user's input

Next to these product folders, there are compressed (gzipped) Hierarchical Data Format version 5 (HDF5) data cubes which are the, respectively, chronologically stacked data of all of the respective folder's .fits files that correspond to the .csv (comma seperated value) files generated with all the remaining good times of the downloaded .fits files.

In case of an interuption or running the program seperately on successive time periods, the user can choose the last day of time entry as the point to resume the program.
Program is designed to ensure that the same exact times are picked up for the same entered time periods so the results are exactly reproducible.

A log file would contain the names of files with holes or missing data.  
Clickable urls are available for EIT and LASCO images that would enable the user to download such images.
For MDI calibrated images from JSOC obtained with DRMS, in order to view those images reported in the log file as having holes:

```python
out_dir = 'PATH_SOME_DIRECTORY'
hole_MDI = client.export("A_NAME_IN_THE_MDI_96m holes_list_THAT'S_IN_YOUR_LOG_FILE")
hole_MDI_file = hole_MDI.download(out_dir)
```
 
The log file will also contain url links to LASCO C2 and LASCO C3 images that are probably comets or planet transients, based on applying the Probabilistic Hough Transform (https://scikit-image.org/docs/dev/auto_examples/edges/plot_line_hough_transform.html) on an edge-filtered image, that was obtained from using the Canny algorithm, to detect straight lines which are signatures of comets. Horizontal straight lines are signatures of planet transients and a subset of comets and so we use our filter to detect and remove C2 and C3 images containing planetary transients which are expected to affect the machine learning system. Comet detection would require some further rules for identifying lines which are not only horizontal and which are not due to concentric bright rings on the coronograph itself. Cosmic rays (CR) can also leave traces that have a small apparent size and resemble some planetary transient signatures on the C2 and C3 instruments. We apply the Laplacian of Gaussian (https://scikit-image.org/docs/dev/auto_examples/features_detection/plot_blob.html) on downsampled C2 and C3 images to detect blobs which correspond to CRs. 

Please check the Planet_Comet_examples folder provided for examples together with the Planet_Comet_Detector.ipynb Jupyter notebook.
Note that the FITS images in the Jupyter notebook plot are currently the vertically flipped images of those rendered by the ds9 program that reads FITS files.

__Note__:
Query Quota:
* SDAC has limit of 10k returns per query
* DRMS has limit of 100 GB per export request
As a result of these query quotas, an internal time increment of 60 days is used for all SDAC queries and for DRMS SDO HMI queries. HMI files are each 14 MB in size compared to MDI files at 1.4 - 1.6 MB each. 
For DRMS SDO AIA queries, with the exception of AIA4500 which has a cadence of once every hour and a file size of 14 MB, the time increment is set to 1 day due to the 12s and 24s cadences of the EUV and UV image products, respecitvely. The EUV and UV file sizes range from 7 - 12 MB. 


Note:
Images from earlier years can appear to have different intensities across the disk as in the case of EIT195 for 1996-01-15-20:51:47.
Perhaps use images from 1999 onwards for all product types.

Two time bars show progress: one time bar is for the completion of a given product type being processed and the other time bar is the progress in terms of all products entered.

SoHO mission data products can be obtained from VSO as follows: for MDI: 1996.05.01 − 2011.04.12, for LASCO: 1995.12.08 till present, for EIT 1996.01.01 to present.
The SDO (Solar Dynamics Observatory) mission provides higher resolution and higher cadence data products from 2010.05.12 AIA (Atmospheric Imaging Assembly) and from 2010.04.30 EVE (Extreme Ultraviolet Variability Experiment) together replaced EIT, from 2010.04.08 HMI (Helioseismic and Magnetic Imager) replaced MDI.
The advantage of using SoHO data is that it has basically covered solar cycles 23 and 24 with all of its products and continues into cyle 25 with most of its products.

Individual product benchmarks for querying all 7 data products with a time window of 6 hours and time span of 01.01.1999 - 12.31.2010.

| SoHO Product | No. Files | Time (hrs) | holes | planetary transients in C2 and C3 removed|
| ------------ | --------- | ---------- | ------| -------------------------------------------|
| LASCO_C2 | 15383 | ~10* | 489 | 7511  |
| LASCO_C3 | 13354 | ~18* | 699  | 23901  |
| MDI_96m JSOC  | 15456 | ~25** | 66  | |
| MDI_96m SDAC | 15428 | ~7*** | 290  | |
| EIT195   | 14439 | ~6 | 1391  | |
| EIT171   | 10234 | ~4 | 1045  | |
| EIT284   | 9601  | ~4 | 895  | |
| EIT304   | 9439  | ~4 | 936  | |

Total: ~58 hours

*These increased run-times are a result of using the Probabilistic Hough Transform for transient detection in LASCO C2 and LASCO C3 via scikit-image's probabilistic\_hough_line().  
OpenCV is faster with its cv2.HoughLinesP() implementation. 
However, OpenCV is not consistent with the current virtual environment settings which might break some of the existing code. 
An experienced user can investigate which lines of code (if any) would need to be updated if OpenCV is to be used instead.
Without any filters, LASCO\_C2 had 18.737 files and took ~7 hrs and LASCO\_C3 had 18.196 files and took ~5 hrs.

**The current MDI\_96m images are obtained from JSOC with DRMS which has a longer run time than using wget to obtain MDI_96m images from SDAC - these are no longer available.
The reported time of ~25 hrs is with the option of not having extended FITS headers included. Including this metadata will contribute to a longer run time. 
A better solution is to retroactively seed the MDI and HMI data cubes with their respective metadata. That is selecting 'N/n' for `fits_headers' and after the datacube has been made and the 
csv file with the times generated, then to use the `Retroactive metadata seeding' Jupyter notebook or eponymous script to include the metadata as an HDF5 attribute. The original HDF5 data cube
is copied and not gzipped so that the metadata can be efficiently added to it.

***Previously, used wget to obtain MDI_96m from VSO SDAC with 18.426 files and took ~7 hrs for 01.01.1996 - 01.05.2011.


__Note__:
2011-05-01 is recommended as the latest finish date only in the event that the user wishes to include MDI images in the analysis.
As described above, EIT and LASCO products continue to the present day.

With all .fits files present, a trick if necessary, to form the HDF5 data cube is to run the program on one date (start and finish date set equal) for which there was no data and the cube will be formed from all pre-existing .fits files in the folder.

If interuption occurs (e.g., SSL connection lost from using wget) just enter the next day after the last .fits file in the folder and the prev_time_resumer function will fill in the previous day.
When the program has been interuped and one is resuming on the next day, there will be another .csv file generated accoriding to the new date range that is input.
Since an interuption did occur, not all of the times of the .fits files will be contained in the .csv file. However, all the actual .fits files will be there.
The function product_retriever currently does not perfectly handle non-zero exit status as it retries once after 15 minutes and if this second try doesn't work an error is obtained which ends the program.
However, this method appears to be sufficiently robust.

__Note__: As an alternative to wget one can also use include Fido.fetch().

__Including Solar Dynamics Observatory (SDO) NASA mission data products__
Since our software pipeline uses SunPy Fido, other mission data can be readily obtained from Fido clients such as SDAC, JSOC, EVE, GONG, and other detauled under 'Fido clients' 
https://docs.sunpy.org/en/stable/guide/acquiring_data/fido.html#fido-guide. 
We use DRMS with JSOC to obtain all SDO mission products.
The SDO HMI Line-of-Sight Magnetograms, hmi.M_720s, with a cadence of 720 sec are used in this pipeline. 
To mention, for completenss, JSOC also has hmi.M_45s at a cadence of 45 sec.
SDO AIA seven EUV filters, series aia.lev1_euv_12s, are at a cadence of every 12 sec. EUV wavelengths are at 94, 131, 171, 193, 211, 304, 335 Angstrom. 
SDO AIA two UV band passes, series aia.lev1_uv_24s, are at a cadence of every 24 sec. UV wavelengths are at 1600 and 1700 Angstrom.
SDO AIA continuum data at 4500 Angstrom, series aia.lev1_vis_1h, is at a cadence of 1 hour.
For SDO AIA, both 'images' and 'spikes' (for calibration) are available and so 'images' are explicitely specified.

__Note__: 
It is suggested not to use the explicit export method with "(method='url', protocol='fits')" unless absolutely necessary to havw the FITS header metadata as this slows down the download and JSOC can stop responding after a certain number of export requests have been exceeded. Without the explicit "(method='url', protocol='fits')", there will be no 'hard' export requests for JSOC.


| SDO Product | ETA* (hrs)| 
| -------------| --------|
| HMI 		| ~100    |
| AIA visible  | ~50     |
| AIA UV**  	| ~50     |
| AIA EUV** 	| ~50     |

*Estimated time of arrival for SDO mission products from 2010-2020 at a cadence of every 6 hours with no explicit call to include FITS metadata.
**ETA per AIA EUV/UV product. These estimates are from a combination of the file size and high natural cadence (12s/24s) coupled with the 100 GB DRMS export limit. 
This previosuly resulted in a maximum internal time step of 1 day rather the 60 day time step used for the SOHO mission and for SDO HMI and SDO AIA 4500. 
Now, using @{time_window}h in DRMS to shorten the export.  

## Retroactive_metadata_seeding.py

Takes in original MDI and HMI data cube obtained without using the FITS protocol from JSOC and retroactively obtains metadata from JSOC and seeds the cube with these keywords as an attribute to HDF5.
Also takes in the corresponding csv times file in order to fetch the corresponding metadata usign DRMS query. 
JSOC can take 15 minutes per FITS file seen using the setup here. 

Retroactive\_metadata\_seeding parameters:

| Input 			     | Description |
| ---------------------  | ---------------------------------------------------- |
| -h, --help             | Show this help message and exit. |
| --image\_size\_output  | Image size originally used (e.g., 128x128). |
| --cube\_dir            | MDI or HMI cube path, e.g., "/home/user/Documents/", need "/" in the end. |
| --cube\_name	          | Name of .h5 cube, e.g., 1999-02-02-16:00:00\_to\_2011-01-01-00:00:00\_MDI\_96m\_subsample\_6\_LASCOlev1-N\_SOHO\_128.h5 |
| --product              | MDI\_96m or HMI\_720s |
| --mission 		     | Choose from 'SOHO' or 'SDO'. |

### Input 
- 1999-02-02-16:00:00_to_2011-01-01-00:00:00\_MDI\_96m\_subsample\_6\_LASCOlev1-N\_SOHO\_128.h5
- 1999-01-01\_to\_2010-12-31\_MDI\_96m\_times\_subsample\_6\_LASCOlev1-N\_SOHO\_128.csv
### Output
- 1999-02-02-16:00:00\_to\_2011-01-01-00:00:00\_MDI\_96m\_subsample\_6\_LASCOlev1-N\_SOHO\_128\_metadata.h5

__Metadata__
The FITS header metadata has been transferred to the HDF5 files and has been updated to account for the scale factor downsampling for several variables:
CDELT1, CDELT2, CRPIX1, CRPIX2 are always updated when they are present. 
RSUN\_OBS, R\_SUN, X0, Y0, and CROP_RAD or SOLAR\_R is updated when present along with other pixel dependent quantities contained in downsample\_header and downsample\_header\_local in 
Mission\_utility/\__init\__.py

All seven data cubes now contain metadata and have the \_metadata.h5 and metadata_sync.h5 suffixes at the end of each name corresponding to the cube generation and syncronization steps, respectively.
The suffix \_retroactivemetadata.h5 is appended to the end of the MDI or HMI data cube after the Retroactive\_metadata\_seeding.py script is used.
Running the retroactive script on the SoHo MDI data takes about ~15 minutes with 10 minutes for loading the Pandas dataframe of 83700 entries fetched from JSOC and 5 minutes with creating the metadata dictionary.
An accompanying Jupyter notebook, Header\_Reader\_for\_Gen\_and\_Sync\_datacubes.ipynb, is also provided for reading metadata from datacubes from Gen and Sync steps.

The cube\_metadata has a dictionary like appearance and can be easily converted to an actual dictionary with key,value pairs.
The FITS standard keywords appearing in cube\_metadata have an extra `_SliceNumber' appended to indicate what slice they belong to. 
So CRPIX1\_0 would correspond to CRPIX1 of the first slice and so forth.

It can be accessed in the following way:

```python
import h5py
cube = h5py.File('file_name.h5','r')
data = cube[f'{Mission_Product_name}_{Mission}_{Output_dimension}'][:]
cube_json_serialization = list(cube.attrs.items())
```

Using JSON to decode the metadata (which has been JSON serialized):
```python
with h5py.File('file_name.h5', 'r') as hfile:
    metadata = json.loads(hfile[f'{base}_{mission}_{image_size_output}_metadata'][()])
```

For the datacubes corresponding to Tables 2 and 3 of the paper, the run times to sync are:

| SoHO Experiment (sync step)		| Time (hrs)| 
| -------------		| --------  |
| MDI  				| ~2     	  |
| MDI, EIT 195, LASCO C2 | ~5        |
| All 7 products  		| ~4        |


## Mission_Product_Sync.py
This is the companion script to \Mission\_Data\_Gen.py to synchronize the times between the specified data products once they have all been downloaded from the VSO and pre-processed.
The program outputs an HDF5 compressed (gzipped) data cube and an accompanying .csv file, both containing the "sync" keyword.
The data is the chronologically stacked data of all of the respective folder's .fits files whose time stamps overlap within the specified time_step.
All data is cast as int16 with interger values in the range [-32767, 32767] for all pixels. This should be sufficient to deal with all pixel values present in any product.
GPUs perform best on images when lower precision is assigned to pixel values.
The corresponding synchronized times are contained in the companion .csv files. Times across products are synchronized when an individual time from a specific product comes within a given time range and such a time exists for all products specified.
The nomenclature reflects original time step and current time step in case have subsampled by an integer multiple of the original time step (e.g., 12 hrs when had 6 hrs originally) but not smaller than the original time step (e.g., can't choose 3 hrs if had 6 hrs originally).
Two time progress counters show the elapsed time in real time.

The algorithm works as follows:

Firstly, two options are provided, which are selected with a yes/no:
- If one just has downlaoded the HDF5 cubes and .csv files without running \Mission_Data_Gen.py locally then the data from HDF5 cubes will be uptaken along with the corresponding times from .csv file for each product specified.
- If have both the HDF5 data cubes along with the corresponding .fits files then the .fits files will be used to obtain both the data and the times since the time stamps are in the names of the .fits files themselves.

Next:

- Find product with shortest time list and base the other times on this product by moving along its time axis.
- Center the time step on this product and set the time range to +/- half the time step.
- Take the closest times from all products to that time. If equal times are found for a product, then always choose the first such time.
- If any time is missing from a product for a particular point on the time axis, then disregard all other product times obtained for that time axis point and move on to next time. Hence, the more products specified, the harder it is for a time axis point to be shared among all these products.

Although the dimensions across the different products are checked to be the same, as described above, the flag (resize strategy used) is left open since the user may have wanted to combine images of the same dimension but obtained with different resize strategies across the different products.

Mission\_Product\_Sync experiment parameters:

| Input 			   | Description |
| ---------------------| ---------------------------------------------------- |
| -h, --help           | Show this help message and exit. |
| --date_start	        | yyyy-mm-dd, 1996-01-01 is earliest start. |
| --date_finish	   | yyyy-mm-dd, 2011-05-01 is recommended latest finish. |
| --time_step		   | Integer time step in hours. |
| --home_dir           | Home directory, e.g., "/home/user/Documents/", need "/" in the end. |
| --option 		   | Are *.fits files present? Y/N or y/n. |
| --products           | Product types to synchronize. Enter all the following or a subset thereof, in any order, seperated by commas: "EIT195, MDI_96m, LASCO_C2, LASCO_C3, EIT171, EIT304, EIT284". |
| --mission 		   | Choose from 'SOHO' or 'SDO'. |

Example usage:
If one had run Mission\_Data\_Gen.py with the following inputs: --products='MDI\_96m, LASCO\_C3, EIT284, EIT195, LASCO\_C2, EIT304, EIT171', --date_start='1996-01-01', --date_finish='2011-05-01', and time_window=6, (the --flag and --image_size_output are not important in this example) then one could do the following to sync a subset of the original products within a subset of the original time range and with a coarser time sampling of 12 hrs instead of 6 hrs:
```python
1. python nohup Mission_Product_Sync.py --date_start=1999-01-01 --date_finish=2011-05-01 --time_step=12 --home_dir=/home/USER/ --option=Y --products='MDI_96m, EIT195, LASCO_C2' --mission='SOHO' > LOG.log
```

Example output:

In addition to the product folders, .h5 files, and .csv files already present, the following new products would be produced:
- 1999-01-01_to_2011-05-01\_MDI\_96m\_SOHO\_3products\_6\_12\_128\_metadata\_sync.h5
- 1999-01-01_to_2011-05-01\_EIT195\_SOHO\_3products\_6\_12\_128\_metadata\_sync.h5
- 1999-01-01_to_2011-05-01\_LASCO\_C2\_SOHO\_3products\_6\_12\_128\_metadata\_sync.h5
- 1999-01-01_to_2011-05-01\_MDI\_96m\_SOHO\_3products\_6\_12\_128\_times\_sync.csv
- 1999-01-01_to_2011-05-01\_EIT195\_SOHO\_3products\_6\_12\_128\_times\_sync.csv
- 1999-01-01_to_2011-05-01\_LASCO\_C2\_SOHO\_3products\_6\_12\_128\_times\_sync.csv
with the following set as well which is synced with the LASCO difference images to subtract the Fcorona:
- 1999-01-01_to_2011-05-01\_MDI\_96m\_SOHO\_3products\_Fcorona\_6\_12\_128\_metadata\_sync.h5
- 1999-01-01_to_2011-05-01\_EIT195\_SOHO\_3products\_Fcorona\_6\_12\_128\_metadata\_sync.h5
- 1999-01-01_to_2011-05-01\_LASCO\_C2\_SOHO\_3products\_Fcorona\_6\_12\_128\_metadata\_sync.h5
- 1999-01-01_to_2011-05-01\_MDI\_96m\_SOHO\_3products\_Fcorona\_6\_12\_128\_times\_sync.csv
- 1999-01-01_to_2011-05-01\_EIT195\_SOHO\_3products\_Fcorona\_6\_12\_128\_times\_sync.csv
- 1999-01-01_to_2011-05-01\_LASCO\_SOHO\_C2\_3products\_Fcorona\_6\_12\_128\_times\_sync.csv

Run time: ~10 min.

These final *\_metadata\sync.h5 data cubes are now ready for input into an ML architecture.
