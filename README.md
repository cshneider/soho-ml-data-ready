# SOHO ML Data Ready

## SOHO_DATA_GEN.py
This program focuses entirely on data products obtained from the NASA SOHO (Solar and Heliospheric Observatory) mission: MDI (Michelson Doppler Imager) @ 96 minute cadence, 
LASCO (Large Angle and Spectrometric Coronagraph): C2 (1.5 - 6 solar radii) and C3 (3.7 to 30 solar radii), and EIT (Extreme ultraviolet Imaging Telescope) @ 171, 195, 284, 304 Angstroms).
This program returns up to 7 folders, as specified by the user, named after their resective products which are queried with SunPy's Fido (Federated Internet Data Obtainer), in this case, specifically from the VSO's (Virtual Solar Observatory) SDAC (Solar Data Analysis Center, NASA/Goddard) provider. 
Each folder contains .fits (FITS = flexible Image Transport System) files which have been sieved for: appropriate image size, filtered for image intergrity (e.g., an absence of thresholded value of pixels in the original image), and time stepped according to the user's input.
Next to these product folders, there are compressed (gzipped) HDF5 data cubes which are the, respectively, chronologically stacked data of all of the respective folder's .fits files that correspond to the .csv (comma seperated value) files generated with all the remaining good times of the downloaded .fits files. 
In case of an interuption or running the program seperately on successive time periods, the user can choose the last day of time entry as the point to resume the program. 
Program is designed to ensure that the same exact times are picked up for the same entered time periods so the results are exactly reproducible.  
A log file would contain the names of files with holes or missing data, including a clickable url path for these images that would enable the user to download such images. 
Due to the VSO limit of 10k returns per query, an internal time increment of 60 days is used. 

Note:
Images from earlier years can appear to have different intensities across the disk as in the case of EIT195 for 1996-01-15-20:51:47.
Perhaps use images from 1999 onwards for all product types.

Two time bars show progress: one time bar is for the completion of a given product type being processed and the other time bar is the progress in terms of all products entered.

SOHO mission data products can be obtained from VSO as follows: for MDI: 1996.05.01 − 2011.04.12, for LASCO: 1995.12.08 till present, for EIT 1996.01.01 to present. 
The SDO (Solar Dynamics Observatory) mission provides higher resolution and higher cadence data products from 2010.05.12 AIA (Atmospheric Imaging Assembly) and from 2010.04.30 EVE (Extreme Ultraviolet Variability Experiment) together replaced EIT, from 2010.04.08 HMI (Helioseismic and Magnetic Imager) replaced MDI. 
The advantage of using SOHO data is that it has basically covered solar cycles 23 and 24 with all of its products and continues into cyle 25 with most of its products.       

Invidual product benchmarks for querying all 7 data products with a time window of 6 hours and time span of 01.01.1996 - 01.05.2011.

| SOHO Product | No. Files | Time (hrs) |
| ------------ | --------- | ---------- | 
| LASCO_C2 | 18.737 | ~7 |
| LASCO_C3 | 18.196 | ~5 |
| MDI_96m  | 18.426 | ~7 |
| EIT195   | 17.521 | ~6 |
| EIT171   | 12.174 | ~4 |
| EIT284   | 11.438 | ~4 |
| EIT304   | 11.279 | ~4 |

Total: ~41 hours

SOHO_DATA_GEN experiment parameters:

| Input 			   | Description | 
| ---------------------| ---------------------------------------------------- |
| -h, --help           |	Show this help message and exit. |  
| --date_start         |	yyyy-mm-dd, 1996-01-01 is earliest start. |
| --date_finish        |	yyyy-mm-dd, 2011-05-01 is recommended latest finish. |
| --target_dimension   |	Image size (e.g., 128x128). |
| --time_window	   | Integer time step in hours |
| --flag 			   | Resize strategy. Choose from either "subsample", "interp", "minpool", or "maxpool". |
| --home_dir           | Home directory, e.g., "/home/user/Documents/", need "/" in the end. |
| --products           | Product types. Enter all the following or a subset thereof, in any order, seperated by commas: "EIT195, MDI_96m, LASCO_C2, LASCO_C3, EIT171, EIT304, EIT284". |

Note:
2011-05-01 is recommended as the latest finish date only in the event that the user wishes to include MDI images in the analysis. 
As described above, EIT and LASCO products continue to the present day. 

Example usages: 

```python 
1. python nohup SOHO_ML_DATA_GEN.py --products='EIT195, MDI_96m, LASCO_C3' --date_start='1996-01-01' --date_finish='2011-05-01' --target_dimension=128 --time_window=6 
--flag=subsample --home_dir=/home/USER/ > LOG.log
```

```python 
2. python nohup SOHO_ML_DATA_GEN.py --products='MDI_96m' --date_start='1999-04-04' --date_finish='1999-04-06' --target_dimension=128 --time_window=6 
--flag=subsample --home_dir=/home/USER/ > LOG.log
```

Example output for MDI_96m with subsample resize strategy to arrive at a final image size of 128x128:
- /home/USER/1999-04-04-00:00:02_to_1999-04-06-22:24:02_MDI_96m_subsample_6_128.h5 --> all .fits files found, chronologically ordered with start time (first slice of cube) to finish time (last slice of cube). nomenclature contains exact time of start and finish .fits files composing the cube together with the product type of the cube, the downscaling strategy and the final image dimension.
- /home/USER/1999-04-04_to_1999-04-06_MDI_96m_times_subsample_6_128.csv --> contains initial times of all .fits files; may have duplicate times present.
- /home/USER/MDI_96m --> contains all *.fits files. All .fits files are unique.

In case of this time range being run as two seperate time ranges as: --date_start='1999-04-04' --date_finish='1999-04-05' followed by --date_start='1999-04-05' --date_finish='1999-04-06' one would obtain: 
- /home/USER/1999-04-04-00:00:02_to_1999-04-05-20:48:02_MDI_96m_subsample_6_128.h5 --> individual HDF5 file for that run
- /home/USER/1999-04-04-00:00:02_to_1999-04-06-03:12:02_MDI_96m_subsample_6_128.h5 --> individual HDF5 file for that run
- /home/USER/1999-04-04-00:00:02_to_1999-04-06-22:24:02_MDI_96m_subsample_6_128.h5 --> combined HDF5 file from both days when second run has finished as this picks up all the .fits files present in the folder.
- /home/USER/1999-04-04_to_1999-04-05_MDI_96m_times_subsample_6_128.csv --> individual CSV file for that run
- /home/USER/1999-04-05_to_1999-04-06_MDI_96m_times_subsample_6_128.csv --> individual CSV file for that run
- /home/USER/MDI_96m --> all .fits files from both runs.

Example name of a .fits file: /home/USER/MDI_96m/SOHO_MDI_96m_19990406031202_128.fits. The date and time information is combined in the file name (i.e., 1999-04-06 03:12:02 --> 19990406031202 ).
NOTE: CSV files pruduced from a split time range are not merged. CSV files are provided as more of a check for the USER. They are not used by any successive programs as all the time information is contained in the individual nomenclature  of each .fits file.  

With all .fits files present, a trick if necessary, to form the HDF5 data cube is to run the program on one date (start and finish date set equal) for which there was no data and the cube will be formed from all pre-existing .fits files in the folder.

If interuption occurs (e.g., SSL connection lost from using wget) just enter the next day after the last .fits file in the folder and the prev_time_resumer function will fill in the previous day. 
When the program has been interuped and one is resuming on the next day, there will be another .csv file generated accoriding to the new date range that is input. 
Since an interuption did occur, not all of the times of the .fits files will be contained in the .csv file. However, all the actual .fits files will be there.
The function product_retriever currently does not perfectly handle non-zero exit status as it retries once after 15 minutes and if this second try doesn't work an error is obtained which ends the program.
However, this method appears to be sufficiently robust.

Note: As an alternative to wget one can also use include Fido.fetch().

## SOHO_PRODUCT_SYNC.py
This is the companion script to SOHO_DATA_GEN.py to synchronize the times between the specified data products once they have all been downloaded from the VSO and pre-processed.
The program outputs an HDF5 compressed (gzipped) data cube and an accompanying .csv file, both containing the "sync" keyword. 
The data is the chronologically stacked data of all of the respective folder's .fits files whose time stamps overlap within the specified time_step. 
All data is cast as int16 with interger values in the range [-32767, 32767] for all pixels. This should be sufficient to deal with all pixel values present in any product. 
GPUs perform best on images when lower precision is assigned to pixel values.
The corresponding synchronized times are contained in the companion .csv files. Times across products are synchronized when an individual time from a specific product comes within a given time range and such a time exists for all products specified.   
The nomenclature reflects original time step and current time step in case have subsampled by an integer multiple of the original time step (e.g., 12 hrs when had 6 hrs originally) but not smaller than the original time step (e.g., can't choose 3 hrs if had 6 hrs originally).
Two time progress counters show the elapsed time in real time. 

The algorithm works as follows:
 
Firstly, two options are provided, which are selected with a yes/no:
- If one just has downlaoded the HDF5 cubes and .csv files without running SOHO_DATA_GEN.py locally then the data from HDF5 cubes will be uptaken along with the corresponding times from .csv file for each product specified. 
- If have both the HDF5 data cubes along with the corresponding .fits files then the .fits files will be used to obtain both the data and the times since the time stamps are in the names of the .fits files themselves.

Next:

- Find product with shortest time list and base the other times on this product by moving along its time axis.
- Center the time step on this product and set the time range to +/- half the time step. 
- Take the closest times from all products to that time. If equal times are found for a product, then always choose the first such time.
- If any time is missing from a product for a particular point on the time axis, then disregard all other product times obtained for that time axis point and move on to next time. Hence, the more products specified, the harder it is for a time axis point to be shared among all these products.

Although the dimensions across the different products are checked to be the same, as described above, the flag (resize strategy used) is left open since the user may have wanted to combine images of the same dimension but obtained with different resize strategies across the different products.

SOHO_PRODUCT_SYNC experiment parameters:

| Input 			   | Description | 
| ---------------------| ---------------------------------------------------- |
| -h, --help           | Show this help message and exit. |  
| --date_start	        | yyyy-mm-dd, 1996-01-01 is earliest start. |
| --date_finish	   | yyyy-mm-dd, 2011-05-01 is recommended latest finish. |
| --time_step		   | Integer time step in hours. |
| --home_dir           | Home directory, e.g., "/home/user/Documents/", need "/" in the end. |
| --option 		   | Are *.fits files present? Y/N or y/n. |
| --products           | Product types to synchronize. Enter all the following or a subset thereof, in any order, seperated by commas: "EIT195, MDI_96m, LASCO_C2, LASCO_C3, EIT171, EIT304, EIT284". |

Example usage:
If one had run SOHO_DATA_GEN.py with the following inputs: --products='MDI_96m, LASCO_C3, EIT284, EIT195, LASCO_C2, EIT304, EIT171', --date_start='1996-01-01', --date_finish='2011-05-01', and time_window=6, (the --flag and --target_dimension are not important in this example) then one could do the following to sync a subset of the original products within a subset of the original time range and with a coarser time sampling of 12 hrs instead of 6 hrs:
```python
1. python nohup SOHO_PRODUCT_SYNC.py --date_start=1999-01-01 --date_finish=2011-05-01 --time_step=12 --home_dir=/home/USER/ --option=Y --products='MDI_96m, EIT195, LASCO_C2' > LOG.log
```

Example output:

In addition to the product folders, h5 files, and .csv files already present, the following new products would be produced:
- 1999-01-01_to_2011-05-01_MDI_96m_6_12_128_sync.h5
- 1999-01-01_to_2011-05-01_EIT195_6_12_128_sync.h5
- 1999-01-01_to_2011-05-01_LASCO_C2_6_12_128_sync.h5
- 1999-01-01_to_2011-05-01_MDI_96m_6_12_128_times_sync.csv
- 1999-01-01_to_2011-05-01_EIT195_6_12_128_times_sync.csv
- 1999-01-01_to_2011-05-01_LASCO_C2_6_12_128_times_sync.csv

Run time: ~8.5 min. (A bit slower than it takes light to reach Earth from the Sun). 
Each of the *sync.csv files contains 6.336 entries.

These final *sync.h5 data cubes are now ready for input into an ML architecture. 
