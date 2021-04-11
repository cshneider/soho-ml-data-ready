import numpy as np

import os
from os import listdir
from os.path import isfile, join

import shlex, subprocess

from skimage.transform import rescale
from skimage.measure import block_reduce

from dateutil import parser
from datetime import datetime, date, time, timedelta
import time

import drms
from sunpy.net import Fido, attrs as a #from sunpy.net.vso import attrs as avso
from sunpy.time import TimeRange

import astropy.units as u
from astropy.io import fits

import h5py

import csv

from skimage.transform import probabilistic_hough_line
from skimage.feature import canny, blob_log

"""
Checks that a fits file contains data otherwise returns unequal axes which results in such a file being discarded. 
"""
def readfits(filename):
    hdrs = []
    try:
        ft = fits.open(filename, memmap=False)
        ft.verify('fix')
        for i in range(len(ft)):
            hdrs.append(ft[i].header)
            hdr = ft[i].header
            data = ft[i].data
            #print(np.shape(data))
            
        axis1 = hdr['naxis1']
        axis2 = hdr['naxis2']
        axisnum = hdr['naxis']
        ft.close()
        
    except ValueError:
        axis1 = 1
        axis2 = 2
        axisnum = 3
        data = None
        hdrs = None
            
    return axis1, axis2, data, hdrs, axisnum

"""
Writes a fits file while preserving the header metadata
"""
def writefits(filename, data, hdrs, home_dir):
    if not os.path.exists(f'{home_dir}{filename}.fits'):
        hdrprimary = fits.PrimaryHDU(header=hdrs[len(hdrs)-1])
        fitsimg = fits.ImageHDU(data)
        hdul = fits.HDUList([hdrprimary, fitsimg])
        hdul.writeto(f'{home_dir}{filename}.fits', 'silentfix+ignore')
        hdul.close() #added this line

"""
Checks original image for missing pixel data.
"""
def holes(filename,base,mission):
    
    #filename = str(filename)
    ft = fits.open(filename, memmap=False)  
    ft.verify('fix')  
    
    if ('MDI' in base) or (mission == 'SDO'):
        hdr = ft[1].header
        data = ft[1].data
                
    else:
        hdr = ft[0].header
        data = ft[0].data
    
    ft.close()

    try:
        x_coord = hdr['CRPIX1']
        y_coord = hdr['CRPIX2']
        comments = hdr['COMMENT'][-1]
        missing_vals = int(hdr['MISSVALS'])
        blank_val = int(hdr['BLANK'])
        
        if 'N_MISSING_BLOCKS' in comments:
            missing_blocks = int(str(comments).split('=')[1])
    
    except KeyError:
        x_coord = hdr['naxis1'] / 2.
        y_coord = hdr['naxis2'] / 2.
        missing_blocks = 0
        missing_vals = 0
        blank_val = -100000

    y_ind,x_ind = np.indices((hdr['naxis1'],hdr['naxis2']))
    rsquared = (x_ind - x_coord)**2 + (y_ind - y_coord)**2
    
    matches_MDI_HMI_AIA = ['mdi', '96m', 'hmi', '720s', 'aia'] #'MDI' #in original downloaded files prior to renaming
    
    if ('efz' in filename): #good for all SOHO EIT and think that this also should work for SDO AIA products 
        rad = x_coord*np.sqrt(2)
        indices = np.where(rsquared.flatten() < rad**2)[0]
        zeros_ind = np.where(data.flatten()[indices] == 0.)[0]
        zeros_ind_len = len(zeros_ind)
        nan_ind = np.where(data.flatten()[indices] != data.flatten()[indices])[0]
        nan_ind_len = len(nan_ind)

        if (zeros_ind_len > 100) or (missing_blocks > 0) or (nan_ind_len > 100) or (missing_vals > 0):
            return True #so image not useable as there are holes
        else:
            return False #can use this image
    
    elif any([x in filename for x in matches_MDI_HMI_AIA]):
        rad1 = float(x_coord)
        rad2 = 0.6*float(x_coord)
        indices_rad1 = np.where(rsquared.flatten() < rad1**2)[0]
        indices_rad2 = np.where(rsquared.flatten() < rad2**2)[0]
        zeros_ind = np.where(data.flatten()[indices_rad1] == 0.)[0]
        nan_ind = np.where(data.flatten()[indices_rad2] != data.flatten()[indices_rad2])[0]
        blank_ind = np.where(data.flatten()[indices_rad2] == blank_val)[0]
        zeros_nan_ind_len = len(list(zeros_ind) + list(nan_ind))
        
        if (len(nan_ind) > 100) or (missing_vals > 0) or (len(blank_ind) > 100): #only nan_ind for calibrated JSOC MDI images which can have many zeros #worked for uncalibrated images: zeros_nan_ind_len > 100:
            return True #so image not useable as there are holes
        else:
            return False #can use this image

    elif 'LASCO_C3' in filename: 
        rad = 0.8*x_coord
        indices = np.where(rsquared.flatten() < rad**2)[0]
        zeros_ind = np.where(data.flatten()[indices] == 0.)[0]
        zeros_ind_len = len(zeros_ind)  
        nan_ind = np.where(data.flatten()[indices] != data.flatten()[indices])[0]
        nan_ind_len = len(nan_ind)
                
        if (zeros_ind_len > 100) or (nan_ind_len > 100): ###### For the case if holes are nans in lev1 images. Holes are 0 in lev0 images. ######
            return True #so image not useable as there are holes
        else:
            return False #can use this image   

    
    elif 'LASCO_C2' in filename:
        rad1 = 160 #this seems good
        #print('rad1:', rad1)
        rad2 = int(x_coord)
        indices = np.where((rad2**2 > rsquared.flatten()) & (rsquared.flatten() > rad1**2))[0]
        zeros_ind = np.where(data.flatten()[indices] == 0.)[0]
        zeros_ind_len = len(zeros_ind)
        nan_ind = np.where(data.flatten()[indices] != data.flatten()[indices])[0]
        nan_ind_len = len(nan_ind)
     
        if (zeros_ind_len > 100) or (nan_ind_len > 100): ###### For the case if holes are nans in lev1 images. Holes are 0 in lev0 images. ######
            return True #so image not useable as there are holes
        else:
            return False #can use this image
        
"""
Reduces the image dimensions via one of the following four methods: subsampling, interpolation, min pooling, or max pooling.
"""
def data_reducer(data,flag,image_size_output,axis1_shape):
    scale_factor = int(axis1_shape/image_size_output)
    
    if flag == 'subsample':
        reduced_data = data[::scale_factor].T[::scale_factor].T #subsampling image; every other row,column
    elif flag == 'interp': #linear interpolation with anti_aliasing and range preserving
        reduced_data = rescale(data, (1/scale_factor), order=1, anti_aliasing=True, preserve_range=True)
    elif flag == 'minpool': #min pooling each block
        reduced_data = block_reduce(data, block_size=(scale_factor,scale_factor), func=np.min)
    elif flag == 'maxpool': #max pooling each block
        reduced_data = block_reduce(data, block_size=(scale_factor,scale_factor), func=np.max)
    
    return reduced_data

"""
Resumes adding fits files if had previous interuption of program. Looks up times present and continues to add files at specified time interval.
"""
def prev_time_resumer(home_dir, base, time_range_orig, date_time_end, mission): 
#CAN RE-RUN PROGRAM FROM THE LAST DATE ON WHICH STOPPED; WILL PICK UP THE TIMES THAT ARE PRESENT AND CONTINUE! For both resuming on same day and next day.
### CHECKS WHETHER THE START DAY THAT ENTERED IS ALREADY CONTAINED IN THE FILES OF PREVIOUS DAY AND WILL SET START_DATE FROM THAT EXACT TIME! 
### ALSO WORKS IF START ON A NEW DAY AND ARE LOOKING BACK ON THE PREVIOUS DAY
    
    filepath = home_dir + base + f'_{mission}' + '/'

    data_files_pre_pre = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    data_files_pre = [f for f in data_files_pre_pre if 'fits' in f] 
    data_files = np.sort(data_files_pre)
    
    if len(data_files) != 0:
        prev_time_pre = data_files[-1]
        
        if ('EIT' in str(prev_time_pre)) or ('AIA' in str(prev_time_pre)): 
            prev_time = [prev_time_pre.split('_')[2]]
        else:
            prev_time = [prev_time_pre.split('_')[3]]
             
        time_orig_pre = str(time_range_orig.start)
        time_orig = ''.join(time_orig_pre.split(' ')[0].split('-'))
        
        if str(prev_time[0][0:8]) == time_orig:
            time_begin = prev_time[0]
            time_range = TimeRange(time_begin, date_time_end)
        else:
            time_range = time_range_orig            
    
    elif len(data_files) == 0:
        prev_time = []
        time_range = time_range_orig   
    
    return prev_time, time_range


"""
Assists with final time_start_name_new and time_finish_name_new names below
"""
def date_name_maker(date_name):

    date_name_chunks = [date_name[i:i+2] for i in range(0,len(date_name),2)]
    date_name_new = ''.join(date_name_chunks[0:2])+'-'+'-'.join(date_name_chunks[2:4])+'-'+':'.join(date_name_chunks[4:7])
    
    return date_name_new


"""
Specified date and time in the naming of both the h5py and csv files
"""
def data_name_selector(home_dir, base, date_start, date_finish, mission):

    filepath = home_dir + base + f'_{mission}' + '/'

    data_files_pre_pre = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    data_files_pre = [f for f in data_files_pre_pre if 'fits' in f] #to ensure that only FITS files are collected, just in case
    data_files = np.sort(data_files_pre)
    
    if len(data_files) != 0: 
        time_start_name_pre = data_files[0]
        time_finish_name_pre = data_files[-1]         
        
        if ('EIT' in str(time_start_name_pre)) or ('AIA' in str(time_start_name_pre)):
            time_start_name = str(time_start_name_pre.split('_')[2])
            time_finish_name = str(time_finish_name_pre.split('_')[2])        
        else:
            time_start_name = str(time_start_name_pre.split('_')[3])
            time_finish_name = str(time_finish_name_pre.split('_')[3])
        
        time_start_name_new = date_name_maker(time_start_name)
        time_finish_name_new = date_name_maker(time_finish_name)
        
    else:
        time_start_name_new = date_start 
        time_finish_name_new = date_finish   
        
    return time_start_name_new, time_finish_name_new

"""
Generated a compressed h5py data cube from all fits files present in a product folder
"""
def data_cuber(home_dir, base, date_start, date_finish, flag, time_window, image_size_output, lev1_LASCO, mission):

    filepath = home_dir + base + f'_{mission}' + '/'

    data_files_pre_pre = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    data_files_pre = [f for f in data_files_pre_pre if 'fits' in f]
    data_files = np.sort(data_files_pre) #to have chronological order and to sink order with list of individual product times

    data_content_list = []
    for elem in data_files:
        axdim1,axdim2,data_content,header_content,axisnum = readfits(f'{filepath}{elem}')
        if f'{mission}' in elem:
            data_content_list.append(data_content)

    if data_content_list:
        data_content_stack = np.stack(data_content_list)
    else:
        data_content_stack = []
                  
    time_start_name_new, time_finish_name_new = data_name_selector(home_dir, base, date_start, date_finish, mission)
        
    data_cube = h5py.File(f'{home_dir}{time_start_name_new}_to_{time_finish_name_new}_{base}_{flag}_{time_window}_LASCOlev1-{lev1_LASCO}_{mission}_{image_size_output}.h5', 'w')
    data_cube.create_dataset(f'{base}_{mission}_{image_size_output}', data=data_content_stack, compression="gzip")
    data_cube.close()
                            
    return data_cube

"""
Used SunPy's Fido to search for the user specified products
"""
def product_search(base,time_range,client,mission,time_window):

    #drms_product_matches = ['MDI', 'HMI', 'AIA']
    euv_matches = ['94', '131', '171', '193', '211', '304', '335'] 
    uv_matches = ['1600', '1700']
    
    if 'EIT' in base:
        wavelen = int(base[3:6])
        print('wavelen:', wavelen)
        #a.vso.Time(time_range,date_time_start) works with VSO but not with JSOC; so putting everything to form: a.vso.Time(time_range.start,time_range.end)
        product_results = Fido.search(a.vso.Time(time_range.start,time_range.end), a.vso.Source(f'{mission}'), a.vso.Instrument('EIT'), a.vso.Provider('SDAC'), a.vso.Wavelength(wavelen * a.vso.u.Angstrom, wavelen * a.vso.u.Angstrom))
        #print('product_results:', product_results)
    
    elif ('MDI' in base) or (mission == 'SDO'): #elif any([x in base for x in drms_product_matches]):
        
        ts = '_'.join(str(time_range.start).split(' '))+'_TAI'
        tf = '_'.join(str(time_range.end).split(' '))+'_TAI'
        
        if 'MDI' in base:
            product_results = client.export(f'mdi.fd_M_96m_lev182[{ts}-{tf}]')
        
        elif 'HMI' in base:
            product_results = client.export(f'hmi.M_720s[{ts}-{tf}]')        
        
        elif 'AIA' in base:
            wavelen = int(base[3:6])
            print('wavelen:', wavelen)
            
            if any([x in base for x in euv_matches]):
                product_results = client.export(f'aia.lev1_euv_12s[{ts}-{tf}@{time_window}h][{wavelen}]{{image}}')
            
            elif any([x in base for x in uv_matches]):
                product_results = client.export(f'aia.lev1_uv_24s[{ts}-{tf}@{time_window}h][{wavelen}]{{image}}')
            
            elif wavelen == 4500:
                product_results = client.export(f'aia.lev1_vis_1h[{ts}-{tf}][{wavelen}]{{image}}')                                                               
                
        ### CAN'T WORK WITH ACTUAL QUERY SINCE THE NUMBER OF ACTUAL RETURNED ITEMS BY EXPORT GETS RID OF THE GHOST FILES; 19990101-19990301; 886 FILES -> 314 FILES!!! 
        #######product_results_pre = client.query(f'mdi.fd_M_96m_lev182[{ts}-{tf}]', key = 'T_REC') #DRMS works with TAI times which do not make use of any leap seconds.
        
        ### THIS WAS FOR WHEN ACCESSING NON-GHOST FILES WITH FIDO. 
        #'mdi.fd_M_96m_lev182' is the calibration MDI series
        #product_results = Fido.search(a.jsoc.Time(time_range.start,time_range.end), a.jsoc.Series('mdi.fd_M_96m_lev182'))        
        
        ### THIS WAS FOR WHEN MDI PRODUCTS WERE AVAILABLE VIA SDAC
        #product_results = Fido.search(a.vso.Time(time_range.start,time_range.end), a.vso.Source('SOHO'), a.vso.Instrument('MDI'), a.vso.Provider('SDAC'), a.vso.Physobs('LOS_MAGNETIC_FIELD'))
    
    elif 'LASCO' in base: #LASCO is a SOHO instrument so hardcoded in a.vso.Source('SOHO') so as not to get 'SDO' in here!
        detector = base.split('_')[1]
        #a.vso.Time(time_range,date_time_start) works with VSO but not with JSOC; so putting everything to form: a.vso.Time(time_range.start,time_range.end)
        product_results = Fido.search(a.vso.Time(time_range.start,time_range.end), a.vso.Provider('SDAC'), a.vso.Source('SOHO'), a.vso.Instrument('LASCO'), a.vso.Detector(detector))
    
    return product_results, client

"""
Check JSOC DRMS for ghost files of 0 MB in query. These ghost files are returned by client.query() but client.export() will fail if only ghost files of 0 MB size are present.
THIS IS FROM THE OLD METHOD WHERE HAD QUERY OBJECT AND EXPORT OBJECT SEPERATE. SINCE QUERY OBJECT CONTAINS ENTRIES OF GHOST FILES AND THE EXPORT OBJECT DOES NOT, THE EXPORT OBJECT SHOULD BE USED TO 
EXTRACT ALL INFORMATION INCLUDING THE TIMES. 
"""

def ghost_file_check(product_results): #(client,time_range):
    
    #ts = '_'.join(str(time_range.start).split(' '))+'_TAI'
    #tf = '_'.join(str(time_range.end).split(' '))+'_TAI'    
    #client_export_MDI = client.export(f'mdi.fd_M_96m_lev182[{ts}-{tf}]')
    
    return product_results.has_failed() #client_export_MDI.has_failed()

"""
Returns the indices corresponding to the proper sizes of each of the data products to ensure that they are single, two dimensional images.
"""
def index_of_sizes(base,product_results,fits_headers,lev1_LASCO,mission):
    
    matches_EIT = ['171', '304', '284']
    
    if 'EIT195' in base:
        size_list = [elem['size'] for elem in product_results.get_response(0)[:]]
        ind_2059 = np.where(np.array(size_list) == 2059)[0] #this is for 1024 x 1024 pixels
        ind_523 = np.where(np.array(size_list) == 523)[0] ##this is for 512 x 512 pixels
        ind = np.sort(list(ind_2059) + list(ind_523)) #important to sort here since combining two lists!
        
    elif ('MDI' in base) or (mission == 'SDO'):
        calib_file_num = product_results.data.count()['record'] #product_results.count()['T_REC'] when had a DRMS query object instead #product_results.file_num works for Fido queries
        
        if calib_file_num == 0:
            ind = []
        else:
            ind = np.arange(calib_file_num)
            
            #size_list = [elem['size'] for elem in product_results.get_response(0)[:]]
            #print(np.unique(size_list), len(size_list))
            #ind = np.where(np.array(size_list) == 4115.0)[0] ### old method when MDI uncalibrated was served on VSO SDAC
            #print(len(ind))        
        
    elif 'LASCO' in base:
        
        if (lev1_LASCO == 'Y') or (lev1_LASCO == 'y'):
            size_list = [elem['size'] for elem in product_results.get_response(0)[:]]
            ind = np.where(np.array(size_list) == 4106.0)[0] #4106 is for calibrated level-1.0 #2100.0 was for uncalibrated level-0.5 
        
        elif (lev1_LASCO == 'N') or (lev1_LASCO == 'n'):
            size_list = [int(np.ceil(elem['size'] / 100.0))*100 for elem in product_results.get_response(0)[:]] 
            ind = np.where(np.array(size_list) == 2100.0)[0] #4106 is for calibrated level-1.0 #2100.0 was for uncalibrated level-0.5         
        
        #[int(np.ceil(elem['size'] / 100.0))*100 for elem in product_results.get_response(0)[:]] #this was for uncalibrated level-0.5
        
    elif (any([x in base for x in matches_EIT])) and (mission == 'SOHO'): #'SOHO' here because have overlap in wavelength between missions at 171 Angstrom.
        size_list = [elem['size'] for elem in product_results.get_response(0)[:]]
        #print(np.unique(size_list), len(size_list))
        ind = np.where(np.array(size_list) == 2059)[0]        
        #print(len(ind))
    
    fits_headers = fits_headers    
    
    return ind, fits_headers
   
"""
From the appropriate file size image indices, this function picks up those times and their corresponding indices that are seperated by the user defined time window. 
The indices returned by this function will be used in the product object returned by Fido search.
"""
def fetch_indices(base,ind,product_results,time_window,look_ahead, prev_time, mission):
    
    all_size_sieved_times_pre = [] #local list to populate at each loop
    all_time_window_sieved_times_product_times = []  #local list to populate at each loop
    
    if ('MDI' in base) or (mission == 'SDO'):
        
        drms_export_list = product_results.data.record  #product_results.data.T.get_values()[0] -> this works for drms 0.5.5 on my laptop but not for drms 0.5.7 on Juniper
        all_size_sieved_times_pre = [drms.to_datetime(elem.split('[')[1].split(']')[0]).strftime('%Y%m%d%H%M%S') for i,elem in enumerate(drms_export_list)] #only way to get to times from export object
        
        #THIS IS WHEN HAD A QUERY OBJECT: DRMS_MDI_UTC_time_stamps = drms.to_datetime(product_results.T_REC) 
        ### converting TAI DRMS times to UTC times. TAI times do not contain the leap second in contrast to UTC times.
        #all_size_sieved_times_pre = [elem.strftime('%Y%m%d%H%M%S') for i,elem in enumerate(DRMS_MDI_UTC_time_stamps)] #can bypass use of the 'ind' list on this line since not searching by file size in JSOC.
        
        ### OLDER METHOD WHEN FIDO WAS BEING USED TO OBTAIN MDI IMAGES FROM JSOC
        #for value in ind:
            #TAI_entity_str =  product_results.get_response(0)[int(value)]['T_REC']
            #TAI_entity_time = ''.join(str(TAI_entity_str).split('_TAI')[0].split('_')[0].split('.')[0:3] + str(TAI_entity_str).split('_TAI')[0].split('_')[1].split(':')[0:3])     
            #all_size_sieved_times_pre.append(TAI_entity_time)
            #print('all_size_sieved_times_pre:', all_size_sieved_times_pre)
                    
    else:
        
        for value in ind:
            all_size_sieved_times_pre.append(product_results.get_response(0)[int(value)]['time']['start'])
            
    all_size_sieved_times = list(np.unique(all_size_sieved_times_pre))
    all_size_sieved_times_aug = prev_time + all_size_sieved_times #prev_time = [] for the very first loop and [last best time from previous loop] for subsequent loops.

    if all_size_sieved_times_aug:
        if prev_time:
            local_time_range = TimeRange(all_size_sieved_times_aug[0],timedelta(hours=time_window)).next() #next() is the important difference here.
        else:
            local_time_range = TimeRange(all_size_sieved_times_aug[0],timedelta(hours=time_window))       
        for i,time_value in enumerate(all_size_sieved_times_aug):
            if time_value in local_time_range:
                all_time_window_sieved_times_product_times.append(time_value)
                local_time_range = TimeRange(time_value,timedelta(hours=time_window)).next() #important distinction between this local_time_range and the intializing one is the presence of time_value          
            elif parser.parse(time_value) > local_time_range.end: 
                all_time_window_sieved_times_product_times.append(time_value)
                local_time_range = TimeRange(time_value,timedelta(hours=time_window)).next()
            else:
                continue        
    new_inds = [np.where(np.array(all_size_sieved_times_pre) == entry)[0][0] for entry in all_time_window_sieved_times_product_times]
    
    fetch_indices_product = ind[new_inds]
    
    return all_size_sieved_times_pre, fetch_indices_product
    
"""
Using Fido's returned query object that has now been sieved for proper data size in the case of EIT and LASCO products and user specified time window, the fileid is extracted from the accompanying Fido dictionary and wget is used to retrieve the product. Fido.fetch() method is used to obtain calibrated MDI products.
"""
def product_retriever(base,product_results,indiv_ind,url_prefix,home_dir,email,all_size_sieved_times_pre, fits_headers, client, mission): #added email, all_size_sieved_times_pre, fits_headers
    
    if ('MDI' in base) or (mission == 'SDO'): 
    
        out_dir = f'{home_dir}{base}_{mission}/'
        
        time_range_for_DRMS = TimeRange(all_size_sieved_times_pre[indiv_ind], timedelta(minutes = 1)) #minutes hardcoded to 1 since just fetching one JSOC image at a time and want small time around it.
        ts_for_DRMS = '_'.join(str(time_range_for_DRMS.start).split(' '))+'_TAI'
        tf_for_DRMS = '_'.join(str(time_range_for_DRMS.end).split(' '))+'_TAI'
        
        euv_matches = ['94', '131', '171', '193', '211', '304', '335'] 
        uv_matches = ['1600', '1700']
                
        if (fits_headers == 'Y') or (fits_headers == 'y'):
            
            if 'MDI' in base:
                client_export_drms = client.export(f'mdi.fd_M_96m_lev182[{ts_for_DRMS}-{tf_for_DRMS}]', method='url', protocol='fits')
            
            elif 'HMI' in base:
                client_export_drms = client.export(f'hmi.M_720s[{ts_for_DRMS}-{tf_for_DRMS}]', method='url', protocol='fits')        
            
            elif 'AIA' in base:
                wavelen = int(base[3:6])
                
                if any([x in base for x in euv_matches]):
                    client_export_drms = client.export(f'aia.lev1_euv_12s[{ts_for_DRMS}-{tf_for_DRMS}][{wavelen}]{{image}}', method='url', protocol='fits')
            
                elif any([x in base for x in uv_matches]):
                    client_export_drms = client.export(f'aia.lev1_uv_24s[{ts_for_DRMS}-{tf_for_DRMS}][{wavelen}]{{image}}', method='url', protocol='fits')           
                
                elif '4500' in base:
                    client_export_drms = client.export(f'aia.lev1_vis_1h[{ts_for_DRMS}-{tf_for_DRMS}][{wavelen}]{{image}}', method='url', protocol='fits')                
        
        elif (fits_headers == 'N') or (fits_headers == 'n'):
            
            if 'MDI' in base:            
                client_export_drms = client.export(f'mdi.fd_M_96m_lev182[{ts_for_DRMS}-{tf_for_DRMS}]')
                                    
            elif 'HMI' in base:
                client_export_drms = client.export(f'hmi.M_720s[{ts_for_DRMS}-{tf_for_DRMS}]')                    
            
            elif 'AIA' in base:
                wavelen = int(base[3:6])
                
                if any([x in base for x in euv_matches]):
                    client_export_drms = client.export(f'aia.lev1_euv_12s[{ts_for_DRMS}-{tf_for_DRMS}][{wavelen}]{{image}}')
            
                elif any([x in base for x in uv_matches]):
                    client_export_drms = client.export(f'aia.lev1_uv_24s[{ts_for_DRMS}-{tf_for_DRMS}][{wavelen}]{{image}}')           
                
                elif '4500' in base:
                    client_export_drms = client.export(f'aia.lev1_vis_1h[{ts_for_DRMS}-{tf_for_DRMS}][{wavelen}]{{image}}') 

                        
            #if ghost_MDI_file_check(client,time_range_for_DRMS) == False:             
                #client_export_drms = client.export(f'mdi.fd_M_96m_lev182[{ts_for_DRMS}-{tf_for_DRMS}]')                    
                #query_result_pre = client_export_drms.download(out_dir,0) #always the first elemement since just downloading one MDI image at a time
                #query_result = list(query_result_pre.download)        
            #else:
                #query_result = []

        if client_export_drms.status == 0:
            query_result_pre = client_export_drms.download(out_dir,0) #always the first elemement since just downloading one MDI image at a time
            query_result = list(query_result_pre.download)
        
        elif client_export_drms.status != 0:
            query_result_pre = []
            query_result = []
                            
        
        if (all_size_sieved_times_pre != []) and (query_result == []):
            print('sleep for 15 minutes and then retry DRMS download')
            time.sleep(900) 
            query_result_pre = client_export_drms.download(out_dir,0)
            query_result = list(query_result_pre.download)
            
                    
        ### OLDER METHOD WHEN FIDO WAS BEING USED TO OBTAIN MDI IMAGES FROM JSOC
        #time_range_fido = TimeRange(all_size_sieved_times_pre[indiv_ind], timedelta(minutes = 5)) #minutes hardcoded arbitrarily to 5 since just fetching one MDI at a time and want small time around it.
        #product_results_fido = Fido.search(a.jsoc.Time(time_range_fido.start, time_range_fido.end), a.jsoc.Series('mdi.fd_M_96m_lev182'),a.jsoc.Notify(email))
        #query_result_pre = [Fido.fetch(product_results_fido, path = f'{home_dir}{base}/')]
        #query_result = query_result_pre[0]
        #print('query_result insdie product_retriever MDI:', query_result)
        
        #if (all_size_sieved_times_pre != []) and (query_result == []):
            #print('sleep for 15 minutes and then retry Fido fetch command')
            #time.sleep(900) 
            #query_result_pre = [Fido.fetch(product_results_fido, path = f'{home_dir}{base}/')]
            #query_result = query_result_pre[0]
    
    elif ('MDI' not in base) and (mission == 'SOHO'):
        fileid = product_results.get_response(0)[int(indiv_ind)]['fileid']
        item_wget =  url_prefix + fileid
        cmd = 'wget' + ' ' + '--retry-connrefused' + ' ' + '--waitretry=1' + ' ' + '--read-timeout=20' + ' ' + '--timeout=15' + ' ' + '-t' + ' ' + '0' + ' ' + '--continue' + ' ' + item_wget + ' ' + '-P' + ' ' + f'{home_dir}{base}_{mission}'     
        args = shlex.split(cmd)    
    
        try: 
            wget_output = subprocess.check_output(args, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as err:
            print('Error output:\n', err.output) 
            print('sleep for 15 minutes and then retry command')
            time.sleep(900)
            wget_output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    
        downloaded_fileid = fileid.split('/')[-1]
        query_result = [f'{home_dir}{base}_{mission}/{downloaded_fileid}']
    
    return query_result
    
"""
Planet and comet transients filter. Uses Probabilistic Hough Transform on Canny edge processed images to detect straight lines. 
The lines object returned by probabilistic_hough_line() is a list.
Lines == False means can use image.
"""
def planet_comet_transient_filter(base, data_product): #lev1_LASCO, #only C2, C3, provide URLS: for C3 lev 1: use theta explicitely for only the horizontal angles since pile on is very well apparent!!! 

    detector = base.split('_')[1]
    data_product_log10 = np.log10(data_product)
    edges_data_product_log10 = canny(data_product_log10,1)
    
    if detector == 'C3':
        #if (lev1_LASCO == 'Y') or (lev1_LASCO == 'y'):
        lines = probabilistic_hough_line(edges_data_product_log10, threshold=80, line_length=5, line_gap=0, seed=1, theta=np.array([-np.pi/2, np.pi/2]))
    
        #elif (lev1_LASCO == 'N') or (lev1_LASCO == 'n'):
            #lines_pre = probabilistic_hough_line(edges_planet, threshold=80, line_length=5, line_gap=0, seed=1)

        #if len(lines) >= 3: #so as not to pick up lines from the pile-on or occulter
                #lines = True
                    
    elif detector == 'C2':
        lines = probabilistic_hough_line(edges_data_product_log10, threshold=80, line_length=5, line_gap=0, seed=1, theta=np.array([-np.pi/2, np.pi/2]))
    
    if len(lines) > 5:
        return True #can't use this image as straight lines detected
    
    elif not lines:
        return False #can use this image as no straight lines detected



"""
Cosmic ray filter on downsized LASCO C2 images; doesn't capture cosmic rays on C3. Uses Laplacian of Gaussian to detect blobs. False means can use image.
"""
def cosmic_ray_filter(base, reduced_product_data): #only works for C2 #i think don't need to exclusively distinguish between C2 and C3 here

    blobs_log = blob_log(reduced_product_data, max_sigma=30, num_sigma=10, threshold=.1)

    ''' ### BLOB DETECTOR SHOULD BE APPLIED LATER ON THE COLLECTED IMAGES SIEVED BY THE FILTERS HOLES AND TRANSIENTS; BLOB FILTER NOT WORKING LIKE THIS! Has 23 FP and 0 TP.
    if len(blobs_log) > 30:
        return True #can't use this image as blobs detected
    
    else:
        return False #can use this image as no blobs detected

    '''
    return False        

"""
Checks the downloaded fits files for holes and discards them if holes are found. Repeats procedure at each time as long as an image contained missing pixels. 
"""
def product_distiller(fetch_indices_product_orig, base, all_size_sieved_times_pre, ind, product_results, look_ahead, time_window, url_prefix, flag, image_size_output, home_dir, email, fits_headers, lev1_LASCO, client, mission):

    holes_product_list = []
    blobs_product_list = []
    transients_product_list = []
    
    unreadable_file_ids_product_list = []

    all_time_window_sieved_times_product_times = [] 
    all_time_window_sieved_times_product_times_inds_list = [] 

    fetch_indices_product = fetch_indices_product_orig.copy()
    
    for i,elem in enumerate(fetch_indices_product): #the i index retains the original number of members of the original fetch_indices_product. The fetch_indices_product list is modified when holes occur.
        if (i > len(fetch_indices_product)-1): #fetch_indices_product is modified in the program to account for times corresponding to holes in images. When all fitting times exhausted then break out of loop.
            break        
        indiv_ind = fetch_indices_product[i] # equivalently int(elem)
        query_result = product_retriever(base,product_results,indiv_ind,url_prefix,home_dir,email,all_size_sieved_times_pre,fits_headers,client,mission)
        axis1_product, axis2_product, data_product, header_product, axisnum_product = readfits(query_result[0])
            
        if (data_product is not None) and (axis1_product == axis2_product) and (axisnum_product == 2):

            if not holes(query_result[0],base,mission): #so if not True; so no holes; can use image
                reduced_product_data = data_reducer(data_product,flag,image_size_output,axis1_product)
                
                if ('MDI' in base) or (mission == 'SDO'):
                    time_data = all_size_sieved_times_pre[indiv_ind]                                        
                    writefits(f'{base}_{mission}/{mission}_{base}_{time_data}_{image_size_output}', reduced_product_data, header_product, home_dir)
                    os.remove(query_result[0]) #delete original downloaded file
                    all_time_window_sieved_times_product_times.append(time_data)
                    all_time_window_sieved_times_product_times_inds_list.append(indiv_ind)  

                elif 'EIT' in base:
                    time_data = product_results.get_response(0)[int(indiv_ind)]['time']['start']
                    writefits(f'{base}_{mission}/{mission}_{base}_{time_data}_{image_size_output}', reduced_product_data, header_product, home_dir)
                    os.remove(query_result[0]) #delete original downloaded file
                    all_time_window_sieved_times_product_times.append(time_data)
                    all_time_window_sieved_times_product_times_inds_list.append(indiv_ind)
                
                elif 'LASCO' in base:
                    time_data = product_results.get_response(0)[int(indiv_ind)]['time']['start']                    
                    if (not planet_comet_transient_filter(base, data_product)) and (not cosmic_ray_filter(base, reduced_product_data)): #if both line list and blob lost is empty then can use LASCO image.

                        writefits(f'{base}_{mission}/{mission}_{base}_{time_data}_{image_size_output}', reduced_product_data, header_product, home_dir)
                        os.remove(query_result[0]) #delete original downloaded file
                        all_time_window_sieved_times_product_times.append(time_data)
                        all_time_window_sieved_times_product_times_inds_list.append(indiv_ind)   
                        
                    elif (planet_comet_transient_filter(base, data_product)) or (cosmic_ray_filter(base, reduced_product_data)): #if either line list or blob list is not emplty then can't use LASCO image
                        
                        if not cosmic_ray_filter(base, reduced_product_data): #so only planet_comet transition
                            transient_loc = url_prefix + product_results.get_response(0)[int(indiv_ind)]['fileid'] 
                            transients_product_list.append((transient_loc, str(time_data)))
                            #transients_time_val = str(time_data)
                            #ind_timespickup = np.where(np.array(all_size_sieved_times_pre) == transients_time_val)[0][0]
                            #zoomed_time_range = TimeRange(transients_time_val,timedelta(hours=time_window))                        
                        
                        elif not planet_comet_transient_filter(base, data_product):
                            blob_loc = url_prefix + product_results.get_response(0)[int(indiv_ind)]['fileid'] 
                            blobs_product_list.append((blob_loc, str(time_data)))
                            #blobs_time_val = str(time_data)
                            #ind_timespickup = np.where(np.array(all_size_sieved_times_pre) == blobs_time_val)[0][0]
                            #zoomed_time_range = TimeRange(blobs_time_val,timedelta(hours=time_window))                          
                        
                        zoomed_time_range = TimeRange(time_data,timedelta(hours=time_window))                                                  
                        ind_timespickup = np.where(np.array(all_size_sieved_times_pre) == str(time_data))[0][0]                        
                        os.remove(query_result[0]) #delete original downloaded file
                        
                        fetch_inds_to_try_list = [] 
                        for time_val in all_size_sieved_times_pre[ind_timespickup+1: ind_timespickup + look_ahead]:
                            if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                                ind_next_good_time = np.where(np.array(all_size_sieved_times_pre) == time_val)[0][0]
                                fetch_indices_next_good = ind[ind_next_good_time]
                                fetch_inds_to_try_list.append(fetch_indices_next_good)

                        for index in fetch_inds_to_try_list:
                            query_result_next = product_retriever(base,product_results,index,url_prefix,home_dir,email,all_size_sieved_times_pre,fits_headers,client,mission)
                            axis1_next_good,axis2_next_good,data_next_good,header_next_good,axisnum_next_good = readfits(query_result_next[0])

                            if (data_next_good is not None) and (axis1_next_good == axis2_next_good) and (axisnum_next_good == 2):
                                reduced_product_data = data_reducer(data_next_good,flag,image_size_output,axis1_next_good)                                    
                                
                                if (not holes(query_result_next[0],base,mission)) and (not planet_comet_transient_filter(base, data_next_good)) and (not cosmic_ray_filter(base, reduced_product_data)):
                                    
                                    time_data = product_results.get_response(0)[int(index)]['time']['start']
                                    writefits(f'{base}_{mission}/{mission}_{base}_{time_data}_{image_size_output}', reduced_product_data, header_next_good, home_dir)
                                    os.remove(query_result_next[0]) #delete original downloaded file                                    
                                    all_time_window_sieved_times_product_times.append(time_data) #(time_val) #unsorted time location
                                    all_time_window_sieved_times_product_times_inds_list.append(index)
                                                                
                                    indiv_ind_modified_list = []
                                    localized_time_range = TimeRange(str(time_data),timedelta(hours=time_window)).next()
                                    for tval in all_size_sieved_times_pre:
                                        if parser.parse(tval) < localized_time_range.start:
                                            continue 
                                        elif tval in localized_time_range:
                                            ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                            indiv_ind_modified_new = ind[ind_time_new]
                                            indiv_ind_modified_list.append(indiv_ind_modified_new)
                                            localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()                                
                                        else:
                                            ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                            indiv_ind_modified_new = ind[ind_time_new]                                    
                                            next_orig_index = np.where(np.array(fetch_indices_product_orig) == np.array(indiv_ind_modified_new))[0]
                                            if len(next_orig_index) != 0:
                                                indiv_ind_modified_new = fetch_indices_product_orig[next_orig_index[0]]
                                                ind_next_index = np.where(np.array(ind) == indiv_ind_modified_new)[0][0]
                                                tval = all_size_sieved_times_pre[ind_next_index]
                                                indiv_ind_modified_list.append(indiv_ind_modified_new)
                                                localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()
                                            else:
                                                ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                                indiv_ind_modified_new = ind[ind_time_new]
                                                indiv_ind_modified_list.append(indiv_ind_modified_new)
                                                localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()
                                                                        
                                    if indiv_ind_modified_list:
                                        fetch_indices_product = list(np.zeros(i+1)) + list(indiv_ind_modified_list)                            
                                    else:
                                        fetch_indices_product = list(np.zeros(i+1))                            
                                    break

                                elif (holes(query_result_next[0],base,mission)) or (planet_comet_transient_filter(base, data_next_good)) or (cosmic_ray_filter(base, reduced_product_data)): #so if True, if there are holes
                                    
                                    time_data = product_results.get_response(0)[int(index)]['time']['start']                            
                                    
                                    if (not planet_comet_transient_filter(base, data_next_good)) and (not cosmic_ray_filter(base, reduced_product_data)): #so no transient
                                        hole_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                                        holes_product_list.append((hole_loc, str(time_data)))
                                    elif planet_comet_transient_filter(base, data_next_good):
                                        transient_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                                        transients_product_list.append((transient_loc, str(time_data)))                                    
                                    elif cosmic_ray_filter(base, reduced_product_data):
                                        blob_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                                        blobs_product_list.append((blob_loc, str(time_data)))                                      
                                    
                                    os.remove(query_result_next[0])
                                    continue 

                            elif (data_next_good is None) or (axis1_next_good != axis2_next_good) or (axisnum_next_good != 2):
                                    
                                unreadable_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                                unreadable_file_ids_product_list.append(unreadable_loc)
                                
                                os.remove(query_result_next[0])
                                continue
                
            elif holes(query_result[0],base,mission): #so if True, if there are holes
                if ('MDI' in base) or (mission == 'SDO'):
                    time_data = all_size_sieved_times_pre[indiv_ind]                                                        
                    hole_loc = product_results.data.record[indiv_ind].split('{')[0] 
                    #hole_loc = product_results.T_REC[indiv_ind] #product_results.get_response(0)[int(indiv_ind)]['T_REC'] #old Fido method to obtain TAI time from JSOC                                      
                
                else:
                    time_data = product_results.get_response(0)[int(indiv_ind)]['time']['start'] 
                    hole_loc = url_prefix + product_results.get_response(0)[int(indiv_ind)]['fileid']                       
                
                holes_product_list.append((hole_loc, str(time_data)))
                hole_time_val = str(time_data)
                os.remove(query_result[0]) #delete original downloaded file
                ind_timespickup = np.where(np.array(all_size_sieved_times_pre) == hole_time_val)[0][0]
                zoomed_time_range = TimeRange(str(hole_time_val),timedelta(hours=time_window))
                
                fetch_inds_to_try_list = [] 
                #the zeroth entry didn't have it so that's why plus 1 in the brackets
                for time_val in all_size_sieved_times_pre[ind_timespickup+1: ind_timespickup + look_ahead]:
                    if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                        ind_next_good_time = np.where(np.array(all_size_sieved_times_pre) == time_val)[0][0]
                        fetch_indices_next_good = ind[ind_next_good_time]
                        fetch_inds_to_try_list.append(fetch_indices_next_good)

                for index in fetch_inds_to_try_list:
                    query_result_next = product_retriever(base,product_results,index,url_prefix,home_dir,email,all_size_sieved_times_pre,fits_headers,client,mission)
                    axis1_next_good,axis2_next_good,data_next_good,header_next_good,axisnum_next_good = readfits(query_result_next[0])

                    if (data_next_good is not None) and (axis1_next_good == axis2_next_good) and (axisnum_next_good == 2):
                        
                        reduced_product_data = data_reducer(data_next_good,flag,image_size_output,axis1_next_good)
                        
                        if not holes(query_result_next[0],base,mission): #and (not planet_comet_transient_filter(base, data_product)) and (not cosmic_ray_filter(base, reduced_product_data)):
                            
                            if ('MDI' in base) or (mission == 'SDO'):
                                time_data = all_size_sieved_times_pre[index]                            
                                writefits(f'{base}_{mission}/{mission}_{base}_{time_data}_{image_size_output}', reduced_product_data, header_next_good, home_dir)
                                all_time_window_sieved_times_product_times.append(time_data) #(time_val) #unsorted time location
                                all_time_window_sieved_times_product_times_inds_list.append(index)
                                os.remove(query_result_next[0]) #delete original downloaded file
                            
                                indiv_ind_modified_list = []
                                localized_time_range = TimeRange(str(time_data),timedelta(hours=time_window)).next()
                                for tval in all_size_sieved_times_pre:
                                    if parser.parse(tval) < localized_time_range.start:
                                        continue 
                                    elif tval in localized_time_range:
                                        ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                        indiv_ind_modified_new = ind[ind_time_new]
                                        indiv_ind_modified_list.append(indiv_ind_modified_new)
                                        localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()                                
                                    else:
                                        ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                        indiv_ind_modified_new = ind[ind_time_new]                                    
                                        next_orig_index = np.where(np.array(fetch_indices_product_orig) == np.array(indiv_ind_modified_new))[0]
                                        if len(next_orig_index) != 0:
                                            indiv_ind_modified_new = fetch_indices_product_orig[next_orig_index[0]]
                                            ind_next_index = np.where(np.array(ind) == indiv_ind_modified_new)[0][0]
                                            tval = all_size_sieved_times_pre[ind_next_index]
                                            indiv_ind_modified_list.append(indiv_ind_modified_new)
                                            localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()
                                        else:
                                            ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                            indiv_ind_modified_new = ind[ind_time_new]
                                            indiv_ind_modified_list.append(indiv_ind_modified_new)
                                            localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()
                            
                                if indiv_ind_modified_list:
                                    fetch_indices_product = list(np.zeros(i+1)) + list(indiv_ind_modified_list)                            
                                else:
                                    fetch_indices_product = list(np.zeros(i+1))                            
                                break                            
                            
                            elif 'EIT' in base:
                                time_data = product_results.get_response(0)[int(index)]['time']['start']
                                writefits(f'{base}_{mission}/{mission}_{base}_{time_data}_{image_size_output}', reduced_product_data, header_next_good, home_dir)
                                all_time_window_sieved_times_product_times.append(time_data) #(time_val) #unsorted time location
                                all_time_window_sieved_times_product_times_inds_list.append(index)
                                os.remove(query_result_next[0]) #delete original downloaded file
                            
                                indiv_ind_modified_list = []
                                localized_time_range = TimeRange(str(time_data),timedelta(hours=time_window)).next()
                                for tval in all_size_sieved_times_pre:
                                    if parser.parse(tval) < localized_time_range.start:
                                        continue 
                                    elif tval in localized_time_range:
                                        ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                        indiv_ind_modified_new = ind[ind_time_new]
                                        indiv_ind_modified_list.append(indiv_ind_modified_new)
                                        localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()                                
                                    else:
                                        ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                        indiv_ind_modified_new = ind[ind_time_new]                                    
                                        next_orig_index = np.where(np.array(fetch_indices_product_orig) == np.array(indiv_ind_modified_new))[0]
                                        if len(next_orig_index) != 0:
                                            indiv_ind_modified_new = fetch_indices_product_orig[next_orig_index[0]]
                                            ind_next_index = np.where(np.array(ind) == indiv_ind_modified_new)[0][0]
                                            tval = all_size_sieved_times_pre[ind_next_index]
                                            indiv_ind_modified_list.append(indiv_ind_modified_new)
                                            localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()
                                        else:
                                            ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                            indiv_ind_modified_new = ind[ind_time_new]
                                            indiv_ind_modified_list.append(indiv_ind_modified_new)
                                            localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()                                
                                
                                if indiv_ind_modified_list:
                                    fetch_indices_product = list(np.zeros(i+1)) + list(indiv_ind_modified_list)                            
                                else:
                                    fetch_indices_product = list(np.zeros(i+1))                            
                                break                             
                            
                            elif 'LASCO' in base:    
                                time_data = product_results.get_response(0)[int(index)]['time']['start']                            
                                if (not planet_comet_transient_filter(base, data_next_good)) and (not cosmic_ray_filter(base, reduced_product_data)):
                                    writefits(f'{base}_{mission}/{mission}_{base}_{time_data}_{image_size_output}', reduced_product_data, header_next_good, home_dir)
                                    all_time_window_sieved_times_product_times.append(time_data) #(time_val) #unsorted time location
                                    all_time_window_sieved_times_product_times_inds_list.append(index)
                                    os.remove(query_result_next[0]) #delete original downloaded file
                            
                                    indiv_ind_modified_list = []
                                    localized_time_range = TimeRange(str(time_data),timedelta(hours=time_window)).next()
                                    for tval in all_size_sieved_times_pre:
                                        if parser.parse(tval) < localized_time_range.start:
                                            continue 
                                        elif tval in localized_time_range:
                                            ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                            indiv_ind_modified_new = ind[ind_time_new]
                                            indiv_ind_modified_list.append(indiv_ind_modified_new)
                                            localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()                                
                                        else:
                                            ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                            indiv_ind_modified_new = ind[ind_time_new]                                    
                                            next_orig_index = np.where(np.array(fetch_indices_product_orig) == np.array(indiv_ind_modified_new))[0]
                                            if len(next_orig_index) != 0:
                                                indiv_ind_modified_new = fetch_indices_product_orig[next_orig_index[0]]
                                                ind_next_index = np.where(np.array(ind) == indiv_ind_modified_new)[0][0]
                                                tval = all_size_sieved_times_pre[ind_next_index]
                                                indiv_ind_modified_list.append(indiv_ind_modified_new)
                                                localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()
                                            else:
                                                ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                                indiv_ind_modified_new = ind[ind_time_new]
                                                indiv_ind_modified_list.append(indiv_ind_modified_new)
                                                localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()                                
                                
                                    if indiv_ind_modified_list:
                                        fetch_indices_product = list(np.zeros(i+1)) + list(indiv_ind_modified_list)                            
                                    else:
                                        fetch_indices_product = list(np.zeros(i+1))                            
                                    break                                                         
                                
                                
                                elif (planet_comet_transient_filter(base, data_product)) or (cosmic_ray_filter(base, reduced_product_data)):                                
                                    
                                    if not planet_comet_transient_filter(base, data_product):
                                        blob_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                                        blobs_product_list.append((blob_loc, str(time_data)))
                                    elif not cosmic_ray_filter(base, reduced_product_data):
                                        transient_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid'] 
                                        transients_product_list.append((transient_loc, str(time_data)))
                                    
                                    os.remove(query_result_next[0]) #these were indented one outside on 27.03.21
                                    continue    #these were indented one outside on 27.03.21                      
                                
                        elif holes(query_result_next[0],base,mission): #so if True, if there are holes
                            
                            if ('MDI' in base) or (mission == 'SDO'):
                                time_data = all_size_sieved_times_pre[index]
                                hole_loc = product_results.data.record[index].split('{')[0]
                                #hole_loc = product_results.T_REC[index] #product_results.get_response(0)[int(index)]['T_REC']                            
                            
                            else:
                                time_data = product_results.get_response(0)[int(index)]['time']['start']
                                hole_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                            
                            holes_product_list.append((hole_loc, str(time_data)))
                            
                            os.remove(query_result_next[0])
                            continue 

                    elif (data_next_good is None) or (axis1_next_good != axis2_next_good) or (axisnum_next_good != 2):
                        
                        if ('MDI' in base) or (mission == 'SDO'):
                            unreadable_file_ids_product_list.append(product_results.data.record[index].split('{')[0]) 
                            #unreadable_file_ids_product_list.append(product_results.T_REC[index]) #product_results.get_response(0)[int(index)]['T_REC'])                        
                        
                        else:
                            unreadable_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                            unreadable_file_ids_product_list.append(unreadable_loc)
                        
                        os.remove(query_result_next[0])
                        continue


        elif (data_product is None) or (axis1_product != axis2_product) or (axisnum_product != 2):
            
            if ('MDI' in base) or (mission == 'SDO'):
                unreadable_file_ids_product_list.append(product_results.data.record[indiv_ind].split('{')[0])               
                #unreadable_file_ids_product_list.append(product_results.T_REC[indiv_ind]) #product_results.get_response(0)[int(indiv_ind)]['T_REC'])
                bad_time_val = all_size_sieved_times_pre[indiv_ind]       
            
            else:
                unreadable_loc = url_prefix + product_results.get_response(0)[int(indiv_ind)]['fileid']
                unreadable_file_ids_product_list.append(unreadable_loc)
                bad_time_val = product_results.get_response(0)[int(indiv_ind)]['time']['start']
            
            os.remove(query_result[0]) #delete original downloaded file
            ind_timespickup = np.where(np.array(all_size_sieved_times_pre) == bad_time_val)[0][0]
            zoomed_time_range = TimeRange(str(bad_time_val),timedelta(hours=time_window))

            fetch_inds_to_try_list = [] #gets reset for each new item
            for time_val in all_size_sieved_times_pre[ind_timespickup+1: ind_timespickup + look_ahead]:
                if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                    ind_next_good_time = np.where(np.array(all_size_sieved_times_pre) == time_val)[0][0]
                    fetch_indices_next_good = ind[ind_next_good_time]
                    fetch_inds_to_try_list.append(fetch_indices_next_good)

            for index in fetch_inds_to_try_list:
                query_result_next = product_retriever(base,product_results,index,url_prefix,home_dir,email,all_size_sieved_times_pre,fits_headers,client,mission)
                axis1_next_good,axis2_next_good,data_next_good,header_next_good,axisnum_next_good = readfits(query_result_next[0])

                if (data_next_good is not None) and (axis1_next_good == axis2_next_good) and (axisnum_next_good == 2):

                    if not holes(query_result_next[0],base,mission): #so if not True; so no holes; can use image
                        reduced_product_data = data_reducer(data_next_good,flag,image_size_output,axis1_next_good)
                        
                        if ('MDI' in base) or (mission == 'SDO'):
                            time_data = all_size_sieved_times_pre[index]                      
                            writefits(f'{base}_{mission}/{mission}_{base}_{time_data}_{image_size_output}', reduced_product_data, header_next_good, home_dir)
                            all_time_window_sieved_times_product_times.append(time_data) #(time_val) #unsorted time location
                            all_time_window_sieved_times_product_times_inds_list.append(index)
                            os.remove(query_result_next[0])                        
                        
                            indiv_ind_modified_list = []
                            localized_time_range = TimeRange(str(time_data),timedelta(hours=time_window)).next()
                            for tval in all_size_sieved_times_pre:
                                if parser.parse(tval) < localized_time_range.start:
                                    continue 
                                elif tval in localized_time_range:
                                    ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                    indiv_ind_modified_new = ind[ind_time_new]
                                    indiv_ind_modified_list.append(indiv_ind_modified_new)
                                    localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()                                
                                else:
                                    ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                    indiv_ind_modified_new = ind[ind_time_new]                                
                                    next_orig_index = np.where(np.array(fetch_indices_product_orig) == np.array(indiv_ind_modified_new))[0]
                                    if len(next_orig_index) != 0:
                                        indiv_ind_modified_new = fetch_indices_product_orig[next_orig_index[0]]
                                        ind_next_index = np.where(np.array(ind) == indiv_ind_modified_new)[0][0]
                                        tval = all_size_sieved_times_pre[ind_next_index]
                                        indiv_ind_modified_list.append(indiv_ind_modified_new)
                                        localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()
                                    else:
                                        ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                        indiv_ind_modified_new = ind[ind_time_new]
                                        indiv_ind_modified_list.append(indiv_ind_modified_new)
                                        localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()

                            if indiv_ind_modified_list:
                                fetch_indices_product = list(np.zeros(i+1)) + list(indiv_ind_modified_list)
                            else:
                                fetch_indices_product = list(np.zeros(i+1))
                            break


                        elif 'EIT' in base:
                            time_data = product_results.get_response(0)[int(index)]['time']['start']
                            writefits(f'{base}_{mission}/{mission}_{base}_{time_data}_{image_size_output}', reduced_product_data, header_next_good, home_dir)
                            all_time_window_sieved_times_product_times.append(time_data) #(time_val) #unsorted time location
                            all_time_window_sieved_times_product_times_inds_list.append(index)
                            os.remove(query_result_next[0])                        
                        
                            indiv_ind_modified_list = []
                            localized_time_range = TimeRange(str(time_data),timedelta(hours=time_window)).next()
                            for tval in all_size_sieved_times_pre:
                                if parser.parse(tval) < localized_time_range.start:
                                    continue 
                                elif tval in localized_time_range:
                                    ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                    indiv_ind_modified_new = ind[ind_time_new]
                                    indiv_ind_modified_list.append(indiv_ind_modified_new)
                                    localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()                                
                                else:
                                    ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                    indiv_ind_modified_new = ind[ind_time_new]                                
                                    next_orig_index = np.where(np.array(fetch_indices_product_orig) == np.array(indiv_ind_modified_new))[0]
                                    if len(next_orig_index) != 0:
                                        indiv_ind_modified_new = fetch_indices_product_orig[next_orig_index[0]]
                                        ind_next_index = np.where(np.array(ind) == indiv_ind_modified_new)[0][0]
                                        tval = all_size_sieved_times_pre[ind_next_index]
                                        indiv_ind_modified_list.append(indiv_ind_modified_new)
                                        localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()
                                    else:
                                        ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                        indiv_ind_modified_new = ind[ind_time_new]
                                        indiv_ind_modified_list.append(indiv_ind_modified_new)
                                        localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()

                            if indiv_ind_modified_list:
                                fetch_indices_product = list(np.zeros(i+1)) + list(indiv_ind_modified_list)
                            else:
                                fetch_indices_product = list(np.zeros(i+1))
                            break                            
                            
                        elif 'LASCO' in base:
                            if (not planet_comet_transient_filter(base, data_product)) and (not cosmic_ray_filter(base, reduced_product_data)): #if line list empty then can use LASCO image.
                                time_data = product_results.get_response(0)[int(index)]['time']['start']
                                writefits(f'{base}_{mission}/{mission}_{base}_{time_data}_{image_size_output}', reduced_product_data, header_next_good, home_dir)
                                all_time_window_sieved_times_product_times.append(time_data) #(time_val) #unsorted time location
                                all_time_window_sieved_times_product_times_inds_list.append(index)
                                os.remove(query_result_next[0])                        
                        
                                indiv_ind_modified_list = []
                                localized_time_range = TimeRange(str(time_data),timedelta(hours=time_window)).next()
                                for tval in all_size_sieved_times_pre:
                                    if parser.parse(tval) < localized_time_range.start:
                                        continue 
                                    elif tval in localized_time_range:
                                        ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                        indiv_ind_modified_new = ind[ind_time_new]
                                        indiv_ind_modified_list.append(indiv_ind_modified_new)
                                        localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()                                
                                    else:
                                        ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                        indiv_ind_modified_new = ind[ind_time_new]                                
                                        next_orig_index = np.where(np.array(fetch_indices_product_orig) == np.array(indiv_ind_modified_new))[0]
                                        if len(next_orig_index) != 0:
                                            indiv_ind_modified_new = fetch_indices_product_orig[next_orig_index[0]]
                                            ind_next_index = np.where(np.array(ind) == indiv_ind_modified_new)[0][0]
                                            tval = all_size_sieved_times_pre[ind_next_index]
                                            indiv_ind_modified_list.append(indiv_ind_modified_new)
                                            localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()
                                        else:
                                            ind_time_new = np.where(np.array(all_size_sieved_times_pre) == tval)[0][0]
                                            indiv_ind_modified_new = ind[ind_time_new]
                                            indiv_ind_modified_list.append(indiv_ind_modified_new)
                                            localized_time_range = TimeRange(str(tval),timedelta(hours=time_window)).next()

                                if indiv_ind_modified_list:
                                    fetch_indices_product = list(np.zeros(i+1)) + list(indiv_ind_modified_list)
                                else:
                                    fetch_indices_product = list(np.zeros(i+1))
                                break                                                      
                    
                            elif (planet_comet_transient_filter(base, data_product)) or (cosmic_ray_filter(base, reduced_product_data)): #if line list not empty then can't use LASCO image.                    
                                time_data = product_results.get_response(0)[int(index)]['time']['start']
                                
                                if not cosmic_ray_filter(base, reduced_product_data):
                                    transient_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                                    transients_product_list.append((transient_loc, str(time_data)))
                                elif not planet_comet_transient_filter(base, data_product):
                                    blob_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                                    blobs_product_list.append((blob_loc, str(time_data)))                                
                                
                                os.remove(query_result_next[0])                        
                                continue                     
                    
                    
                    elif holes(query_result_next[0],base,mission): #so if True, if there are holes
                        
                        if ('MDI' in base) or (mission == 'SDO'):
                            time_data = all_size_sieved_times_pre[index]                        
                            hole_loc = product_results.data.record[index].split('{')[0]
                            #hole_loc = product_results.T_REC[index] #product_results.get_response(0)[int(index)]['T_REC']                        
                        
                        else:
                            time_data = product_results.get_response(0)[int(index)]['time']['start']
                            hole_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                        
                        holes_product_list.append((hole_loc, str(time_data)))
                        os.remove(query_result_next[0])                        
                        continue 

                elif (data_next_good is None) or (axis1_next_good != axis2_next_good) or (axisnum_next_good != 2):
                    
                    if ('MDI' in base) or (mission == 'SDO'):
                        unreadable_file_ids_product_list.append(product_results.data.record[index].split('{')[0])                                            
                        #unreadable_file_ids_product_list.append(product_results.T_REC[index]) #product_results.get_response(0)[int(index)]['T_REC'])                    
                    
                    else:
                        unreadable_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                        unreadable_file_ids_product_list.append(unreadable_loc)
                    
                    os.remove(query_result_next[0])                    
                    continue
    
    all_time_window_sieved_times_product_times_modified = all_time_window_sieved_times_product_times

    return all_time_window_sieved_times_product_times_modified, holes_product_list, transients_product_list, blobs_product_list, unreadable_file_ids_product_list 


"""
The times corresponding to all fits files that passed all tests are written to csv files.
"""
def csv_writer(base, home_dir, date_start, date_finish, flag, time_window, image_size_output, all_time_window_sieved_times_sorted, lev1_LASCO, mission):
    with open(f'{home_dir}{date_start}_to_{date_finish}_{base}_times_{flag}_{time_window}_LASCOlev1-{lev1_LASCO}_{mission}_{image_size_output}.csv', 'a') as f: #appending lines so not overwriting the file
        writer = csv.writer(f, delimiter='\n')
        writer.writerow(all_time_window_sieved_times_sorted)

