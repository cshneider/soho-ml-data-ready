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
from astropy.io.fits import Header

import h5py

import csv


import json

from tqdm import tqdm

import Mission_utility.sdo_mdi as sdo
import Mission_utility.soho_other as soho

def product_search(time_range, client, time_window):
    ts = '_'.join(str(time_range.start).split(' '))+'_TAI'
    tf = '_'.join(str(time_range.end).split(' '))+'_TAI'
        
    print(ts, tf)
    product_results = client.export(f'mdi.fd_M_96m_lev182[{ts}-{tf}]')
    
    return product_results, client








"""
Encoder to solve TypeError: Object of type 'int64' is not JSON serializable from stackoverflow 
"""
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)

"""
Checks that a fits file contains data otherwise returns unequal axes which results in such a file being discarded. 
"""
def readfits(filename):
    hdrs = []
    try:
        ft = fits.open(filename, memmap=False)
        ft.verify('silentfix+ignore') #ft.verify('fix')-> this gives Unfixable error for card
        
        ##hdr = ft[0].header #[0] for the original large file, #when file written, there is a second header [1] which contains the new axes info automatically
        ##data = ft[0].data        
        for i in range(len(ft)):
            hdrs.append(ft[i].header)
            
            hdr = ft[i].header
            data = ft[i].data
            
        axis1 = hdr['naxis1'] #hdr will be hdr[0] for original large file and hdr[1] for downsampled file since comes from the for loop with i=1 after downsampling!
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
def writefits(filename, data, hdrs, home_dir): #hdrs
    if not os.path.exists(f'{home_dir}{filename}.fits'):
        hdrprimary = fits.PrimaryHDU(header=hdrs[len(hdrs)-1]) #header=hdrs[len(hdrs)-1]) #header=hdr
        fitsimg = fits.ImageHDU(data)
        hdul = fits.HDUList([hdrprimary, fitsimg])
        hdul.writeto(f'{home_dir}{filename}.fits', 'silentfix+ignore') #output_verify='ignore' i think would also work
        
        hdul.close() #added this line

"""
Checks original image for missing pixel data.
"""
def holes(filename,BaseClass):
   
    ft = fits.open(filename, memmap=False)  
    ft.verify('silentfix+ignore') #ft.verify('fix')-> this gives Unfixable error for card 
 
    hdr = ft[BaseClass.holes_method_num].header
    data = ft[BaseClass.holes_method_num].data
    
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
    
    else:
        return BaseClass.are_holes(data, x_coord, filename, rsquared, blank_val, missing_vals)
        
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
def prev_time_resumer(home_dir, time_range_orig, date_time_end, BaseClass): 
#CAN RE-RUN PROGRAM FROM THE LAST DATE ON WHICH STOPPED; WILL PICK UP THE TIMES THAT ARE PRESENT AND CONTINUE! For both resuming on same day and next day.
### CHECKS WHETHER THE START DAY THAT ENTERED IS ALREADY CONTAINED IN THE FILES OF PREVIOUS DAY AND WILL SET START_DATE FROM THAT EXACT TIME! 
### ALSO WORKS IF START ON A NEW DAY AND ARE LOOKING BACK ON THE PREVIOUS DAY
    
    filepath = home_dir + BaseClass.base_full + '_' + BaseClass.mission + '/'
    #filepath = home_dir + 'MDI_96m' + '_' + 'SOHO' + '/'

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
def data_name_selector(home_dir, date_start, date_finish, BaseClass):

    filepath = home_dir + BaseClass.base_full + '_' + BaseClass.mission + '/'

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
Downsample FITS header using the FITS convention correspomnding to the downsampling of image
Based on Ajay's code for metadata seeding HDF5 cubes
"""
def downsample_header(header_content, image_size_output, BaseClass):

    header_content_new = header_content[0].copy() #header_content[0] because header_content is a list of two headers, [0] is original file and [1] is the downsampled file header
    
    orig_img_size = BaseClass.orig_img_size
    
    rescale_factor = int(orig_img_size / image_size_output)
    print('rescale_factor:', rescale_factor)
    
    if (('MDI' in BaseClass.base_full) or ('HMI' in BaseClass.base_full)) and ((BaseClass.fits_headers == 'N') or (BaseClass.fits_headers == 'n')):
        header_content_new['COMMENT'] = f'Zeros outside solar disk for {BaseClass.base_full}'
    
    else:
        header_content_new.update(CDELT1 = header_content_new['CDELT1']*rescale_factor, CDELT2 = header_content_new['CDELT2']*rescale_factor)
        header_content_new.update(CRPIX1 = header_content_new['CRPIX1']/rescale_factor, CRPIX2 = header_content_new['CRPIX2']/rescale_factor)
        header_content_new['COMMENT'] = f'testing that can comment on this {BaseClass.base_full}'
        
        try: 
            header_content_new.update(RSUN_OBS = header_content_new['RSUN_OBS']/rescale_factor, R_SUN = header_content_new['R_SUN']/rescale_factor)
            header_content_new.update(X0 = header_content_new['X0']/rescale_factor, Y0 = header_content_new['Y0']/rescale_factor)
            header_content_new.update(CROP_RAD = header_content_new['CROP_RAD']/rescale_factor)
            header_content_new.update(SOLAR_R = header_content_new['SOLAR_R']/rescale_factor)
                
        except KeyError:
            pass            

    return header_content_new

"""
Downsample FITS header using the FITS convention correspomnding to the downsampling of image 
for the case of retroactively seeding MDI and HMI data cubes with metadata since fits protocol is presently time consuming on JSOC 
"""
def downsample_header_local(image_size_output, query, mag_keys, BaseClass):
    
    orig_img_size = BaseClass.orig_img_size

    rescale_factor = int(orig_img_size / image_size_output)

    for key in mag_keys:
        if (key == 'CDELT1') or (key == 'CDELT2'):
            query[key] = query[key]*rescale_factor #this updates the original data frame, #originally had query[key][0] with [0] here and on every variable here below
        elif (key == 'CRPIX1') or (key == 'CRPIX2'):
            query[key] = query[key]/rescale_factor
        
        try: 
            query['RSUN_OBS'] = query['RSUN_OBS']/rescale_factor 
            query['R_SUN'] = query['R_SUN']/rescale_factor
            query['X0'] = query['X0']/rescale_factor 
            query['Y0'] = query['Y0']/rescale_factor 
            query['CROP_RAD'] = query['CROP_RAD']/rescale_factor
            query['SOLAR_R'] = query['SOLAR_R']/rescale_factor 
             
        except KeyError:
            pass       
    
    return query

"""
Generated a compressed h5py data cube from all fits files present in a product folder
"""
def data_cuber(home_dir, date_start, date_finish, flag, time_window, image_size_output, BaseClass):

    filepath = home_dir + BaseClass.base_full + '_' + BaseClass.mission + '/'

    data_files_pre_pre = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    data_files_pre = [f for f in data_files_pre_pre if 'fits' in f]
    data_files = np.sort(data_files_pre) #to have chronological order and to sink order with list of individual product times

    data_content_list = []
    
    header_down_list = []
    
    for i,elem in tqdm(enumerate(data_files)):
        axdim1,axdim2,data_content,header_content,axisnum = readfits(f'{filepath}{elem}')
        hdr_down = downsample_header(header_content, image_size_output, BaseClass) # hdr_down_keys, hdr_down_vals
        if BaseClass.mission in elem:
            data_content_list.append(data_content)
            header_down_list.append(hdr_down)

    if data_content_list:
        data_content_stack = np.stack(data_content_list)
    else:
        data_content_stack = [] 
                  
    time_start_name_new, time_finish_name_new = data_name_selector(home_dir, date_start, date_finish, BaseClass)
        
    data_cube = h5py.File(f'{home_dir}{time_start_name_new}_to_{time_finish_name_new}_{BaseClass.base_full}_{flag}_{time_window}_LASCOlev1-{BaseClass.lev1_lasco}_{BaseClass.mission}_{image_size_output}_metadata.h5', 'w')
    ##########data_cube_group = data_cube.create_group(f'{base}_{mission}_{image_size_output}')
    data_cube.create_dataset(f'{BaseClass.base_full}_{BaseClass.mission}_{image_size_output}', data=data_content_stack, compression="gzip")
    
    print('len(header_down_list):', len(header_down_list))
    ##########data_cube_group.create_dataset('header', data=header_down_stack.tostring(), compression="gzip")
    
    meta_data_dict = {}
    for i, head in tqdm(enumerate(header_down_list)): #range(len(header_down_list)):
        for j,key in enumerate(list(header_down_list[i].keys())): #list(header_content_new.keys())):
            if (key == 'COMMENT') or (key == 'HISTORY'): #since 'COMMENTS' and 'HISTORY' can occur multiple times so modifying with 'key_counter' to make them unique per image slice
                key = f'{key}{i}' #{key_counter}                  
            meta_data_dict[f'{key}_{i}'] = list(header_down_list[i].values())[j]
    ###data_cube.create_dataset(f'{base}_{mission}_{image_size_output}_header', data=header_down_stack, compression="gzip") #first try approach
    
    data_cube.create_dataset(f'{BaseClass.base_full}_{BaseClass.mission}_{image_size_output}_metadata', data=json.dumps(meta_data_dict, cls=NpEncoder))
    data_cube.attrs['NOTE'] = 'JSON serialization'    
    
    data_cube.close()
    
    return data_cube


   
"""
From the appropriate file size image indices, this function picks up those times and their corresponding indices that are seperated by the user defined time window. 
The indices returned by this function will be used in the product object returned by Fido search.
"""
def fetch_indices(ind,product_results,time_window,prev_time, BaseClass):
    
    all_size_sieved_times_pre = [] #local list to populate at each loop
    all_time_window_sieved_times_product_times = []  #local list to populate at each loop
    
    if BaseClass.class_type == 'SDO_MDI':
        
        drms_export_list = product_results.data.record  
        all_size_sieved_times_pre = [drms.to_datetime(elem.split('[')[1].split(']')[0]).strftime('%Y%m%d%H%M%S') for i,elem in enumerate(drms_export_list)] #only way to get to times from export object
        
                    
    else:
        
        for value in ind:
            all_size_sieved_times_pre.append(product_results[0,:][int(value)]['time']['start'])
            
    all_size_sieved_times = list(np.unique(all_size_sieved_times_pre))
    all_size_sieved_times_aug = prev_time + all_size_sieved_times #prev_time = [] for the very first loop and [last best time from previous loop] for subsequent loops.

    if all_size_sieved_times_aug:
        if prev_time:
            local_time_range = TimeRange(all_size_sieved_times_aug[0],timedelta(hours=time_window)).next() #next() is the important difference here.
        else:
            local_time_range = TimeRange(all_size_sieved_times_aug[0],timedelta(hours=time_window))       
        for i,time_value in enumerate(all_size_sieved_times_aug):
            if time_value in local_time_range or parser.parse(time_value) > local_time_range.end:
                all_time_window_sieved_times_product_times.append(time_value)
                local_time_range = TimeRange(time_value,timedelta(hours=time_window)).next() #important distinction between this local_time_range and the intializing one is the presence of time_value 
            else:
                continue        
    new_inds = [np.where(np.array(all_size_sieved_times_pre) == entry)[0][0] for entry in all_time_window_sieved_times_product_times]
    
    fetch_indices_product = ind[new_inds]
    
    return all_size_sieved_times_pre, fetch_indices_product


def distiller_with_holes(i, query_result_next,fetch_indices_product_orig, all_size_sieved_times_pre, ind, index,time_window, all_time_window_sieved_times_product_times, all_time_window_sieved_times_product_times_inds_list, time_data,image_size_output, reduced_product_data, header_next_good, home_dir, BaseClass):    
    writefits(f'{BaseClass.base_full}_{BaseClass.mission}/{BaseClass.mission}_{BaseClass.base_full}_{time_data}_{image_size_output}', reduced_product_data, header_next_good, home_dir)
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
    return fetch_indices_product                          


"""
Checks the downloaded fits files for holes and discards them if holes are found. Repeats procedure at each time as long as an image contained missing pixels. 
"""
def product_distiller(fetch_indices_product_orig, all_size_sieved_times_pre, ind, product_results, look_ahead, time_window, url_prefix, flag, image_size_output, home_dir, email, client, BaseClass):

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
        query_result = BaseClass.product_retriever(product_results,indiv_ind,url_prefix,home_dir,email,all_size_sieved_times_pre,client)
        axis1_product, axis2_product, data_product, header_product, axisnum_product = readfits(query_result[0])
            
        if (data_product is not None) and (axis1_product == axis2_product) and (axisnum_product == 2):

            if not holes(query_result[0],BaseClass): #so if not True; so no holes; can use image
                reduced_product_data = data_reducer(data_product,flag,image_size_output,axis1_product)
                
                if BaseClass.class_type == 'SDO_MDI':
                    time_data = all_size_sieved_times_pre[indiv_ind]
                    writefits(f'{BaseClass.base_full}_{BaseClass.mission}/{BaseClass.mission}_{BaseClass.base_full}_{time_data}_{image_size_output}', reduced_product_data, header_product, home_dir)
                    os.remove(query_result[0]) #delete original downloaded file
                    all_time_window_sieved_times_product_times.append(time_data)
                    all_time_window_sieved_times_product_times_inds_list.append(indiv_ind)

                elif 'EIT' in BaseClass.base_full:
                    time_data = product_results[0,:][int(indiv_ind)]['time']['start']
                    writefits(f'{BaseClass.base_full}_{BaseClass.mission}/{BaseClass.mission}_{BaseClass.base_full}_{time_data}_{image_size_output}', reduced_product_data, header_product, home_dir)
                    os.remove(query_result[0]) #delete original downloaded file
                    all_time_window_sieved_times_product_times.append(time_data)
                    all_time_window_sieved_times_product_times_inds_list.append(indiv_ind)
                
                
                elif 'LASCO' in BaseClass.base_full:
                    time_data = product_results[0,:][int(indiv_ind)]['time']['start']                    
                    if (not BaseClass.planet_comet_transient_filter(data_product)): #if both line list and blob lost is empty then can use LASCO image.

                        writefits(f'{BaseClass.base_full}_{BaseClass.mission}/{BaseClass.mission}_{BaseClass.base_full}_{time_data}_{image_size_output}', reduced_product_data, header_product, home_dir)
                        os.remove(query_result[0]) #delete original downloaded file
                        all_time_window_sieved_times_product_times.append(time_data)
                        all_time_window_sieved_times_product_times_inds_list.append(indiv_ind)   
                        
                    else:
                        
                        transient_loc = url_prefix + product_results[0,:][int(indiv_ind)]['fileid'] 
                        transients_product_list.append((transient_loc, str(time_data)))  
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
                            query_result_next = BaseClass.product_retriever(product_results,index,url_prefix,home_dir,email,all_size_sieved_times_pre,client)
                            axis1_next_good,axis2_next_good,data_next_good,header_next_good,axisnum_next_good = readfits(query_result_next[0])

                            if (data_next_good is not None) and (axis1_next_good == axis2_next_good) and (axisnum_next_good == 2):
                                reduced_product_data = data_reducer(data_next_good,flag,image_size_output,axis1_next_good)                                    
                                
                                if (not holes(query_result_next[0],BaseClass)) and (not BaseClass.planet_comet_transient_filter(data_next_good)):
                                    
                                    time_data = product_results[0,:][int(index)]['time']['start']
                                    fetch_indices_product, ind_time_new, indiv_ind_modified_new, indiv_ind_modified_list, localized_time_range = distiller_with_holes(i, query_result_next,fetch_indices_product_orig, all_size_sieved_times_pre, ind, index,time_window, all_time_window_sieved_times_product_times, all_time_window_sieved_times_product_times_inds_list, time_data,image_size_output, reduced_product_data, header_next_good, home_dir, BaseClass)
                                                              
                                    break

                                elif (holes(query_result_next[0],BaseClass)) or (BaseClass.planet_comet_transient_filter(data_next_good)): #so if True, if there are holes
                                    
                                    time_data = product_results[0,:][int(index)]['time']['start']                            
                                    
                                    if (not BaseClass.planet_comet_transient_filter(data_next_good)): #so no transient
                                        hole_loc = url_prefix + product_results[0,:][int(index)]['fileid']
                                        holes_product_list.append((hole_loc, str(time_data)))
                                    else:
                                        transient_loc = url_prefix + product_results[0,:][int(index)]['fileid']
                                        transients_product_list.append((transient_loc, str(time_data)))                                  
                                    
                                    os.remove(query_result_next[0])
                                    continue 

                            elif (data_next_good is None) or (axis1_next_good != axis2_next_good) or (axisnum_next_good != 2):
                                    
                                unreadable_loc = url_prefix + product_results[0,:][int(index)]['fileid']
                                unreadable_file_ids_product_list.append(unreadable_loc)
                                
                                os.remove(query_result_next[0])
                                continue
                
            elif holes(query_result[0],BaseClass): #so if True, if there are holes
                if BaseClass.class_type == 'SDO_MDI':
                    time_data = all_size_sieved_times_pre[indiv_ind]                                                        
                    hole_loc = product_results.data.record[indiv_ind].split('{')[0]                                      
                
                else:
                    time_data = product_results[0,:][int(indiv_ind)]['time']['start'] 
                    hole_loc = url_prefix + product_results[0,:][int(indiv_ind)]['fileid']                       
                
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
                    query_result_next = BaseClass.product_retriever(product_results,index,url_prefix,home_dir,email,all_size_sieved_times_pre,client)
                    axis1_next_good,axis2_next_good,data_next_good,header_next_good,axisnum_next_good = readfits(query_result_next[0])

                    if (data_next_good is not None) and (axis1_next_good == axis2_next_good) and (axisnum_next_good == 2):
                        
                        reduced_product_data = data_reducer(data_next_good,flag,image_size_output,axis1_next_good)
                        
                        if not holes(query_result_next[0],BaseClass): #and (not planet_comet_transient_filter(base, data_product)) and (not cosmic_ray_filter(base, reduced_product_data)):
                            
                            if BaseClass.class_type == 'SDO_MDI':
                                time_data = all_size_sieved_times_pre[index]     
                                fetch_indices_product = distiller_with_holes(i, query_result_next,fetch_indices_product_orig, all_size_sieved_times_pre, ind, index,time_window, all_time_window_sieved_times_product_times, all_time_window_sieved_times_product_times_inds_list, time_data,image_size_output, reduced_product_data, header_next_good, home_dir, BaseClass)
                                                           
                                break                            
                            
                            elif 'EIT' in BaseClass.base_full:
                                time_data = product_results[0,:][int(index)]['time']['start']
                                fetch_indices_product = distiller_with_holes(i, query_result_next,fetch_indices_product_orig, all_size_sieved_times_pre, ind, index,time_window, all_time_window_sieved_times_product_times, all_time_window_sieved_times_product_times_inds_list, time_data,image_size_output, reduced_product_data, header_next_good, home_dir, BaseClass)                            
                                break                             
                            
                            elif 'LASCO' in BaseClass.base_full:    
                                time_data = product_results[0,:][int(index)]['time']['start']                            
                                if (not BaseClass.planet_comet_transient_filter(data_next_good)):
                                    fetch_indices_product = distiller_with_holes(i, query_result_next,fetch_indices_product_orig, all_size_sieved_times_pre, ind, index,time_window, all_time_window_sieved_times_product_times, all_time_window_sieved_times_product_times_inds_list, time_data,image_size_output, reduced_product_data, header_next_good, home_dir, BaseClass)                           
                                    break                                                         
                                
                                
                                elif (BaseClass.planet_comet_transient_filter(data_product)):                                

                                    transient_loc = url_prefix + product_results[0,:][int(index)]['fileid'] 
                                    transients_product_list.append((transient_loc, str(time_data)))
                                    
                                    os.remove(query_result_next[0]) #these were indented one outside on 27.03.21
                                    continue    #these were indented one outside on 27.03.21                      
                                
                        elif holes(query_result_next[0],BaseClass): #so if True, if there are holes
                            
                            if BaseClass.class_type == 'SDO_MDI':
                                time_data = all_size_sieved_times_pre[index]
                                hole_loc = product_results.data.record[index].split('{')[0]                          
                            
                            else:
                                time_data = product_results[0,:][int(index)]['time']['start']
                                hole_loc = url_prefix + product_results[0,:][int(index)]['fileid']
                            
                            holes_product_list.append((hole_loc, str(time_data)))
                            
                            os.remove(query_result_next[0])
                            continue 

                    elif (data_next_good is None) or (axis1_next_good != axis2_next_good) or (axisnum_next_good != 2):
                        
                        if BaseClass.class_type == 'SDO_MDI':
                            unreadable_file_ids_product_list.append(product_results.data.record[index].split('{')[0])                      
                        
                        else:
                            unreadable_loc = url_prefix + product_results[0,:][int(index)]['fileid']
                            unreadable_file_ids_product_list.append(unreadable_loc)
                        
                        os.remove(query_result_next[0])
                        continue


        elif (data_product is None) or (axis1_product != axis2_product) or (axisnum_product != 2):
            
            if BaseClass.class_type == 'SDO_MDI':
                unreadable_file_ids_product_list.append(product_results.data.record[indiv_ind].split('{')[0])
                bad_time_val = all_size_sieved_times_pre[indiv_ind]       
            
            else:
                unreadable_loc = url_prefix + product_results[0,:][int(indiv_ind)]['fileid']
                unreadable_file_ids_product_list.append(unreadable_loc)
                bad_time_val = product_results[0,:][int(indiv_ind)]['time']['start']
            
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
                query_result_next = BaseClass.product_retriever(product_results,index,url_prefix,home_dir,email,all_size_sieved_times_pre,client)
                axis1_next_good,axis2_next_good,data_next_good,header_next_good,axisnum_next_good = readfits(query_result_next[0])

                if (data_next_good is not None) and (axis1_next_good == axis2_next_good) and (axisnum_next_good == 2):

                    if not holes(query_result_next[0],BaseClass): #so if not True; so no holes; can use image
                        reduced_product_data = data_reducer(data_next_good,flag,image_size_output,axis1_next_good)
                        
                        if BaseClass.class_type == 'SDO_MDI':
                            time_data = all_size_sieved_times_pre[index]                      
                            fetch_indices_product = distiller_with_holes(i, query_result_next,fetch_indices_product_orig, all_size_sieved_times_pre, ind, index,time_window, all_time_window_sieved_times_product_times, all_time_window_sieved_times_product_times_inds_list, time_data,image_size_output, reduced_product_data, header_next_good, home_dir, BaseClass)
                            break


                        elif 'EIT' in BaseClass.base_full:
                            time_data = product_results[0,:][int(index)]['time']['start']
                            fetch_indices_product = distiller_with_holes(i, query_result_next,fetch_indices_product_orig, all_size_sieved_times_pre, ind, index,time_window, all_time_window_sieved_times_product_times, all_time_window_sieved_times_product_times_inds_list, time_data,image_size_output, reduced_product_data, header_next_good, home_dir, BaseClass)
                            break                            
                            
                        elif 'LASCO' in BaseClass.base_full:
                            if (not BaseClass.planet_comet_transient_filter(data_product)): #if line list empty then can use LASCO image.
                                time_data = product_results[0,:][int(index)]['time']['start']
                                fetch_indices_product = distiller_with_holes(i, query_result_next,fetch_indices_product_orig, all_size_sieved_times_pre, ind, index,time_window, all_time_window_sieved_times_product_times, all_time_window_sieved_times_product_times_inds_list, time_data,image_size_output, reduced_product_data, header_next_good, home_dir, BaseClass)
                                break                                                      
                    
                            elif (BaseClass.planet_comet_transient_filter(data_product)): #if line list not empty then can't use LASCO image.                    
                                time_data = product_results[0,:][int(index)]['time']['start']
                                
                                transient_loc = url_prefix + product_results[0,:][int(index)]['fileid']
                                transients_product_list.append((transient_loc, str(time_data)))                              
                                
                                os.remove(query_result_next[0])                        
                                continue                     
                    
                    
                    elif holes(query_result_next[0],BaseClass): #so if True, if there are holes
                        
                        if BaseClass.class_type == 'SDO_MDI':
                            time_data = all_size_sieved_times_pre[index]                        
                            hole_loc = product_results.data.record[index].split('{')[0]                     
                        
                        else:
                            time_data = product_results[0,:][int(index)]['time']['start']
                            hole_loc = url_prefix + product_results[0,:][int(index)]['fileid']
                        
                        holes_product_list.append((hole_loc, str(time_data)))
                        os.remove(query_result_next[0])                        
                        continue 

                elif (data_next_good is None) or (axis1_next_good != axis2_next_good) or (axisnum_next_good != 2):
                    
                    if BaseClass.class_type == 'SDO_MDI':
                        unreadable_file_ids_product_list.append(product_results.data.record[index].split('{')[0])                   
                    
                    else:
                        unreadable_loc = url_prefix + product_results[0:,][int(index)]['fileid']
                        unreadable_file_ids_product_list.append(unreadable_loc)
                    
                    os.remove(query_result_next[0])                    
                    continue
    
    all_time_window_sieved_times_product_times_modified = all_time_window_sieved_times_product_times

    return all_time_window_sieved_times_product_times_modified, holes_product_list, transients_product_list, blobs_product_list, unreadable_file_ids_product_list 


"""
The times corresponding to all fits files that passed all tests are written to csv files.
"""
def csv_writer(home_dir, date_start, date_finish, flag, time_window, image_size_output, all_time_window_sieved_times_sorted, BaseClass):
    with open(f'{home_dir}{date_start}_to_{date_finish}_{BaseClass.base_full}_times_{flag}_{time_window}_LASCOlev1-{BaseClass.lev1_lasco}_{BaseClass.mission}_{image_size_output}.csv', 'a') as f: #appending lines so not overwriting the file
        writer = csv.writer(f, delimiter='\n')
        writer.writerow(all_time_window_sieved_times_sorted)

