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

from sunpy.net import Fido
from sunpy.net.vso import attrs as avso
from sunpy.time import TimeRange

import astropy.units as u
from astropy.io import fits

import h5py

import csv

"""
Checks that a fits file contains data otherwise returns unequal axes which results in such a file being discarded. 
"""
def readfits(filename):
    try:
        ft = fits.open(filename, memmap=False)
        hdr = ft[0].header
        data = ft[0].data
        axis1 = hdr['naxis1']
        axis2 = hdr['naxis2']
        axisnum = hdr['naxis']
        ft.close()
    
    except ValueError:
        axis1 = 1
        axis2 = 2
        axisnum = 3
        data = None
            
    return axis1, axis2, data, axisnum

"""
Writes a fits file
"""
def writefits(filename, data, home_dir):
    if not os.path.exists(f'{home_dir}{filename}.fits'):
        fitsname = fits.PrimaryHDU(data)
        fitsname.writeto(f'{home_dir}{filename}.fits')

"""
Checks original image for missing pixel data.
"""
def holes(filename):
    filename = str(filename)
    
    ft = fits.open(filename, memmap=False)
    hdr = ft[0].header
    data = ft[0].data
    ft.close()

    try:
        x_coord = hdr['CRPIX1']
        y_coord = hdr['CRPIX2']
    
    except KeyError:
        x_coord = hdr['naxis1'] / 2.
        y_coord = hdr['naxis2'] / 2.

    y_ind,x_ind = np.indices((hdr['naxis1'],hdr['naxis2']))
    rsquared = (x_ind - x_coord)**2 + (y_ind - y_coord)**2
    
    matches = ['96m', 'MDI']
    
    if 'efz' in filename: #good for all EIT products 
        rad = x_coord*np.sqrt(2)
        indices = np.where(rsquared.flatten() < rad**2)[0]
        zeros_ind = np.where(data.flatten()[indices] == 0.)[0]
        zeros_ind_len = len(zeros_ind)

        if zeros_ind_len > 100:
            return True #so image not useable as there are holes
        else:
            return False #can use this image
    
    elif any([x in filename for x in matches]):
        rad1 = float(x_coord)
        rad2 = 0.6*float(x_coord)
        indices_rad1 = np.where(rsquared.flatten() < rad1**2)[0]
        indices_rad2 = np.where(rsquared.flatten() < rad2**2)[0]
        zeros_ind = np.where(data.flatten()[indices_rad1] == 0.)[0]
        nan_ind = np.where(data.flatten()[indices_rad2] != data.flatten()[indices_rad2])[0]
        zeros_nan_ind_len = len(list(zeros_ind) + list(nan_ind))
        
        if zeros_nan_ind_len > 100:
            return True #so image not useable as there are holes
        else:
            return False #can use this image

    elif 'LASCO_C3' in filename:
        #print('LASCO_C3')
        rad = 0.8*x_coord
        indices = np.where(rsquared.flatten() < rad**2)[0]
        zeros_ind = np.where(data.flatten()[indices] == 0.)[0]
        zeros_ind_len = len(zeros_ind)  
        
        if zeros_ind_len > 100:
            return True #so image not useable as there are holes
        else:
            return False #can use this image   

    
    elif 'LASCO_C2' in filename:
        #print('LASCO_C2')
        rad1 = 160 #this seems good
        #print('rad1:', rad1)
        rad2 = int(x_coord)
        indices = np.where((rad2**2 > rsquared.flatten()) & (rsquared.flatten() > rad1**2))[0]
        zeros_ind = np.where(data.flatten()[indices] == 0.)[0]
        zeros_ind_len = len(zeros_ind)
     
        if zeros_ind_len > 100:
            return True #so image not useable as there are holes
        else:
            return False #can use this image
        
"""
Reduces the image dimensions via one of the following four methods: subsampling, interpolation, min pooling, or max pooling.
"""
def data_reducer(data,flag,target_dimension,axis1_shape):
    scale_factor = int(axis1_shape/target_dimension)
    
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
def prev_time_resumer(home_dir, base, time_range_orig, date_time_end): 
#CAN RE-RUN PROGRAM FROM THE LAST DATE ON WHICH STOPPED; WILL PICK UP THE TIMES THAT ARE PRESENT AND CONTINUE! For both resuming on same day and next day.
### CHECKS WHETHER THE START DAY THAT ENTERED IS ALREADY CONTAINED IN THE FILES OF PREVIOUS DAY AND WILL SET START_DATE FROM THAT EXACT TIME! 
### ALSO WORKS IF START ON A NEW DAY AND ARE LOOKING BACK ON THE PREVIOUS DAY
    
    print('base:', base)
    filepath = home_dir + base + '/'

    data_files_pre = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    data_files = np.sort(data_files_pre)
    
    if len(data_files) != 0:
        prev_time_pre = data_files[-1]
        if 'EIT' in str(prev_time_pre): 
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
def data_name_selector(home_dir, base, date_start, date_finish):

    print('base:', base)
    filepath = home_dir + base + '/'

    data_files_pre = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    data_files = np.sort(data_files_pre)
    print('len(data_files):', len(data_files)) 
    
    if len(data_files) != 0: 
        time_start_name_pre = data_files[0]
        time_finish_name_pre = data_files[-1]         
        
        if 'EIT' in str(time_start_name_pre):
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
def data_cuber(home_dir, base, date_start, date_finish, flag, target_dimension):

    print('base:', base)
    filepath = home_dir + base + '/'

    data_files_pre = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    data_files = np.sort(data_files_pre) #to have chronological order and to sink order with list of individual product times
    print('len(data_files):', len(data_files))

    data_content_list = []
    for elem in data_files:
        axdim1,axdim2,data_content,axisnum = readfits(f'{filepath}{elem}')
        if (axdim1 == axdim2) and ('SOHO' in elem) and (axisnum == 2):
            data_content_list.append(data_content)
        data_content_list.append(data_content)

    if data_content_list:
        data_content_stack = np.stack(data_content_list)
    else:
        data_content_stack = []
                  
    time_start_name_new, time_finish_name_new = data_name_selector(home_dir, base, date_start, date_finish)
        
    data_cube = h5py.File(f'{home_dir}{time_start_name_new}_to_{time_finish_name_new}_{base}_{flag}_{target_dimension}.h5', 'w')
    data_cube.create_dataset(f'{base}_{target_dimension}', data=data_content_stack, compression="gzip")
    data_cube.close()
                            
    return data_cube

"""
Used SunPy's Fido to search for the user specified products
"""
def product_search(base,time_range,date_time_start):
    if 'EIT' in base:
        wavelen = int(base[3:6])
        product_results = Fido.search(avso.Time(time_range,date_time_start),avso.Source('SOHO'),avso.Instrument('EIT'),avso.Provider('SDAC'),avso.Wavelength(wavelen * avso.u.Angstrom, wavelen * avso.u.Angstrom))
    
    elif 'MDI' in base:
        product_results = Fido.search(avso.Time(time_range,date_time_start),avso.Source('SOHO'),avso.Instrument('MDI'),avso.Provider('SDAC'),avso.Physobs('LOS_MAGNETIC_FIELD'))
    
    elif 'LASCO' in base:
        detector = base.split('_')[1]
        product_results = Fido.search(avso.Time(time_range,date_time_start),avso.Provider('SDAC'),avso.Source('SOHO'),avso.Instrument('LASCO'),avso.Detector(detector))
    
    return product_results

"""
Returns the indices corresponding to the proper sizes of each of the data products to ensure that they are single, two dimensional images.
"""
def index_of_sizes(base,product_results):
    
    matches = ['171', '304', '284']
    
    if 'EIT195' in base:
        size_list = [elem['size'] for elem in product_results.get_response(0)[:]]
        print(np.unique(size_list), len(size_list))
        ind_2059 = np.where(np.array(size_list) == 2059)[0]
        ind_523 = np.where(np.array(size_list) == 523)[0]
        print(len(ind_2059))
        print(len(ind_523))
        ind = np.sort(list(ind_2059) + list(ind_523)) #important to sort here since combining two lists!
        print(len(ind))
        
    elif 'MDI' in base:
        size_list = [elem['size'] for elem in product_results.get_response(0)[:]]
        print(np.unique(size_list), len(size_list))
        ind = np.where(np.array(size_list) == 4115.0)[0]
        print(len(ind))        
        
    elif 'LASCO' in base:
        size_list = [int(np.ceil(elem['size'] / 100.0))*100 for elem in product_results.get_response(0)[:]]
        print(np.unique(size_list), len(size_list))
        ind = np.where(np.array(size_list) == 2100.0)[0] 
        print(len(ind))
        
    elif any([x in base for x in matches]):
        size_list = [elem['size'] for elem in product_results.get_response(0)[:]]
        print(np.unique(size_list), len(size_list))
        ind = np.where(np.array(size_list) == 2059)[0]        
        print(len(ind))
        
    return ind
   
"""
From the appropriate file size image indices, this function picks up those times and their corresponding indices that are seperated by the user defined time window. 
The indices returned by this function will be used in the product object returned by Fido search.
"""
def fetch_indices(base,ind,product_results,time_window,look_ahead, prev_time):
    
    all_size_sieved_times_pre = [] #local list to populate at each loop
    all_time_window_sieved_times = [] #local list to populate at each loop

    for value in ind:
        all_size_sieved_times_pre.append(product_results.get_response(0)[int(value)]['time']['start'])
    all_size_sieved_times = list(np.unique(all_size_sieved_times_pre))
    all_size_sieved_times_aug = prev_time + all_size_sieved_times #prev_time = [] for the very first loop and [last best time from previous loop] for subsequent loops.

    for i,time_value in enumerate(all_size_sieved_times_aug):
        local_time_range = TimeRange(str(time_value),timedelta(hours=time_window))

        local_list = []
        for k,time_val in enumerate(all_size_sieved_times_aug[i:i+look_ahead]):
            if time_val in local_time_range:
                local_list.append(time_val)
        if local_list:
            for entry in local_list[1:]:
                all_size_sieved_times_aug.remove(entry)
            all_time_window_sieved_times.append(local_list[0])

    all_time_window_sieved_times_product_times = list(np.unique(all_time_window_sieved_times)) #np.unique() does np.array() and np.sort()

    if not prev_time: #so if no prev_time (i.e., prev_time = [] at the start of the first loop)    
        new_inds = [np.where(np.array(all_size_sieved_times_pre) == entry)[0][0] for entry in all_time_window_sieved_times_product_times]      
    
    elif prev_time:
        new_inds = [np.where(np.array(all_size_sieved_times_pre) == entry)[0][0] for entry in all_time_window_sieved_times_product_times[1:]]                
    
    fetch_indices_product = ind[new_inds]
    
    return all_size_sieved_times_pre, fetch_indices_product
    
"""
Using Fido's returned query object that has now been sieved for proper data size and user specified time window, the fileid is extracted from the accompanying Fido dictionary and wget is used to retrieve the product. Fido has a command Fido.fetch() which can be used to. This method appears to time out after it has been called a certain number of times across the network. wget seems more robust for this purpose. 
"""
def product_retriever(base,product_results,indiv_ind,url_prefix,home_dir):
    
    fileid = product_results.get_response(0)[int(indiv_ind)]['fileid']
    item_wget =  url_prefix + fileid
    cmd = 'wget' + ' ' + item_wget + ' ' + '-P' + ' ' + f'{home_dir}{base}' #OBTAIN TIMEOUT ISSUE WITH FIDO FETCH! SEEMS THAT WITH WGET CAN CIRCUMVENT IT.
    args = shlex.split(cmd)    
    wget_output = str(subprocess.check_output(args)).strip('b')
    
    while wget_output != "''":
        print(f'Encountered wget error with exit status {wget_output} for {item_wget}')
        cmd = 'wget' + ' ' + '--retry-connrefused' + ' ' + '--waitretry=1' + ' ' + '--read-timeout=20' + ' ' + '--timeout=15' + ' ' + '-t' + ' ' + '0' + ' ' + '--continue' + ' ' + item_wget + ' ' + '-P' + ' ' + f'{home_dir}{base}'    
        args = shlex.split(cmd)
        wget_output = str(subprocess.check_output(args)).strip('b')
        if wget_output == "''":
            break
        time.sleep(1) 
    
    downloaded_fileid = fileid.split('/')[-1]
    query_result = [f'{home_dir}{base}/{downloaded_fileid}']
    
    return query_result

"""
Checks the downloaded fits files for holes and discards them if holes are found. Repeats procedure at each time as long as an image contained missing pixels. 
"""
def product_distiller(fetch_indices_product_orig, base, all_size_sieved_times_pre, ind, product_results, look_ahead, time_window, url_prefix, flag, target_dimension, home_dir):

    holes_product_list = []
    unreadable_file_ids_product_list = []

    all_time_window_sieved_times_product_times = [] 
    all_time_window_sieved_times_product_times_inds_list = [] 

    fetch_indices_product = fetch_indices_product_orig.copy()
    
    for i,elem in enumerate(fetch_indices_product): #the i index retains the original number of members of the original fetch_indices_product. The fetch_indices_product list is modified when holes occur.
        if (i > len(fetch_indices_product)-1): #fetch_indices_product is modified in the program to account for times corresponding to holes in images. When all fitting times exhausted then break out of loop.
            break        
        indiv_ind = fetch_indices_product[i]
        query_result = product_retriever(base,product_results,indiv_ind,url_prefix,home_dir)
        axis1_product, axis2_product, data_product, axisnum_product = readfits(query_result[0])
            
        if (data_product is not None) and (axis1_product == axis2_product) and (axisnum_product == 2):

            if not holes(query_result[0]): #so if not True; so no holes; can use image
                reduced_product_data = data_reducer(data_product,flag,target_dimension,axis1_product)
                time_data = product_results.get_response(0)[int(indiv_ind)]['time']['start']
                writefits(f'{base}/SOHO_{base}_{time_data}_{target_dimension}', reduced_product_data, home_dir)
                os.remove(query_result[0]) #delete original downloaded file
                all_time_window_sieved_times_product_times.append(time_data)
                all_time_window_sieved_times_product_times_inds_list.append(indiv_ind)                

            elif holes(query_result[0]): #so if True, if there are holes
                time_data = product_results.get_response(0)[int(indiv_ind)]['time']['start'] 
                hole_loc = url_prefix + product_results.get_response(0)[int(indiv_ind)]['fileid']                       
                holes_product_list.append((hole_loc, str(time_data)))
                hole_time_val = product_results.get_response(0)[int(indiv_ind)]['time']['start']
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
                    query_result_next = product_retriever(base,product_results,index,url_prefix,home_dir)
                    axis1_next_good,axis2_next_good,data_next_good,axisnum_next_good = readfits(query_result_next[0])

                    if (data_next_good is not None) and (axis1_next_good == axis2_next_good) and (axisnum_next_good == 2):

                        if not holes(query_result_next[0]): #so if not True; so no holes; can use image
                            reduced_product_data = data_reducer(data_next_good,flag,target_dimension,axis1_next_good)
                            time_data = product_results.get_response(0)[int(index)]['time']['start']
                            writefits(f'{base}/SOHO_{base}_{time_data}_{target_dimension}', reduced_product_data, home_dir)

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
                                                                        
                            print('indiv_ind_modified_list:', indiv_ind_modified_list)
                            
                            if indiv_ind_modified_list:
                                fetch_indices_product = list(np.zeros(i+1)) + list(indiv_ind_modified_list)                            
                            else:
                                fetch_indices_product = list(np.zeros(i+1))                            
                            break

                        elif holes(query_result_next[0]): #so if True, if there are holes
                            time_data = product_results.get_response(0)[int(index)]['time']['start']
                            hole_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                            holes_product_list.append((hole_loc, str(time_data)))
                            os.remove(query_result_next[0])
                            continue 

                    elif (data_next_good is None) or (axis1_next_good != axis2_next_good) or (axisnum_next_good != 2):
                        unreadable_file_ids_product_list.append(product_results.get_response(0)[int(index)]['fileid'])
                        os.remove(query_result_next[0])
                        continue


        elif (data_product is None) or (axis1_product != axis2_product) or (axisnum_product != 2):
            unreadable_file_ids_product_list.append(product_results.get_response(0)[int(indiv_ind)]['fileid'])
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
                query_result_next = product_retriever(base,product_results,index,url_prefix,home_dir)
                axis1_next_good,axis2_next_good,data_next_good,axisnum_next_good = readfits(query_result_next[0])

                if (data_next_good is not None) and (axis1_next_good == axis2_next_good) and (axisnum_next_good == 2):

                    if not holes(query_result_next[0]): #so if not True; so no holes; can use image
                        reduced_product_data = data_reducer(data_next_good,flag,target_dimension,axis1_next_good)
                        time_data = product_results.get_response(0)[int(index)]['time']['start']
                        writefits(f'{base}/SOHO_{base}_{time_data}_{target_dimension}', reduced_product_data, home_dir)

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

                        print('indiv_ind_modified_list:', indiv_ind_modified_list)
                        
                        if indiv_ind_modified_list:
                            fetch_indices_product = list(np.zeros(i+1)) + list(indiv_ind_modified_list)
                        else:
                            fetch_indices_product = list(np.zeros(i+1))
                        break

                    elif holes(query_result_next[0]): #so if True, if there are holes
                        time_data = product_results.get_response(0)[int(index)]['time']['start']
                        hole_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                        holes_product_list.append((hole_loc, str(time_data)))
                        os.remove(query_result_next[0])                        
                        continue 

                elif (data_next_good is None) or (axis1_next_good != axis2_next_good) or (axisnum_next_good != 2):
                    unreadable_file_ids_product_list.append(product_results.get_response(0)[int(index)]['fileid'])
                    os.remove(query_result_next[0])                    
                    continue
    
    all_time_window_sieved_times_product_times_modified = all_time_window_sieved_times_product_times

    return all_time_window_sieved_times_product_times_modified, holes_product_list, unreadable_file_ids_product_list 


"""
The times corresponding to all fits files that passed all tests are written to csv files.
"""
def csv_writer(base,home_dir,date_start,date_finish,flag,target_dimension, all_time_window_sieved_times_sorted):
    with open(f'{home_dir}{date_start}_to_{date_finish}_{base}_times_{flag}_{target_dimension}.csv', 'a') as f: #appending lines so not overwriting the file
        writer = csv.writer(f, delimiter='\n')
        writer.writerow(all_time_window_sieved_times_sorted)

"""
The unique times corresponding to all fits files that passed all tests are written to csv files to control for duplicate times in the output csv in order to match the exact number of fits files.
"""
def csv_time_uniq_writer(base,home_dir,date_start,date_finish,flag,target_dimension):
    if isfile(f'{home_dir}{date_start}_to_{date_finish}_{base}_times_{flag}_{target_dimension}.csv'):
        with open(f'{home_dir}{date_start}_to_{date_finish}_{base}_times_{flag}_{target_dimension}.csv', 'r') as ff:
            csv_reader = csv.reader(ff, delimiter='\n')
            csv_data = [line for line in csv_reader]
        csv_uniq_times = list(np.unique(csv_data))
        with open(f'{home_dir}{date_start}_to_{date_finish}_{base}_times_{flag}_{target_dimension}_new.csv', 'w') as g:
            writer_new = csv.writer(g, delimiter='\n')
            writer_new.writerow(csv_uniq_times)
