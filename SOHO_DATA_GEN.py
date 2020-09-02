import numpy as np

import os
from os import listdir
from os.path import isfile, join

import shlex, subprocess

from skimage.transform import rescale
from skimage.measure import block_reduce

from datetime import datetime, date, time, timedelta
from dateutil import parser

from sunpy.net import Fido
from sunpy.net.vso import attrs as avso
from sunpy.time import TimeRange #parse_time,

import astropy.units as u
from astropy.io import fits

import h5py

import csv


def readfits(filename):
    try:
        ft = fits.open(filename, memmap=False)
        hdr = ft[0].header
        data = ft[0].data
        axis1 = hdr['naxis1']
        axis2 = hdr['naxis2']
        ft.close()
    
    except ValueError:
        axis1 = 1
        axis2 = 2
        data = None
            
    return axis1,axis2,data

def writefits(filename, data, home_dir):
    if not os.path.exists(f'{home_dir}{filename}.fits'):
        fitsname = fits.PrimaryHDU(data)
        fitsname.writeto(f'{home_dir}{filename}.fits')

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


def data_cuber(home_dir, base, date_start, date_finish, flag, target_dimension):

    print('base:', base)
    filepath = home_dir + base + '/'

    data_files_pre = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    data_files = np.sort(data_files_pre) #to have chronological order and to sink order with list of individual product times
    print('len(data_files):', len(data_files))
    
    data_content_list = []
    for elem in data_files:
        axdim1,axdim2,data_content = readfits(f'{filepath}{elem}')
        if (axdim1 == axdim2) and ('SOHO' in elem):
            data_content_list.append(data_content)
        data_content_list.append(data_content)

    if data_content_list:
        data_content_stack = np.stack(data_content_list)
    else:
        data_content_stack = []
                  
    data_cube = h5py.File(f'{home_dir}{date_start}_to_{date_finish}_{base}_{flag}_{target_dimension}.h5', 'w')
    data_cube.create_dataset(f'{base}_{target_dimension}', data=data_content_stack, compression="gzip")
    data_cube.close()
                            
    return data_cube


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
   

def fetch_indices(base,ind,product_results,time_window,look_ahead):
    
    all_size_sieved_times = [] #local list to populate at each loop
    all_2hr_sieved_times = [] #local list to populate at each loop

    for value in ind:
        all_size_sieved_times.append(product_results.get_response(0)[int(value)]['time']['start'])
    all_size_sieved_times_copy = all_size_sieved_times.copy()

    for i,time_value in enumerate(all_size_sieved_times_copy):
        local_time_range = TimeRange(str(time_value),timedelta(hours=time_window))

        local_list = []
        for k,time_val in enumerate(all_size_sieved_times_copy[i:i+look_ahead]):
            if time_val in local_time_range:
                local_list.append(time_val)
        if local_list:
            for entry in local_list[1:]:
                all_size_sieved_times_copy.remove(entry)
            all_2hr_sieved_times.append(local_list[0])

    all_2hr_sieved_times_product_times = list(np.unique(all_2hr_sieved_times)) #np.unique() does np.array() and np.sort()

    all_2hr_sieved_times_product_times_inds_list_pre = [np.where(np.array(all_size_sieved_times) == item)[0][0] for item in all_2hr_sieved_times_product_times]
    
    if all_2hr_sieved_times_product_times_inds_list_pre:
        all_2hr_sieved_times_product_times_inds_list = list(np.hstack(all_2hr_sieved_times_product_times_inds_list_pre))

    fetch_indices_product = ind[all_2hr_sieved_times_product_times_inds_list] #THESE ARE THE INDICES TO FETCH!
    
    return all_size_sieved_times, all_2hr_sieved_times_product_times, all_2hr_sieved_times_product_times_inds_list, fetch_indices_product
    

def product_retriever(base,product_results,indiv_ind,url_prefix,home_dir):
    
    fileid = product_results.get_response(0)[int(indiv_ind)]['fileid']
    item_wget =  url_prefix + fileid
    cmd = 'wget' + ' ' + item_wget + ' ' + '-P' + ' ' + f'{home_dir}{base}' #OBTAIN TIMEOUT ISSUE WITH FIDO FETCH! SEEMS THAT WITH WGET CAN CIRCUMNAVIGATE IT.
    args = shlex.split(cmd)
    subprocess.check_call(args)
    downloaded_fileid = fileid.split('/')[-1]
    query_result = [f'{home_dir}{base}/{downloaded_fileid}']
    
    return query_result


def product_distiller(base, axis1_product,axis2_product,data_product, all_size_sieved_times, all_2hr_sieved_times_product_times, all_2hr_sieved_times_product_times_inds_list, query_result, ind, indiv_ind, product_results, look_ahead, time_window, url_prefix, flag, target_dimension, home_dir):

    holes_product_list = []
    unreadable_file_ids_product_list = []
    
    if (data_product is not None) and (axis1_product == axis2_product):

        if not holes(query_result[0]): #so if not True; so no holes; can use image
            reduced_product_data = data_reducer(data_product,flag,target_dimension,axis1_product)
            time_data = product_results.get_response(0)[int(indiv_ind)]['time']['start']
            writefits(f'{base}/SOHO_{base}_{time_data}_{target_dimension}', reduced_product_data, home_dir)
            os.remove(query_result[0]) #delete original downloaded file

        elif holes(query_result[0]): #so if True, if there are holes
            time_data = product_results.get_response(0)[int(indiv_ind)]['time']['start'] 
            hole_loc = url_prefix + product_results.get_response(0)[int(indiv_ind)]['fileid']                       
            holes_product_list.append((hole_loc, str(time_data)))
            hole_time_val = product_results.get_response(0)[int(indiv_ind)]['time']['start']
            
            all_2hr_sieved_times_product_times.remove(hole_time_val)
            
            ind_hole_time_val = np.where(np.array(all_size_sieved_times) == hole_time_val)[0][0]
            
            all_2hr_sieved_times_product_times_inds_list.remove(ind_hole_time_val)

            os.remove(query_result[0]) #delete original downloaded file
            ind_timespickup = np.where(np.array(all_size_sieved_times) == hole_time_val)[0][0]
            zoomed_time_range = TimeRange(str(hole_time_val),timedelta(hours=time_window))

            fetch_inds_to_try_list = [] 
            #the zeroth entry didn't have it so that's why plus 1 in the brackets
            for time_val in all_size_sieved_times[ind_timespickup+1: ind_timespickup + look_ahead]:
                if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                    ind_next_good_time = np.where(np.array(all_size_sieved_times) == time_val)[0]
                    fetch_indices_next_good = ind[ind_next_good_time]
                    fetch_inds_to_try_list.append(fetch_indices_next_good)

            for index in fetch_inds_to_try_list:
                query_result_next = product_retriever(base,product_results,index,url_prefix,home_dir)
                axis1_next_good,axis2_next_good,data_next_good = readfits(query_result_next[0])

                if (data_next_good is not None) and (axis1_next_good == axis2_next_good):

                    if not holes(query_result_next[0]): #so if not True; so no holes; can use image
                        reduced_product_data = data_reducer(data_next_good,flag,target_dimension,axis1_next_good)
                        time_data = product_results.get_response(0)[int(index)]['time']['start']
                        writefits(f'{base}/SOHO_{base}_{time_data}_{target_dimension}', reduced_product_data, home_dir)

                        all_2hr_sieved_times_product_times.append(time_data) #(time_val) #unsorted time location
                        all_2hr_sieved_times_product_times_inds_list.append(index)
                        os.remove(query_result_next[0]) #delete original downloaded file
                        break

                    elif holes(query_result_next[0]): #so if True, if there are holes
                        time_data = product_results.get_response(0)[int(index)]['time']['start']
                        hole_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                        holes_product_list.append((hole_loc, str(time_data)))
                        os.remove(query_result_next[0])
                        continue 

                elif (data_next_good is None) or (axis1_next_good != axis2_next_good):
                    unreadable_file_ids_product_list.append(product_results.get_response(0)[int(index)]['fileid'])
                    os.remove(query_result_next[0])
                    continue #continue the for loop


    elif (data_product is None) or (axis1_product != axis2_product):
        unreadable_file_ids_product_list.append(product_results.get_response(0)[int(indiv_ind)]['fileid'])
        bad_time_val = product_results.get_response(0)[int(indiv_ind)]['time']['start']
        all_2hr_sieved_times_product_times.remove(bad_time_val)
        ind_bad_time_val = np.where(np.array(all_size_sieved_times) == bad_time_val)[0][0]
        all_2hr_sieved_times_product_times_inds_list.remove(ind_bad_time_val)
        os.remove(query_result[0]) #delete original downloaded file
        ind_timespickup = np.where(np.array(all_size_sieved_times) == bad_time_val)[0][0]
        zoomed_time_range = TimeRange(str(bad_time_val),timedelta(hours=time_window))

        fetch_inds_to_try_list = [] #gets reset for each new item
        for time_val in all_size_sieved_times[ind_timespickup+1: ind_timespickup + look_ahead]:
            if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                ind_next_good_time = np.where(np.array(all_size_sieved_times) == time_val)[0]
                fetch_indices_next_good = ind[ind_next_good_time]
                fetch_inds_to_try_list.append(fetch_indices_next_good)

        for index in fetch_inds_to_try_list:
            query_result_next = product_retriever(base,product_results,index,url_prefix,home_dir)
            axis1_next_good,axis2_next_good,data_next_good = readfits(query_result_next[0])

            if (data_next_good is not None) and (axis1_next_good == axis2_next_good):

                if not holes(query_result_next[0]): #so if not True; so no holes; can use image
                    reduced_product_data = data_reducer(data_next_good,flag,target_dimension,axis1_next_good)
                    time_data = product_results.get_response(0)[int(index)]['time']['start']
                    writefits(f'{base}/SOHO_{base}_{time_data}_{target_dimension}', reduced_product_data, home_dir)

                    all_2hr_sieved_times_product_times.append(time_data) #(time_val) #unsorted time location
                    all_2hr_sieved_times_product_times_inds_list.append(index)
                    os.remove(query_result_next[0])
                    break

                elif holes(query_result_next[0]): #so if True, if there are holes
                    time_data = product_results.get_response(0)[int(index)]['time']['start']
                    hole_loc = url_prefix + product_results.get_response(0)[int(index)]['fileid']
                    holes_product_list.append((hole_loc, str(time_data)))
                    os.remove(query_result_next[0])
                    continue 

            elif (data_next_good is None) or (axis1_product != axis2_product):
                unreadable_file_ids_product_list.append(product_results.get_response(0)[int(index)]['fileid'])
                os.remove(query_result_next[0])
                continue
    
    all_2hr_sieved_times_product_times_modified = all_2hr_sieved_times_product_times

    return all_2hr_sieved_times_product_times_modified, holes_product_list, unreadable_file_ids_product_list 
    #think whether need this new name or need to feed back, i think is ok


def csv_writer(base,home_dir,date_start,date_finish,flag,target_dimension, all_2hr_sieved_times_sorted):
    with open(f'{home_dir}{date_start}_to_{date_finish}_{base}_times_{flag}_{target_dimension}.csv', 'a') as f: #appending lines so not overwriting the file
        writer = csv.writer(f, delimiter='\n')
        writer.writerow(all_2hr_sieved_times_sorted)


def main(date_start, date_finish, target_dimension, time_increment, time_window, flag, home_dir):
    
    date_time_pre_start = date_start + '-0000'
    date_time_start= parser.parse(date_time_pre_start)
    print('date_time_start:', date_time_start)

    date_time_pre_end = date_finish + '-2359'
    date_time_end = parser.parse(date_time_pre_end)
    print('date_time_end:', date_time_end)

    target_dimension = int(target_dimension)
    print('target_dimension:', target_dimension)

    time_increment = int(time_increment)
    time_window = float(time_window)

    flag = str(flag) 
    print('flag:', flag)

    home_dir = str(home_dir)
    print('home_dir:', home_dir)

    url_prefix = 'https://seal.nascom.nasa.gov/'
    print('url_prefix:', url_prefix)

    look_ahead = int(np.ceil(time_window*60/10)) #should sufficiently cover all 6 products based on their cadence.
    print('look_ahead:', look_ahead)

    diff_start_finish_total_sec = (date_time_end - date_time_start).total_seconds()
    print('diff_start_finish_total_sec:', diff_start_finish_total_sec)

    total_sec = timedelta(days = time_increment).total_seconds()
    print('total_sec:', total_sec)

    num_loops = np.ceil(diff_start_finish_total_sec/total_sec) + 1 #num_loops would be equal to 94 + 1 for 19960101-0000' - '20110501-0000'; discete number of loops so go over rather than under
    print('num_loops:', num_loops)

    for base in ['EIT195', 'MDI_96m','LASCO_C2','LASCO_C3','EIT171','EIT304','EIT284']: #[1:2] can place range to run a subset of the products here
        holes_list = []
        unreadable_file_ids_product_list_global = []
    
        print(f'***{base}***')
        base_dir = home_dir + base
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)    

        time_range = TimeRange(date_time_start, timedelta(days = time_increment)) #time_range re-initialized here
        #print('time_range:', time_range)

        for t_value in np.arange(num_loops): #main workhorse loop
            print('t_value:', t_value)
                
            if time_range.end > date_time_end:
                time_range = TimeRange(time_range.start, date_time_end)  
                       
            product_results = product_search(base,time_range,date_time_start)
            product_results_number = product_results.file_num
            if product_results_number != 0:
                ind = index_of_sizes(base,product_results)
                all_size_sieved_times, all_2hr_sieved_times_product_times, all_2hr_sieved_times_product_times_inds_list, fetch_indices_product = fetch_indices(base,ind,product_results,time_window,look_ahead)
                for item in fetch_indices_product:
                    query_result = product_retriever(base,product_results,item,url_prefix,home_dir)
                    axis1_product,axis2_product,data_product = readfits(query_result[0])
                    all_2hr_sieved_times_product_times_modified, holes_product_list, unreadable_file_ids_product_list_local = product_distiller(base, axis1_product,axis2_product,data_product, all_size_sieved_times, all_2hr_sieved_times_product_times, all_2hr_sieved_times_product_times_inds_list, query_result, ind, item, product_results, look_ahead, time_window, url_prefix, flag, target_dimension, home_dir)
                
                    if holes_product_list:
                        holes_list.append(holes_product_list)
                    
                    if unreadable_file_ids_product_list_local:
                        unreadable_file_ids_product_list_global.append(unreadable_file_ids_product_list_local)
            
                all_2hr_sieved_times_sorted = np.sort(all_2hr_sieved_times_product_times_modified)

                print(f'{base} all_size_sieved_times:', all_size_sieved_times, len(all_size_sieved_times))
                print(f'{base} list(all_2hr_sieved_times_sorted):', list(all_2hr_sieved_times_sorted), len(all_2hr_sieved_times_sorted))
            
                csv_writer(base,home_dir,date_start,date_finish,flag,target_dimension, all_2hr_sieved_times_sorted)

            time_range.next() #Sunpy iterator to go for the next 2 months #also have time_range.previous() to go back. #### UNCOMMENT!    
            #print('time_range next:', time_range)
        
        print(f'{base} holes_list', holes_list)
        print(f'{base} unreadable_file_ids_product_list_global:', unreadable_file_ids_product_list_global)

        data_cuber(home_dir, base, date_start, date_finish, flag, target_dimension)    
    
    

    
if __name__ == '__main__':
    import argparse
    parser_args = argparse.ArgumentParser(description='SOHO ML Data experiment parameters')
    parser_args.add_argument('--date_start', metavar='time', required=True, help='yyyy-mm-dd, 1996-01-01 is earliest start', type = str)
    parser_args.add_argument('--date_finish', metavar='time', required=True, help='yyyy-mm-dd, 2011-05-01 is recommended latest finish, select a minimum range of 2 months', type = str)
    parser_args.add_argument('--target_dimension', metavar='image size', required=True, help='e.g., 128', type = int)
    parser_args.add_argument('--time_increment', metavar='days at a time to loop over', required=False, help='max time span must be around 2 months as there is a 10k limit to VSO return search query', default = 60, type = int)
    parser_args.add_argument('--time_window', metavar='time', required=True, help='time step in hours', type = float)
    parser_args.add_argument('--flag', metavar='resize strategy', required=True, help='choose from either "subsample", "interp", "minpool", or "maxpool" ', type = str)
    parser_args.add_argument('--home_dir', metavar='home directory', required=True, help='str, e.g., "/home/user/Documents/", need "/" in the end', type = str)


    args = parser_args.parse_args()
    main(
        date_start = args.date_start,
        date_finish = args.date_finish,
        target_dimension = args.target_dimension,
        time_increment = args.time_increment,
        time_window = args.time_window,
        flag = args.flag,
        home_dir = args.home_dir)  
