import numpy as np
import os
import warnings
from dateutil import parser
from datetime import datetime, date, time, timedelta
from sunpy.time import TimeRange
from time import process_time
from tqdm import tqdm
from Mission_utility import *


def main(date_start, date_finish, image_size_output, time_window, flag, home_dir, bases, fits_headers, lev1_LASCO, email, mission):
    
    warnings.simplefilter(action='ignore', category=RuntimeWarning)
    
    client = drms.Client(email=email, verbose=False) #True
    
    date_time_pre_start = date_start + '-0000'
    date_time_start= parser.parse(date_time_pre_start)
    print('date_time_start:', date_time_start)

    date_time_pre_end = date_finish + '-2359'
    date_time_end = parser.parse(date_time_pre_end)
    print('date_time_end:', date_time_end)

    time_increment = 60
    print(f'time increment {time_increment} days')
    
    print('image_size_output:', image_size_output)
    print('flag:', flag)
    print('home_dir:', home_dir)

    url_prefix = 'https://seal.nascom.nasa.gov/'
    print('url_prefix:', url_prefix)

    look_ahead = int(np.ceil(time_window*60/10.)) #should sufficiently cover all 7 products based on their cadence.
    print('look_ahead:', look_ahead)

    diff_start_finish_total_sec = (date_time_end - date_time_start).total_seconds()
    print('diff_start_finish_total_sec:', diff_start_finish_total_sec)

    total_sec = timedelta(days = time_increment).total_seconds()
    print('total_sec:', total_sec)

    num_loops = np.ceil(diff_start_finish_total_sec/total_sec) + 1 #num_loops would be equal to 94 + 1 for 19960101-0000' - '20110501-0000'
    print('num_loops:', num_loops)

    base_list = bases.split(',')
    for base in tqdm(base_list):
        start_process_time = process_time() #initialize clock per product type
            
        base = base.strip(' ')
        
        holes_list = []
        transients_list = []
        blobs_list = []
        unreadable_file_ids_product_list_global = []
    
        print(f'{base}')
        base_dir = home_dir + base + f'_{mission}'
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)    

        if ('AIA' in base) and (int(base.split('AIA')[1]) != 4500): #AIA4500 is once an hour. So if not AIA4500, then have AIA cadences at 12s and 24s -> time_increment has to be reduced from 60d -> 1d
            time_increment = 1
            print(f'new time increment {time_increment} days:')
            total_sec = timedelta(days = time_increment).total_seconds()
            print('new total_sec:', total_sec)
            num_loops = np.ceil(diff_start_finish_total_sec/total_sec) + 1 #num_loops would be equal to 94 + 1 for 19960101-0000' - '20110501-0000'
            print('new num_loops:', num_loops)            
            time_range = TimeRange(date_time_start, timedelta(days = time_increment))
            
        else:
            time_range = TimeRange(date_time_start, timedelta(days = time_increment)) #time_range re-initialized here for each new base name 

        prev_time, time_range_modified = prev_time_resumer(home_dir, base, time_range, date_time_end, mission) #time_range_modified.next() is the workshorse that advances time at the end of the time for-loop
        for t_value in tqdm(np.arange(num_loops)): #this is the main workhorse loop of the program
            print('t_value:', t_value)
            print('prev_time:', prev_time)
                
            if time_range_modified.end > date_time_end:
                time_range_modified = TimeRange(time_range_modified.start, date_time_end)  
                       
            product_results, client = product_search(base, time_range_modified, client, mission) 
            
            if ('MDI' in base) or (mission == 'SDO'):
                client_export_failed = ghost_file_check(product_results)
                
                if client_export_failed == True:
                    product_results_number = 0
                
                elif client_export_failed == False:
                    product_results_number = product_results.data.count()['record'] #### old when this object was a query and not an export: product_results.count()['T_REC'] 
            
            else:
                product_results_number = product_results.file_num
                client_export_failed = False
                
            if (product_results_number != 0) and (client_export_failed == False): #product_results.has_failed()
                ind, fits_headers = index_of_sizes(base,product_results, fits_headers, lev1_LASCO, mission) #putting fits_headers here to insert it into __init__.py
                all_size_sieved_times_pre, fetch_indices_product_orig = fetch_indices(base,ind,product_results,time_window,look_ahead, prev_time, mission)
                if len(fetch_indices_product_orig) != 0:
                
                    all_time_window_sieved_times_product_times_modified, holes_product_list, transients_product_list, blobs_product_list, unreadable_file_ids_product_list_local = product_distiller(fetch_indices_product_orig, base, all_size_sieved_times_pre, ind, product_results, look_ahead, time_window, url_prefix, flag, image_size_output, home_dir, email, fits_headers, lev1_LASCO, client, mission)
                    
                    if holes_product_list: #if image had missing regions (e.g., arising from telemetry errors)
                        holes_list.append(holes_product_list)
                    
                    if transients_product_list:
                        transients_list.append(transients_product_list)
                    
                    if blobs_product_list:
                        blobs_list.append(blobs_product_list)
                                            
                    if unreadable_file_ids_product_list_local:
                        unreadable_file_ids_product_list_global.append(unreadable_file_ids_product_list_local)
            
                    all_time_window_sieved_times_sorted = np.unique(all_time_window_sieved_times_product_times_modified)
                                
                    #print(f'{base} np.unique(all_size_sieved_times_pre):', np.unique(all_size_sieved_times_pre), len(np.unique(all_size_sieved_times_pre)))
                    #print(f'{base} list(all_time_window_sieved_times_sorted):', list(all_time_window_sieved_times_sorted), len(all_time_window_sieved_times_sorted))

                    prev_time = [] #reset to empty list
                    if len(all_time_window_sieved_times_sorted) != 0:
                        prev_time.append(all_time_window_sieved_times_sorted[-1]) #append the last good time entry from the previous loop
                            
                    csv_writer(base,home_dir, date_start, date_finish, flag, time_window, image_size_output, all_time_window_sieved_times_sorted, lev1_LASCO, mission)
                

            time_range_modified.next() #Sunpy iterator to go for the next time increment in number of days. There is also time_range_modified.previous() to go backwards in time.    
            #print('time_range_modified next:', time_range_modified)
        
        print(f'{base} holes_list', holes_list)
        print(f'{base} transients_list', transients_list)
        print(f'{base} blobs_list', blobs_list)        
        print(f'{base} unreadable_file_ids_product_list_global:', unreadable_file_ids_product_list_global)

        data_cuber(home_dir, base, date_start, date_finish, flag, time_window, image_size_output, lev1_LASCO, mission)
        
        end_process_time = process_time()
        time_of_process = end_process_time - start_process_time
        print(f'{base} time of process in seconds:', time_of_process)            

    
if __name__ == '__main__':
    import argparse
    parser_args = argparse.ArgumentParser(description='SOHO ML Data experiment parameters')
    parser_args.add_argument('--date_start', metavar='time', required=True, help='yyyy-mm-dd, 1996-01-01 is earliest start', type = str)
    parser_args.add_argument('--date_finish', metavar='time', required=True, help='yyyy-mm-dd, 2011-05-01 is recommended latest finish for SOHO mission, select a max range of 2 months', type = str)
    parser_args.add_argument('--image_size_output', metavar='image size', required=True, help='e.g., 128', type = int)
    parser_args.add_argument('--time_window', metavar='time', required=True, help='time step in hours', type = int)
    parser_args.add_argument('--flag', metavar='resize strategy', required=True, help='choose from either "subsample", "interp", "minpool", or "maxpool" ', type = str)
    parser_args.add_argument('--home_dir', metavar='home directory', required=True, help='str, e.g., "/home/user/Documents/", need "/" in the end', type = str)
    parser_args.add_argument('--products', metavar='product types', required=True, help='str, Enter all the following or a subset thereof, in any order, seperated by commas: for the SOHO mission: "EIT195, MDI_96m, LASCO_C2, LASCO_C3, EIT171, EIT304, EIT284" and for the SDO mission: "HMI_720s, AIA94, AIA131, AIA171, AIA193, AIA211, AIA304, AIA335, AIA1600, AIA1700, AIA4500"', type = str)
    parser_args.add_argument('--fits_headers', metavar='FITS header metadata',required=True, help='Include header metadata in individual FITS files? Y/y or N/n', type = str)
    parser_args.add_argument('--lev1_LASCO', metavar='level 1.0 or level 0.5 from SDAC', required=True, help='If want level-1 LASCO C2 and C3 images? Y/y or N/n', type = str)
    parser_args.add_argument('--email', metavar='email address of user', required=True, help='Required to obtain calibrated MDI images from JSOC', type = str)    
    parser_args.add_argument('--mission', metavar='SOHO or SDO', required=True, help='SOHO or SDO mission', type = str)    
    
    args = parser_args.parse_args()
    main(
        date_start = args.date_start,
        date_finish = args.date_finish,
        image_size_output = args.image_size_output,
        time_window = args.time_window,
        flag = args.flag,
        home_dir = args.home_dir,
        bases = args.products,
        fits_headers = args.fits_headers,
        lev1_LASCO = args.lev1_LASCO,
        email = args.email,
        mission = args.mission)  
