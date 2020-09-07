import numpy as np
import os
from dateutil import parser
from datetime import datetime, date, time, timedelta
from sunpy.time import TimeRange

from SOHO_utility import *

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

        prev_time = []
        for t_value in np.arange(num_loops): #main workhorse loop
            print('t_value:', t_value)
            print('prev_time:', prev_time)
                
            if time_range.end > date_time_end:
                time_range = TimeRange(time_range.start, date_time_end)  
                       
            product_results = product_search(base,time_range,date_time_start)
            product_results_number = product_results.file_num
            if product_results_number != 0:
                ind = index_of_sizes(base,product_results)
                all_size_sieved_times, all_2hr_sieved_times_product_times, all_2hr_sieved_times_product_times_inds_list, fetch_indices_product = fetch_indices(base,ind,product_results,time_window,look_ahead, prev_time)
                for item in fetch_indices_product:
                    query_result = product_retriever(base,product_results,item,url_prefix,home_dir)
                    axis1_product,axis2_product,data_product = readfits(query_result[0])
                    all_2hr_sieved_times_product_times_modified, holes_product_list, unreadable_file_ids_product_list_local = product_distiller(base, axis1_product,axis2_product,data_product, all_size_sieved_times, all_2hr_sieved_times_product_times, all_2hr_sieved_times_product_times_inds_list, query_result, ind, item, product_results, look_ahead, time_window, url_prefix, flag, target_dimension, home_dir)
                
                    if holes_product_list:
                        holes_list.append(holes_product_list)
                    
                    if unreadable_file_ids_product_list_local:
                        unreadable_file_ids_product_list_global.append(unreadable_file_ids_product_list_local)
            
                if not prev_time:
                    all_2hr_sieved_times_sorted = np.unique(all_2hr_sieved_times_product_times_modified)
                elif prev_time:                    
                    all_2hr_sieved_times_sorted = np.unique(all_2hr_sieved_times_product_times_modified[1:])
                                
                print(f'{base} all_size_sieved_times:', all_size_sieved_times, len(all_size_sieved_times))
                print(f'{base} list(all_2hr_sieved_times_sorted):', list(all_2hr_sieved_times_sorted), len(all_2hr_sieved_times_sorted))

                prev_time = [] #reset to empty list
                if len(all_2hr_sieved_times_sorted) != 0:
                    prev_time.append(all_2hr_sieved_times_sorted[-1]) #append the last good time entry from the previous loop
                            
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
