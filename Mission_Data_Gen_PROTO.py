import numpy as np
import os
from dateutil import parser
from datetime import timedelta
from sunpy.time import TimeRange
from time import process_time
from tqdm import tqdm
import drms

import Mission_utility.data_gen_helper as ipy
import Mission_utility.sdo_mdi as sdo
import Mission_utility.soho_other as soho

date_start =  '2001-04-01'#'2002-04-01' #'1999-04-04' # '2010-10-01' '2005-12-01'
date_finish = '2001-04-07'#'2002-04-07' #'1999-04-09' # '2010-10-05' '2005-12-07'
image_size_output = 128
time_window = 6
flag = 'subsample'
home_dir = '/Users/gohawks/Desktop/soho-ml-data/soho-ml-data-ready-martinkus/'
bases = 'MDI_96m' #,LASCO_C3'#MDI_96m #AIA #HMI #'EIT195'
fits_headers = 'N'
lev1_LASCO = 'Y' #CANNOT USE 'N' UNTIL UNIT CONVERSION SORTED OUT
email = 'charlotte.martinkus@noaa.gov'


#def main(date_start, date_finish, image_size_output, time_window, flag, home_dir, bases, fits_headers, lev1_LASCO, email):
    
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

url_prefix = 'https://seal.nascom.nasa.gov'
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
        if ('LASCO' in base) or ('EIT' in base):
            BaseClass = soho.SOHO_no_MDI(base, lev1_LASCO, fits_headers)
            
        else:
            BaseClass = sdo.SDO_MDI(base, lev1_LASCO, fits_headers, time_window)
            
        BaseClass.set_base_dictionary()
        
        holes_list = []
        transients_list = []
        unreadable_file_ids_product_list_global = []
    
        print(f'{base}')
        base_dir = home_dir + base + '_' + BaseClass.mission
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)   
            
        time_range = TimeRange(date_time_start, timedelta(days = time_increment))
        
        
        prev_time, time_range_modified = ipy.prev_time_resumer(home_dir, time_range, date_time_end, BaseClass) #time_range_modified.next() is the workshorse that advances time at the end of the time for-loop
        for t_value in tqdm(np.arange(num_loops)): #this is the main workhorse loop of the program
            #print('t_value:', t_value)
            print('prev_time:', prev_time)
                
            if time_range_modified.end > date_time_end:
                time_range_modified = TimeRange(time_range_modified.start, date_time_end) 
                
            product_results, client = BaseClass.product_search(time_range_modified, client)
            
            if BaseClass.class_type == 'SDO_MDI':
                client_export_failed = product_results.has_failed()
                
                if client_export_failed == True:
                    product_results_number = 0
                
                elif client_export_failed == False:
                    product_results_number = product_results.data.count()['record'] 
            
            else:
                product_results_number = product_results.file_num
                client_export_failed = False
                
            if (product_results_number != 0) and (client_export_failed == False): #product_results.has_failed()
                ind = BaseClass.index_of_sizes(product_results) 
                size_sieved_df, fetch_indices_product_orig = ipy.fetch_indices(ind,product_results,time_window, prev_time, time_range_modified, BaseClass)
                if len(fetch_indices_product_orig) != 0:
                
                    size_sieved_df = ipy.product_distiller(fetch_indices_product_orig, size_sieved_df, date_time_end, product_results, look_ahead, time_window, url_prefix, flag, image_size_output, home_dir, email, client, BaseClass)
                    
                    all_time_window_sieved_times_sorted = np.unique(list(size_sieved_df['time_at_ind']))
                                
                    prev_time = [] #reset to empty list
                    if len(all_time_window_sieved_times_sorted) != 0:
                        prev_time.append(all_time_window_sieved_times_sorted[-1]) #append the last good time entry from the previous loop     
                    size_sieved_df.to_csv(f'{home_dir}{date_start}_to_{date_finish}_{BaseClass.base_full}_times_{flag}_{time_window}_LASCOlev1-{BaseClass.lev1_lasco}_{BaseClass.mission}_{image_size_output}_{t_value}.csv')
                        
                        #ipy.csv_writer(home_dir, date_start, date_finish, flag, time_window, image_size_output, all_time_window_sieved_times_sorted, BaseClass)
                    
    
                time_range_modified.next() #Sunpy iterator to go for the next time increment in number of days. There is also time_range_modified.previous() to go backwards in time.    
                #print('time_range_modified next:', time_range_modified)
            


        ipy.data_cuber(home_dir, date_start, date_finish, flag, time_window, image_size_output,BaseClass)
        
        end_process_time = process_time()
        time_of_process = end_process_time - start_process_time
        print(f'{base} time of process in seconds:', time_of_process)      
    
    
    
    
 
    
 
# if __name__ == '__main__':
#     import argparse
#     parser_args = argparse.ArgumentParser(description='Mission ML Data experiment parameters')
#     parser_args.add_argument('--date_start', metavar='time', required=True, help='yyyy-mm-dd, 1996-01-01 is earliest start', type = str)
#     parser_args.add_argument('--date_finish', metavar='time', required=True, help='yyyy-mm-dd, 2011-05-01 is recommended latest finish for SOHO mission, select a max range of 2 months', type = str)
#     parser_args.add_argument('--image_size_output', metavar='image size', required=True, help='e.g., 128', type = int)
#     parser_args.add_argument('--time_window', metavar='time', required=True, help='time step in hours', type = int)
#     parser_args.add_argument('--flag', metavar='resize strategy', required=True, help='choose from either "subsample", "interp", "minpool", or "maxpool" ', type = str)
#     parser_args.add_argument('--home_dir', metavar='home directory', required=True, help='str, e.g., "/home/user/Documents/", need "/" in the end', type = str)
#     parser_args.add_argument('--products', metavar='product types', required=True, help='str, Enter all the following or a subset thereof, in any order, seperated by commas: for the SOHO mission: "EIT195, MDI_96m, LASCO_C2, LASCO_C3, EIT171, EIT304, EIT284" and for the SDO mission: "HMI_720s, AIA94, AIA131, AIA171, AIA193, AIA211, AIA304, AIA335, AIA1600, AIA1700, AIA4500"', type = str)
#     parser_args.add_argument('--fits_headers', metavar='FITS header metadata',required=True, help='Include header metadata in individual FITS files? Y/y or N/n', type = str)
#     parser_args.add_argument('--lev1_LASCO', metavar='level 1.0 or level 0.5 from SDAC', required=True, help='If want level-1 LASCO C2 and C3 images? Y/y or N/n', type = str)
#     parser_args.add_argument('--email', metavar='email address of user', required=True, help='Required to obtain calibrated MDI images from JSOC', type = str)
    
#     args = parser_args.parse_args()
#     main(
#         date_start = args.date_start,
#         date_finish = args.date_finish,
#         image_size_output = args.image_size_output,
#         time_window = args.time_window,
#         flag = args.flag,
#         home_dir = args.home_dir,
#         bases = args.products,
#         fits_headers = args.fits_headers,
#         lev1_LASCO = args.lev1_LASCO,
#         email = args.email)      
 
    
 
    
 
    
 
    

