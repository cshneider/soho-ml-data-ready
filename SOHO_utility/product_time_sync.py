import numpy as np

import os, fnmatch
from os import listdir
from os.path import isfile, join

from datetime import datetime, date, time, timedelta
from dateutil import parser

from sunpy.time import TimeRange

import h5py
import csv
from tqdm import tqdm

"""
FIND SPECIFIED PATTERN USING FNMATCH 
"""
def pattern_finder(home_dir, pattern):

    for root, dirs, files in os.walk(home_dir):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                name = str(name)
                break
        break
    
    return name

"""
READ IN TIMES FROM FITS FILES PER PRODUCT TYPE IF THESE HAVE BEEN GENERATED LOCALLY BY PREVIOUSLY RUNNING SOHO_DATA_GEN.PY.
FITS FILES ARE ALL UNIQUE AND SORTED TO ENSURE THAT THEY FOLLOW THE ORDER OF THE H5 DATA CUBE THAT HAS BEEN MADE EARLIER BY RUNNING SOHO_DATA_GEN.PY.
"""
def fits_times_reader(home_dir, base, mission): 

    print('base:', base)
    filepath = home_dir + base  + f'_{mission}' + '/'

    data_files_pre_pre = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    data_files_pre = [f for f in data_files_pre_pre if 'fits' in f] #to ensure that only FITS files are collected, just in case    
    data_files = np.sort(data_files_pre)
    
    if 'EIT' in base:
        data_raw_times = [elem.split('_')[2] for elem in data_files]
    else:
        data_raw_times = [elem.split('_')[3] for elem in data_files]            
        
    return data_raw_times
    
"""
READ IN TIMES FROM CSV FILES PER PRODUCT TYPE. TAKING UNIQUE VALUES TO ENSURE UNIQUE TIMES IN CSV FILES. {THIS IS THE CASE EXCEPT FOR LASCO_C2 WHERE THERE APPEARS ONCE TIME DUPLICATE}
"""
def csv_times_reader(home_dir, pattern):

    name = pattern_finder(home_dir, pattern)
    print('name from csv_times_reader:', name)
    with open(f'{home_dir}{name}', 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\n')
        csv_data = [line for line in csv_reader]
    csv_uniq_times = list(np.unique(csv_data))
    print('len(csv_uniq_times):', len(csv_uniq_times))
                
    return csv_uniq_times

"""
CHECKS THAT THE DIMENSION AMONG FITS FILES COMING FROM THE DIFFERENT SPECIFIED PRODUCTS IS INDEED THE SAME.
"""
def dimension_checker_from_fits(home_dir, base_list, mission):

    data_dim_list = []
    for base in base_list:
        base = base.strip(' ')
        filepath = home_dir + base + f'_{mission}' + '/'
        data_file_pre = [f for f in listdir(filepath) if isfile(join(filepath, f))] #[0]
        data_file = [f for f in data_file_pre if 'fits' in f][0]
    
        if 'EIT' in base:
            data_dim = data_file.split('_')[3].split('.')[0]
        else:
            data_dim = data_file.split('_')[4].split('.')[0]        
    
        data_dim_list.append(data_dim)
    
    ind_dim = np.where(data_dim_list[0] == np.array(data_dim_list))[0] #seems that np.array() is required since have str as elements.
    
    if len(ind_dim) == len(base_list):
        return True
    else:
        return False
        
"""
CHECKS THAT THE DIMENSION AMONG THE H5 CUBE AND CSV FILES COMING FROM THE DIFFERENT SPECIFIED PRODUCTS IS INDEED THE SAME.
"""
def dimension_checker_from_h5cube_csv(home_dir, base_list, mission): #assumes all h5 files in the same directory which is one above the product (base) directory

    data_dim_list = []
    for base in base_list:
        base = base.strip(' ')
        
        name = pattern_finder(home_dir, pattern = f'*{base}*{mission}*[!sync].h5')
        print('cube name:', name)
        cube_dim = name.split('_')[-1].split('.')[0]
        print('cube_dim:', cube_dim)
        data_dim_list.append(cube_dim)

    ind_dim = np.where(data_dim_list[0] == np.array(data_dim_list))[0] 
     
    if len(ind_dim) == len(data_dim_list):
        return True
    else: 
        return False        

"""
CONVERTS TIME STRINGS INTO DATETIME OBJECTS. ALLOWS TO SPECIFY A SUBSET OF THE DATE RANGES THAT HAD BEEN USED WHEN RUNNING SOHO_DATA_GEN.PY. 
"""
def times_actualizer(data_raw_times, date_start, date_finish): #produces a subset of the original times if so desired by the user. Adds flexibility in choice of time range!

    date_time_pre_start = date_start + '-0000'
    date_time_start= parser.parse(date_time_pre_start)
    print('date_time_start:', date_time_start)

    date_time_pre_end = date_finish + '-2359'
    date_time_end = parser.parse(date_time_pre_end)
    print('date_time_end:', date_time_end)
    
    data_times = np.array([parser.parse(elem) for elem in data_raw_times]) #this needs to be an array for the where fcn right below not to require array explicitely there!
    
    ind_start = np.where(data_times >= date_time_start)[0] #shifted the second [0] below just in case dates were selected incorrectly
    ind_end = np.where(data_times <= date_time_end)[0] #shifted the second [-1] below just in case dates were selected incorrectly
    
    if (len(ind_start) != 0) and (len(ind_end) != 0): #so very important to have -1 and not 0 here since ind_end[-1] corresponds to the last date
        data_times_revised = data_times[ind_start[0]:ind_end[-1]+1] #a flexible range offered here #+1 since bracket range is one less. #this corresponds to cube_data[ind_start:ind_end+1] in cube_sync_maker.
        print('len(data_times_revised):', len(data_times_revised))
    
    else:
        data_times_revised = data_times.copy()
        raise ValueError("date selected is outside original date start and date finish range") 
        
    return list(data_times_revised), data_times, ind_start, ind_end

"""
FINDS THE TIME_WINDOW USED PREVIOUSLY WHEN RUNNING SOHO_DATA_GEN.PY. IN THIS PROGRAM THE NAME IS CHANGED TO TIME_STEP.
"""
def min_time_step(data_times):
    
    data_times_diff_pre = data_times[1:] - data_times[:-1]
    data_times_diff = [elem for elem in data_times_diff_pre] #.total_seconds()
    
    min_time_diff = np.min(data_times_diff)

    return min_time_diff

"""
READS THE TIME WINDOW USED WHEN RUNNING SOHO_DATA_GEN.PY FROM THE H5 DATA FILE. THIS TIME_WINDOW IS TIME_STEP_PREV HERE.
"""
def time_step_prev_reader(home_dir, pattern):

    name = pattern_finder(home_dir, pattern)
    print('name from time_step_prev_reader:', name)
    time_step_prev = name.split('_')[-4] #was [-2] previously
    print('time_step_prev from fcn:', time_step_prev)
    
    return int(time_step_prev)
    
"""
FIND CORRESPONDING DATA CUBE PER PRODUCT AND EXTRACT ITS DATA AND DIMENSION.
"""
def cube_data_reader(home_dir, base, mission, pattern):

    name = pattern_finder(home_dir, pattern)            
    print('cube name:', name)
    cube_dim = name.split('_')[-1].split('.')[0] #str
    print('cube_dim:', cube_dim)
                
    cube = h5py.File(f'{home_dir}{name}', 'r')
    cube_data = cube[f'{base}_{mission}_{cube_dim}'][:]
    cube.close()

    return cube_data, cube_dim
    

"""
START WITH SHORTEST PRODUCT TIME LIST SINCE THAT'S THE NATURAL BOTTELNECK.
"""
def shortest_prod_list_index_finder(product_list):

     product_list_lengths = [len(i) for i in product_list]
     print('product_list_lengths:', product_list_lengths)
     product_list_lengths_min = np.min(product_list_lengths)
     ind_min_len = np.where(np.array(product_list_lengths) == product_list_lengths_min)[0][0]
     print('product_list_lengths[ind_min_len]:', product_list_lengths[ind_min_len])
     
     return ind_min_len


"""
FINDING THE TIMES AND INDICES THAT ARE SYNCHED BETWEEN THE ENTERED PRODUCTS. MOVES ALONG THE SHORTEST PRODUCT LIST. OVERLAP TIME INTERVAL DETERMINED BY MOVING +/- HALF OF THE ORIGINAL TIME_STEP.
THEN IF THE TIME STEP ENTERED IN THIS CODE IS AN INTEGER MULTIPLE OF THE ORIGINAL TIME STEP, THE SYNCHED TIMES COMPUTED ON THE ORIGINAL TIME STEP ARE SUBSEQUENTLY SUBSAMPLED.
"""
def sync_times_and_inds(product_list, ind_min_len, time_step, time_step_prev): #main engine of the algorithm

    synch_time_list = []
    synch_time_inds_list = []
    
    ratio = int(time_step / time_step_prev) #in order to subsample (e.g., time_step_prev=6 but now want time_step=12, so will need to take every other element found from synching by original time)
    print('ratio:', ratio)
    
    for i,time_val in tqdm(enumerate(product_list[ind_min_len])): #so moving along shortest product list
        time_range = TimeRange(time_val - timedelta(hours=time_step_prev/2), time_val + timedelta(hours=time_step_prev/2)) #heart of the algorithm
        #always synch on original time_window and then can subsample from that if necessary. #time_val is a datetime object so can subtract a timedelta from it.
        temp_ind_list = []
        temp_time_list = []
        for j,product in enumerate(product_list):
            time_range_list = [(item in time_range) for item in product] #boolean #bool
            if any(np.array(time_range_list)): #to ensure that have at least one entry that is True.
                ind_temp = np.where(time_range_list)[0][0] #this second [0] ensures that only the first matching time is taken and hence that the resulting times sieved are uniquely obtained!
                temp_time_pre = product[ind_temp]
                temp_time = ''.join(str(temp_time_pre).split(' ')[0].split('-')) + ''.join(str(temp_time_pre).split(' ')[1].split(':'))
            else:
                ind_temp = np.nan
                temp_time = np.nan
            
            temp_ind_list.append(ind_temp)
            temp_time_list.append(temp_time)
        
        if len(np.where(np.array(temp_ind_list) != np.array(temp_ind_list))[0]) == 0: #this is used to pick only those instances where nan vals are not present.
            synch_time_inds_list.append(temp_ind_list) #is of len(base_list)            
            synch_time_list.append(temp_time_list) #is of len(base_list)     

    return synch_time_inds_list[::ratio], synch_time_list[::ratio]


"""
USE THE TIMES AND INDICES THAT ARE SYNCHED BETWEEN THE ENTERED PRODUCTS AND REORDER BY ORDER OF PRODUCTS FOLLOWING ORDER OF BASE_LIST: USER ENTERED ORDER OF PRODUCTS
"""
def sync_times_and_inds_sort_by_product(synch_time_inds_list, synch_time_list):

     synch_time_inds_list_ravel = np.ravel(synch_time_inds_list, order='F')
     synch_time_inds_list_mod = np.hsplit(synch_time_inds_list_ravel,len(synch_time_inds_list[0])) #len(synch_time_inds_list[0]) should be equal to len(base_list)!
     
     synch_time_list_ravel = np.ravel(synch_time_list, order='F')
     synch_time_list_mod = np.hsplit(synch_time_list_ravel,len(synch_time_list[0])) #len(synch_time_list[0]) should be equal to len(base_list)!
     
     return synch_time_inds_list_mod, synch_time_list_mod

"""
OUTPUTS TIMES AND INDICES TO BE USED FROM LASCO DIFFERENCE IMAGES. ### NEED TO TAKE INTO ACCOUNT THAT COULD HAVE BOTH C2 AND C3 PRESENT SIMULTANEOUSLY!!!
"""
def lasco_diff_times_inds(lasco_sync_times):
     
     #print('lasco_sync_times_internal:', lasco_sync_times)
     synced_lasco_datetimes = [parser.parse(elem) for elem in lasco_sync_times]
     #print('synced_lasco_datetimes:', synced_lasco_datetimes)
     synced_lasco_datetimes_Fcorona_remov_pre = np.array(synced_lasco_datetimes[1:]) - np.array(synced_lasco_datetimes[:-1])
     #print('synced_lasco_datetimes_Fcorona_remov_pre:', synced_lasco_datetimes_Fcorona_remov_pre)
     synced_lasco_datetimes_Fcorona_remov = [np.round(elem.total_seconds()/3600.) for elem in synced_lasco_datetimes_Fcorona_remov_pre]
     #print('synced_lasco_datetimes_Fcorona_remov:', synced_lasco_datetimes_Fcorona_remov)
     lasco_ind_Fcorona_24h = np.where(np.array(synced_lasco_datetimes_Fcorona_remov) <= 24)[0]
     #print('lasco_ind_Fcorona_24h:', lasco_ind_Fcorona_24h)
     print('len(lasco_ind_Fcorona_24h):', len(lasco_ind_Fcorona_24h))     
     
     return lasco_ind_Fcorona_24h

"""
OUTPUTS DATA CUBES FOR EACH SPECIFIED PRODUCT. THESE CUBES ARE THE REDUCED VERSIONS OF THE ORIGINAL ONES SINCE ONLY THE TIME SLICES THAT COME WITHIN THE SPECIFIED TIME_STEP HAVE BEEN RETAINED.
"""
def cube_sync_maker(home_dir, base, base_list_len, cube_data, cube_dim, ind_start, ind_end, synch_time_inds_mod, date_start, date_finish, time_step_prev, time_step, mission, flag_lasco=None):

     if flag_lasco is None:
          cube_data_mod_pre_pre = cube_data[ind_start:ind_end+1]
          cube_data_mod_pre = np.array([cube_data_mod_pre_pre[i] for i in synch_time_inds_mod])
          cube_data_mod = cube_data_mod_pre.astype('int16')
               
          cube_sync = h5py.File(f'{home_dir}{date_start}_to_{date_finish}_{base}_{mission}_{base_list_len}products_{time_step_prev}_{time_step}_{cube_dim}_sync.h5', 'w')
          cube_sync.create_dataset(f'{base}_{mission}_{cube_dim}', data=cube_data_mod) #not compressing images here since images compressed initially in data generation step #compression="gzip"
          cube_sync.close()
     
     else:
          cube_data_mod = cube_data.astype('int16')
          
          cube_sync = h5py.File(f'{home_dir}{date_start}_to_{date_finish}_{base}_{flag_lasco}_{mission}_{base_list_len}products_{time_step_prev}_{time_step}_{cube_dim}_sync.h5', 'w')
          cube_sync.create_dataset(f'{base}_{mission}_{cube_dim}', data=cube_data_mod) #not compressing images here since images compressed initially in data generation step #compression="gzip"
          cube_sync.close()     
     
     return cube_sync
     
"""
OUTPUTS A CSV FILE CONTAINING THE RETAINED TIMES PER SPECIFIED PRODUCT WHICH COINCIDE WITH THE TIMES OF OTHER PRODUCTS WITHIN THE TIME_STEP.
"""
def csv_time_sync_writer(home_dir, base, base_list_len, date_start, date_finish, cube_dim, synch_time_list_mod, time_step_prev, time_step, mission, flag_lasco=None):
     
     if flag_lasco is None:    
          if not isfile(f'{home_dir}{date_start}_to_{date_finish}_{base}_{mission}_{base_list_len}products_{time_step_prev}_{time_step}_{cube_dim}_times_sync.csv'):
               with open(f'{home_dir}{date_start}_to_{date_finish}_{base}_{mission}_{base_list_len}products_{time_step_prev}_{time_step}_{cube_dim}_times_sync.csv', 'a') as f:
                 writer = csv.writer(f, delimiter='\n')
                 writer.writerow(synch_time_list_mod)
     else:
          if not isfile(f'{home_dir}{date_start}_to_{date_finish}_{base}_{mission}_{base_list_len}products_{flag_lasco}_{time_step_prev}_{time_step}_{cube_dim}_times_sync.csv'):
               with open(f'{home_dir}{date_start}_to_{date_finish}_{base}_{mission}_{base_list_len}products_{flag_lasco}_{time_step_prev}_{time_step}_{cube_dim}_times_sync.csv', 'a') as f:
                 writer = csv.writer(f, delimiter='\n')
                 writer.writerow(synch_time_list_mod)     
     
  
  
     
