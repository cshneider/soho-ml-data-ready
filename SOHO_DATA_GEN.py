import numpy as np

import os
import glob

import time
from time import clock

import scipy
from scipy.stats import skew, kurtosis
from scipy import interpolate

from skimage.transform import rescale
from skimage.measure import block_reduce

from datetime import datetime, date, time, timedelta
from dateutil import parser

from sunpy.net import Fido
from sunpy.net.vso import attrs as avso
from sunpy.time import parse_time, TimeRange

import astropy.units as u
from astropy.io import fits

from matplotlib import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import h5py

######## ADD seven_stats_calc and FRACTAL DIMENSION CALCULATION TOO AND NAIVE BAYES STUFF TOO --> put into the second companion code once have generated the data cubes! ########

def readfits(filename):
    ft = fits.open(filename, memmap=False)
    hdr = ft[0].header
    data = ft[0].data
    axis1 = hdr['naxis1']
    axis2 = hdr['naxis2']
    ft.close()
    return axis1,axis2,data

def writefits(filename, data):
    if not os.path.exists(f'{home_dir}{filename}.fits'):
        fitsname = fits.PrimaryHDU(data)
        fitsname.writeto(f'{home_dir}{filename}.fits')
	
def holes(filename): #change --> to hole_finder
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
    
    matches = ['96m', 'MDI', 'LASCO_C3']
    
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
        rad = float(x_coord)
        indices = np.where(rsquared.flatten() < rad**2)[0]
        zeros_ind = np.where(data.flatten()[indices] == 0.)[0]
        zeros_ind_len = len(zeros_ind)  
        
        if zeros_ind_len > 100:
            return True #so image not useable as there are holes
        else:
            return False #can use this image
    
    elif 'LASCO_C2' in filename:
        print('LASCO_C2')
        rad1 = 160 #this seems good
        print('rad1:', rad1)
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
    #print('scale_factor:', scale_factor)
    
    if flag == 'subsample':
        reduced_EIT195_data = data[::scale_factor].T[::scale_factor].T #subsampling image; every other row,column
    elif flag == 'interp': #linear interpolation with anti_aliasing and range preserving
        reduced_EIT195_data = rescale(data, (1/scale_factor), order=1, anti_aliasing=True, preserve_range=True)
    elif flag == 'minpool': #min pooling each block
        reduced_EIT195_data = block_reduce(data, block_size=(scale_factor,scale_factor), func=np.min)
    elif flag == 'maxpool': #max pooling each block
        reduced_EIT195_data = block_reduce(data, block_size=(scale_factor,scale_factor), func=np.max)
    return reduced_EIT195_data


def interpn_fcn(n,data1,data2):
    points = (r_[0, 6], arange(n), arange(n)) #6 for six hours, hardcoded #*******************************
    values = stack((data1, data2))
    xi = rollaxis(mgrid[:n, :n], 0, 3).reshape((n**2, 2))
    xi_t2 = c_[2*ones(n**2), xi]
    xi_t4 = c_[4*ones(n**2), xi]
    values_x_t2 = interpolate.interpn(points, values, xi_t2, method='linear')
    values_x_t2 = values_x_t2.reshape((n, n))

    values_x_t4 = interpolate.interpn(points, values, xi_t4, method='linear')
    values_x_t4 = values_x_t4.reshape((n, n))
    return values_x_t2, values_x_t4
 
def seven_stats_calc(data):
    min_data = np.nanmin(data)
    max_data = np.nanmax(data)
    mean_data = np.nanmean(data)
    median_data = np.nanmedian(data)
    std_data = np.nanstd(data)
    if len(np.where(np.array(data).flatten() != np.array(data).flatten())[0]) > 0:
        skew_data = float(skew(data.flatten(),nan_policy='omit').data)
    else:
        skew_data = skew(data.flatten(),nan_policy='omit')
    kurtosis_data = kurtosis(data.flatten(),nan_policy='omit')
    return np.array([min_data, max_data, mean_data, median_data, std_data, skew_data, kurtosis_data])


all_size_sieved_times_EIT195_grand = []
all_2hr_sieved_times_EIT195_grand = []
all_fileids_EIT195_grand = []
result_query_list_EIT195_grand = []

all_size_sieved_times_MDI_grand = []
all_2hr_sieved_times_MDI_grand = []
all_fileids_MDI_grand = []
result_query_list_MDI_grand = []

all_size_sieved_times_LASCO_C2_grand = []
all_2hr_sieved_times_LASCO_C2_grand = []
all_fileids_LASCO_C2_grand = []
result_query_list_LASCO_C2_grand = []

all_size_sieved_times_LASCO_C3_grand = []
all_2hr_sieved_times_LASCO_C3_grand = []
all_fileids_LASCO_C3_grand = []
result_query_list_LASCO_C3_grand = []

all_size_sieved_times_EIT171_grand = []
all_2hr_sieved_times_EIT171_grand = []
all_fileids_EIT171_grand = []
result_query_list_EIT171_grand = []

all_size_sieved_times_EIT304_grand = []
all_2hr_sieved_times_EIT304_grand = []
all_fileids_EIT304_grand = []
result_query_list_EIT304_grand = []

all_size_sieved_times_EIT284_grand = []
all_2hr_sieved_times_EIT284_grand = []
all_fileids_EIT284_grand = []
result_query_list_EIT284_grand = []

error_list_EIT195 = []
error_list_MDI = []
error_list_LASCO_C2 = []
error_list_LASCO_C3 = []
error_list_EIT171 = []
error_list_EIT304 = []
error_list_EIT284 = []

unreadable_file_ids_EIT195 = []
unreadable_file_ids_MDI = []
unreadable_file_ids_LASCO_C2 = []
unreadable_file_ids_LASCO_C3 = []
unreadable_file_ids_EIT171 = []
unreadable_file_ids_EIT304 = []
unreadable_file_ids_EIT284 = []

#these are the files that had holes
holes_EIT195 = []
holes_MDI = []
holes_LASCO_C2 = []
holes_LASCO_C3 = []
holes_EIT171 = []
holes_EIT304 = []
holes_EIT284 = []

#these are the data lists needed for .h5py construction
reduced_EIT195 = []
reduced_MDI = []
reduced_LASCO_C2 = []
reduced_LASCO_C3 = []
reduced_EIT171 = []
reduced_EIT304 = []
reduced_EIT284 = []

######## CREATE LISTS TO STORE THE ACTUAL DATA FOR THE H5PY CUBES!!! ###################
######## print sizes of all of these lists!

home_dir = '/run/media/carl/Seagate/' #user-input
#url_prefix = 'https://seal.nascom.nasa.gov/'

flag = 'subsample' #user input #options: 'subsample' #'interp' #'minpool', #'maxpool'
print('flag:', flag)
target_dimension = 128 #user input #done this way since some imgs are 1024, others are 512
print('target_dimension:', target_dimension)

####### CHANGE INFO -> DATA!!

days_fixed = 31 #change this to --> 60 #31 in Jupyter #60 is hardcoded to be sure that Fido.search returns less than the allowed VSOclient quota of 10k result returns

time_window = 2 #user-defined #two-hour time window as fundamental unit if want to sample at same cadence as 96m MDI images. Otherwise can use natural image cadence
look_ahead_EIT195 = int(np.ceil(time_window*60/12)) #(i.e., if time_window = 2, then look_ahead = 10 since 12 min cadence is minimum for all data products here!! and time_window is two hours so that don't spend time looking through whole list!)
print('look_ahead_EIT195:', look_ahead_EIT195)
look_ahead_MDI = int(np.ceil(time_window*60/96)) #(i.e., if time_window = 2, then look_ahead = 2 since 96 min cadence is minimum for MDI!!
print('look_ahead_MDI:', look_ahead_MDI)
look_ahead_LASCO_C2 = int(np.ceil(time_window*60/20)) #(i.e., if time_window = 2, then look_ahead = 6 since ~20 min cadence is minimum for LASCO C2!!
print('look_ahead_LASCO_C2:', look_ahead_LASCO_C2)
look_ahead_LASCO_C3 = int(np.ceil(time_window*60/20)) #(i.e., if time_window = 2, then look_ahead = 6 since ~20 min cadence is minimum for LASCO C3!!
print('look_ahead_LASCO_C3:', look_ahead_LASCO_C3)
look_ahead_EIT171 = 2 #since EIT171 has 6 hr cadence; 2 is just to be safe.
print('look_ahead_EIT171:', look_ahead_EIT171)
look_ahead_EIT304 = 2 #since EIT304 has 6 hr cadence; 2 is just to be safe. Sometimes cadence as high as 8 min but for the most part 5-6 hrs.
print('look_ahead_EIT304:', look_ahead_EIT304)
look_ahead_EIT284 = 2 #since EIT284 has 6 hr cadence; 2 is just to be safe. 
print('look_ahead_EIT284:', look_ahead_EIT284)


date_time_pre_start = '19960101-0000' #start time input #make this date default ### start from May 1996 or 1997; because EIT195 images in Jan 1996 are wierd!!! --> wierd
date_time_start= parser.parse(date_time_pre_start)
print(date_time_start)

date_time_pre_end = '20110501-0000' #end time input #make this date default
date_time_end = parser.parse(date_time_pre_end)
print(date_time_end)

diff_start_finish_total_sec = (date_time_end - date_time_start).total_seconds()
print(diff_start_finish_total_sec)

two_months_total_sec = timedelta(days = days_fixed).total_seconds()
print(two_months_total_sec)

#np.ceil here to make sure that always will get the end user specified date using a fixed 60 day increment
num_loops = np.ceil(diff_start_finish_total_sec/two_months_total_sec) #num_loops = 94 for 19960101-0000' - '20110501-0000' 
print(num_loops)

time_range = TimeRange(date_time_start, timedelta(days = days_fixed)) #time_range initialized here
print(time_range)

### CREATE THE DIRECTORIES FOR EACH DATA PRODUCT TYPE
for base in ['EIT195', 'MDI_96m','LASCO_C2','LASCO_C3','EIT171','EIT304','EIT284']:
    base_dir = home_dir + base
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

############# PUT IN MAIN FUNCTION HERE!!! #############
for t_value in np.arange(num_loops)[0:1]: ##remove [0:1] #main workhorse loop
    print('t_value:', t_value)
    
    #EIT195
    #'''
    
    EIT195_results = Fido.search(avso.Time(time_range,date_time_start),avso.Source('SOHO'),avso.Instrument('EIT'),avso.Provider('SDAC'),avso.Wavelength(195 * avso.u.Angstrom,195 * avso.u.Angstrom))
    print(EIT195_results)

    EIT195_results_number = EIT195_results.file_num #need to check that file number is non-zero!
    if EIT195_results_number != 0:
        size_list_EIT195 = [elem['size'] for elem in EIT195_results.get_response(0)[:]]
        print(np.unique(size_list_EIT195), len(size_list_EIT195))

        ind_2059_EIT195 = np.where(np.array(size_list_EIT195) == 2059)[0]
        ind_523_EIT195 = np.where(np.array(size_list_EIT195) == 523)[0]
        print(len(ind_2059_EIT195))
        print(len(ind_523_EIT195))
        ind_2059_523_EIT195 = np.sort(list(ind_2059_EIT195) + list(ind_523_EIT195)) #important to sort here!
        print(ind_2059_523_EIT195)

        all_size_sieved_times_EIT195 = [] #local list to populate at each loop
        all_2hr_sieved_times_EIT195 = [] #local list to populate at each loop

        for value in ind_2059_523_EIT195:
            all_size_sieved_times_EIT195.append(EIT195_results.get_response(0)[int(value)]['time']['start'])
        print(all_size_sieved_times_EIT195, len(all_size_sieved_times_EIT195)) 
        all_size_sieved_times_EIT195_copy = all_size_sieved_times_EIT195.copy() #gets reset each two-month period.

        for i,time_value in enumerate(all_size_sieved_times_EIT195_copy):
            local_time_range = TimeRange(str(time_value),timedelta(hours=time_window))

            local_list = []
            for k,time_val in enumerate(all_size_sieved_times_EIT195_copy[i:i+look_ahead_EIT195]):
                if time_val in local_time_range:
                    local_list.append(time_val)
                print(local_list)
            if len(local_list) >= 1:
                for entry in local_list[1:]:
                    all_size_sieved_times_EIT195_copy.remove(entry)
                all_2hr_sieved_times_EIT195.append(local_list[0])
				
        all_2hr_sieved_times_EIT195_times = list(np.unique(all_2hr_sieved_times_EIT195)) #np.unique() does np.array() and np.sort()
        print('all_2hr_sieved_times_EIT195_times:', all_2hr_sieved_times_EIT195_times, len(all_2hr_sieved_times_EIT195_times))

        all_2hr_sieved_times_EIT195_times_inds_list = list(np.hstack([np.where(np.array(all_size_sieved_times_EIT195) == item)[0] for item in all_2hr_sieved_times_EIT195_times]))
        print(all_2hr_sieved_times_EIT195_times_inds_list, len(all_2hr_sieved_times_EIT195_times_inds_list))

        fetch_indices_EIT195 = ind_2059_523_EIT195[all_2hr_sieved_times_EIT195_times_inds_list] #THESE ARE THE INDICES FOR FIDO TO FETCH!
        print(fetch_indices_EIT195, len(fetch_indices_EIT195))

        all_fileids_EIT195 = []
        result_query_list_EIT195 = []

		
        for item in fetch_indices_EIT195: ##### REMOVE [0:10], [0:2], [2:4]!!!!
            query_result = Fido.fetch(EIT195_results[0,int(item)], path=f'{home_dir}EIT195', progress=False) #once fetched its downloaded!
            print('query_result:', query_result, len(query_result)) 
            if len(query_result) > 1: #means that errors have occured #should be something like this!; perhaps just pick up offending index and refetch again!
                print(item) #the offending index
                query_result = Fido.fetch(EIT195_results[0,int(item)], path=f'{home_dir}EIT195', progress=False) #re-fetch
                if len(query_result) > 1:
                    print('Error file not found') #or raise error?
                    error_list_EIT195.append(EIT195_results.get_response(0)[int(item)]['fileid'])

            result_query_list_EIT195.append(query_result) #so all Fido error free results; prior to fits open check and hole check
            print('result_query_list_EIT195:', result_query_list_EIT195, len(result_query_list_EIT195))
            all_fileids_EIT195.append(EIT195_results.get_response(0)[int(item)]['fileid'])
            print('all_fileids_EIT195:', all_fileids_EIT195, len(all_fileids_EIT195))

            axis1_EIT195,axis2_EIT195,data_EIT195 = readfits(query_result[0]) #queries, the way they're constructed here, is that they return one item at a time.
            print('axis1_EIT195,axis2_EIT195,np.shape(data_EIT195):', axis1_EIT195,axis2_EIT195,np.shape(data_EIT195))

            #data_EIT195 = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!

            if data_EIT195 is not None and axis1_EIT195 == axis2_EIT195:

                if not holes(query_result[0]): #so if not True; so no holes; can use image
                    reduced_EIT195_data = data_reducer(data_EIT195,flag,target_dimension,axis1_EIT195)
                    print(reduced_EIT195_data.shape)
                    time_data = EIT195_results.get_response(0)[int(item)]['time']['start']
                    writefits(f'EIT195/SOHO_EIT195_{time_data}_{target_dimension}', reduced_EIT195_data)
                    reduced_EIT195.append(reduced_EIT195_data) #start building up cube for .h5py
                    os.remove(query_result[0]) #delete original downloaded file

                elif holes(query_result[0]): #so if True, if there are holes
                    holes_EIT195.append(EIT195_results.get_response(0)[int(item)]['fileid'])
                    hole_time_val_EIT195 = EIT195_results.get_response(0)[int(item)]['time']['start']
                    print('hole_time_val_EIT195:', hole_time_val_EIT195)
                    all_2hr_sieved_times_EIT195_times.remove(hole_time_val_EIT195)
                    print('len(all_2hr_sieved_times_EIT195_times):', len(all_2hr_sieved_times_EIT195_times))
                    ind_hole_time_val_EIT195 = np.where(np.array(all_size_sieved_times_EIT195) == hole_time_val_EIT195)[0][0]
                    print('ind_hole_time_val_EIT195:', ind_hole_time_val_EIT195)
                    all_2hr_sieved_times_EIT195_times_inds_list.remove(ind_hole_time_val_EIT195)
                    print('len(all_2hr_sieved_times_EIT195_times_inds_list):', len(all_2hr_sieved_times_EIT195_times_inds_list))
                    
                    os.remove(query_result[0]) #delete original downloaded file

                    ind_timespickup_EIT195 = np.where(np.array(all_size_sieved_times_EIT195) == hole_time_val_EIT195)[0][0]
                    print('ind_timespickup_EIT195:', ind_timespickup_EIT195)
                    zoomed_time_range = TimeRange(str(hole_time_val_EIT195),timedelta(hours=time_window))
                    #the zeroth entry didn't have it so that's why plus 1 in the brackets
            
                    fetch_inds_to_try_list = [] #gets reset for each new item
                    for time_val in all_size_sieved_times_EIT195[ind_timespickup_EIT195+1: ind_timespickup_EIT195 + look_ahead_EIT195]:
                        #ind_local_times_to_try_list = []
                        #fetch_inds_to_try_list = [] #gets reset each time_val
                        print('time_val:', time_val)
                        if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                            ind_next_good_time = np.where(np.array(all_size_sieved_times_EIT195) == time_val)[0]
                            print('ind_next_good_time:', ind_next_good_time)
                            #ind_local_times_to_try_list.append(ind_next_good_time)
                            #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                            fetch_indices_EIT195_next_good = ind_2059_523_EIT195[ind_next_good_time]
                            print('fetch_indices_EIT195_next_good:', fetch_indices_EIT195_next_good)

                            fetch_inds_to_try_list.append(fetch_indices_EIT195_next_good)
                            print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))

                    for index in fetch_inds_to_try_list:
                        print('index:', index)
                        query_result = Fido.fetch(EIT195_results[0,int(index)], path=f'{home_dir}EIT195', progress=False)
                        print('query_result:', query_result, len(query_result)) 
                        axis1_EIT195_next_good,axis2_EIT195_next_good,data_EIT195_next_good = readfits(query_result[0])
                        print('data_EIT195_next_good == None:', data_EIT195_next_good.all() == None)

                        #######data_EIT195_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
    
                        if data_EIT195_next_good is not None and axis1_EIT195_next_good == axis2_EIT195_next_good:

                            if not holes(query_result[0]): #so if not True; so no holes; can use image
                                reduced_EIT195_data = data_reducer(data_EIT195_next_good,flag,target_dimension,axis1_EIT195_next_good)
                                print(reduced_EIT195_data.shape)
                                time_data = EIT195_results.get_response(0)[int(index)]['time']['start']
                                writefits(f'EIT195/SOHO_EIT195_{time_data}_{target_dimension}', reduced_EIT195_data)
                                reduced_EIT195.append(reduced_EIT195_data) #start building up cube for .h5py

                                all_2hr_sieved_times_EIT195_times.append(time_val) #unsorted time location
                                print('len(all_2hr_sieved_times_EIT195_times):', len(all_2hr_sieved_times_EIT195_times))
                                all_2hr_sieved_times_EIT195_times_inds_list.append(index)
                                print('len(all_2hr_sieved_times_EIT195_times_inds_list):', len(all_2hr_sieved_times_EIT195_times_inds_list))
                                all_fileids_EIT195.append(EIT195_results.get_response(0)[int(index)]['fileid'])
                                print('len(all_fileids_EIT195):', len(all_fileids_EIT195))
                                result_query_list_EIT195.append(query_result)
                                print('len(result_query_list_EIT195):', len(result_query_list_EIT195))

                                os.remove(query_result[0]) #delete original downloaded file
                                break #don't need to continue, just exiting the for loop then

                            elif holes(query_result[0]): #so if True, if there are holes
                                holes_EIT195.append(EIT195_results.get_response(0)[int(index)]['fileid'])
                                os.remove(query_result[0])
                                continue 


                        elif data_EIT195_next_good is None:
                            unreadable_file_ids_EIT195.append(EIT195_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue #continue the for loop


            elif data_EIT195 is None:
                unreadable_file_ids_EIT195.append(EIT195_results.get_response(0)[int(item)]['fileid'])

                #remove the time and associated index of bad_time_val_EIT195!
                bad_time_val_EIT195 = EIT195_results.get_response(0)[int(item)]['time']['start']
                print('bad_time_val_EIT195:', bad_time_val_EIT195)
                all_2hr_sieved_times_EIT195_times.remove(bad_time_val_EIT195)
                print('len(all_2hr_sieved_times_EIT195_times):', len(all_2hr_sieved_times_EIT195_times))
                ind_bad_time_val_EIT195 = np.where(np.array(all_size_sieved_times_EIT195) == bad_time_val_EIT195)[0][0]
                print('ind_bad_time_val_EIT195:', ind_bad_time_val_EIT195)
                all_2hr_sieved_times_EIT195_times_inds_list.remove(ind_bad_time_val_EIT195)
                print('len(all_2hr_sieved_times_EIT195_times_inds_list):', len(all_2hr_sieved_times_EIT195_times_inds_list))
                
                os.remove(query_result[0]) #delete original downloaded file

                ind_timespickup_EIT195 = np.where(np.array(all_size_sieved_times_EIT195) == bad_time_val_EIT195)[0][0]
                print('ind_timespickup_EIT195:', ind_timespickup_EIT195)
                zoomed_time_range = TimeRange(str(bad_time_val_EIT195),timedelta(hours=time_window))
                #the zeroth entry didn't have it so that's why plus 1 in the brackets

                #ind_local_times_to_try_list = []
                fetch_inds_to_try_list = [] #gets reset for each new item
                for time_val in all_size_sieved_times_EIT195[ind_timespickup_EIT195+1: ind_timespickup_EIT195 + look_ahead_EIT195]:
                    #ind_local_times_to_try_list = []
                    #fetch_inds_to_try_list = [] #gets reset each time_val
                    print('time_val:', time_val)
                    if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                        ind_next_good_time = np.where(np.array(all_size_sieved_times_EIT195) == time_val)[0]
                        print('ind_next_good_time:', ind_next_good_time)
                        #ind_local_times_to_try_list.append(ind_next_good_time)
                        #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                        fetch_indices_EIT195_next_good = ind_2059_523_EIT195[ind_next_good_time]
                        print('fetch_indices_EIT195_next_good:', fetch_indices_EIT195_next_good)

                        fetch_inds_to_try_list.append(fetch_indices_EIT195_next_good)
                        print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))
  
                for index in fetch_inds_to_try_list:
                    print('index:', index)
                    query_result = Fido.fetch(EIT195_results[0,int(index)], path=f'{home_dir}EIT195', progress=False)
                    print('query_result:', query_result, len(query_result)) 
                    axis1_EIT195_next_good,axis2_EIT195_next_good,data_EIT195_next_good = readfits(query_result[0])
                    print('data_EIT195_next_good == None:', data_EIT195_next_good.all() == None)

                    #data_EIT195_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
  
                    if data_EIT195_next_good is not None and axis1_EIT195_next_good == axis2_EIT195_next_good:
 
                        if not holes(query_result[0]): #so if not True; so no holes; can use image
                            reduced_EIT195_data = data_reducer(data_EIT195_next_good,flag,target_dimension,axis1_EIT195_next_good)
                            print(reduced_EIT195_data.shape)
                            time_data = EIT195_results.get_response(0)[int(index)]['time']['start']
                            writefits(f'EIT195/SOHO_EIT195_{time_data}_{target_dimension}', reduced_EIT195_data)
                            reduced_EIT195.append(reduced_EIT195_data) #start building up cube for .h5py

                            all_2hr_sieved_times_EIT195_times.append(time_val) #unsorted time location
                            print('len(all_2hr_sieved_times_EIT195_times):', len(all_2hr_sieved_times_EIT195_times))
                            all_2hr_sieved_times_EIT195_times_inds_list.append(index)
                            print('len(all_2hr_sieved_times_EIT195_times_inds_list):', len(all_2hr_sieved_times_EIT195_times_inds_list))
                            all_fileids_EIT195.append(EIT195_results.get_response(0)[int(index)]['fileid'])
                            print('len(all_fileids_EIT195):', len(all_fileids_EIT195))
                            result_query_list_EIT195.append(query_result)
                            print('len(result_query_list_EIT195):', len(result_query_list_EIT195))

                            os.remove(query_result[0]) #delete original downloaded file
                            break #don't need to continue, just exiting the for loop then

                        elif holes(query_result[0]): #so if True, if there are holes
                            holes_EIT195.append(EIT195_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue 

                    elif data_EIT195_next_good is None:
                        unreadable_file_ids_EIT195.append(EIT195_results.get_response(0)[int(index)]['fileid'])
                        os.remove(query_result[0])
                        continue #continue the for loop
                
    print('all_size_sieved_times_EIT195:', all_size_sieved_times_EIT195, len(all_size_sieved_times_EIT195))
    print('all_2hr_sieved_times_EIT195_times:', all_2hr_sieved_times_EIT195_times, len(all_2hr_sieved_times_EIT195_times))

    print('all_fileids_EIT195:', all_fileids_EIT195, len(all_fileids_EIT195))
    print('error_list_EIT195:', error_list_EIT195)
    print('unreadable_file_ids_EIT195:', unreadable_file_ids_EIT195)

    print('reduced_EIT195:', len(reduced_EIT195)) #reduced_EIT195
    print('holes_EIT195:', holes_EIT195, len(holes_EIT195))
    
    #'''
    #EIT195
    
    
    ##### MDI #####
    #REMOVE THESE COMMENTED LINES
    #date_time_pre_MDI = '20000501-0000' #month in the middle
    #date_time_start = parser.parse(date_time_pre_MDI)
    #time_range = TimeRange(date_time_start, timedelta(days=14))
    #print(time_range)
    
    #MDI
    #'''
    MDI_results = Fido.search(avso.Time(time_range,date_time_start),avso.Source('SOHO'),avso.Instrument('MDI'),avso.Provider('SDAC'),avso.Physobs('LOS_MAGNETIC_FIELD'))
    print(MDI_results)
	
    MDI_results_number = MDI_results.file_num
    if MDI_results_number != 0:
        size_list_MDI = [elem['size'] for elem in MDI_results.get_response(0)[:]]
        print(np.unique(size_list_MDI), len(size_list_MDI))

        ind_4115_MDI = np.where(np.array(size_list_MDI) == 4115.0)[0] #this 4115 value is very consistent for MDI from '96 - '11.
        print(len(ind_4115_MDI))

        all_size_sieved_times_MDI = [] #local list to populate at each loop
        all_2hr_sieved_times_MDI = [] #local list to populate at each loop
        for value in ind_4115_MDI:
            all_size_sieved_times_MDI.append(MDI_results.get_response(0)[int(value)]['time']['start'])
        print(all_size_sieved_times_MDI, len(all_size_sieved_times_MDI)) 
        all_size_sieved_times_MDI_copy = all_size_sieved_times_MDI.copy() #gets reset each loop too!

        for i,time_value in enumerate(all_size_sieved_times_MDI_copy):
            local_time_range = TimeRange(str(time_value),timedelta(hours=time_window))

            local_list = []
            for k,time_val in enumerate(all_size_sieved_times_MDI_copy[i:i+look_ahead_MDI]):
                if time_val in local_time_range:
                    local_list.append(time_val)
                print(local_list)
            if len(local_list) >= 1:
                for entry in local_list[1:]:
                    all_size_sieved_times_MDI_copy.remove(entry)
                    all_2hr_sieved_times_MDI.append(local_list[0])

        all_2hr_sieved_times_MDI_times = np.unique(all_2hr_sieved_times_MDI) #np.unique() does np.array() and np.sort()
        print('all_2hr_sieved_times_MDI_times:', all_2hr_sieved_times_MDI_times, len(all_2hr_sieved_times_MDI_times))

        all_2hr_sieved_times_MDI_times_inds_list = np.hstack([np.where(np.array(all_size_sieved_times_MDI) == item)[0] for item in all_2hr_sieved_times_MDI_times])
        print(all_2hr_sieved_times_MDI_times_inds_list, len(all_2hr_sieved_times_MDI_times_inds_list))

        fetch_indices_MDI = ind_4115_MDI[all_2hr_sieved_times_MDI_times_inds_list] #THESE ARE THE INDICES FOR FIDO TO FETCH!
        print(fetch_indices_MDI, len(fetch_indices_MDI))

        all_fileids_MDI = []
        result_query_list_MDI = []

        for item in fetch_indices_MDI[0:10]: ##### REMOVE [2:4]!!!!
            query_result = Fido.fetch(MDI_results[0,int(item)], path=f'{home_dir}MDI_96m', progress=False)
            print(query_result, len(query_result)) 
            if len(query_result) > 1: #means that errors have occured #should be something like this!; perhaps just pick up offending index and refetch again!
                print(item) #the offending index
                query_result = Fido.fetch(MDI_results[0,int(item)], path=f'{home_dir}MDI_96m', progress=False) #redo the fetch
                if len(query_result) > 1:
                    print('Error file not found') #or raise error?
                    error_list_MDI.append(MDI_results.get_response(0)[int(item)]['fileid'])
            
            result_query_list_MDI.append(query_result) #so all Fido error free results; prior to fits open check and hole check
            print('result_query_list_MDI:', result_query_list_MDI, len(result_query_list_MDI))
            all_fileids_MDI.append(MDI_results.get_response(0)[int(item)]['fileid'])
            print('all_fileids_MDI:', all_fileids_MDI, len(all_fileids_MDI))
            
            axis1_MDI,axis2_MDI,data_MDI = readfits(query_result[0])
            print('axis1_MDI,axis2_MDI,np.shape(data_MDI):', axis1_MDI,axis2_MDI,np.shape(data_MDI))

            #data_MDI = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!

            if data_MDI is not None and axis1_MDI == axis2_MDI:

                if not holes(query_result[0]): #so if not True; so no holes; can use image
                    reduced_MDI_data = data_reducer(data_MDI,flag,target_dimension,axis1_MDI)
                    print(reduced_MDI_data.shape)
                    time_data = MDI_results.get_response(0)[int(item)]['time']['start']
                    writefits(f'MDI_96m/SOHO_MDI_{time_data}_{target_dimension}', reduced_MDI_data)
                    reduced_MDI.append(reduced_MDI_data) #start building up cube for .h5py
                    os.remove(query_result[0]) #delete original downloaded file

                elif holes(query_result[0]): #so if True, if there are holes
                    holes_MDI.append(MDI_results.get_response(0)[int(item)]['fileid'])
                    hole_time_val_MDI = MDI_results.get_response(0)[int(item)]['time']['start']
                    print('hole_time_val_MDI:', hole_time_val_MDI)
                    all_2hr_sieved_times_MDI_times.remove(hole_time_val_MDI)
                    print('len(all_2hr_sieved_times_MDI_times):', len(all_2hr_sieved_times_MDI_times))
                    ind_hole_time_val_MDI = np.where(np.array(all_size_sieved_times_MDI) == hole_time_val_MDI)[0][0]
                    print('ind_hole_time_val_MDI:', ind_hole_time_val_MDI)
                    all_2hr_sieved_times_MDI_times_inds_list.remove(ind_hole_time_val_MDI)
                    print('len(all_2hr_sieved_times_MDI_times_inds_list):', len(all_2hr_sieved_times_MDI_times_inds_list))
                    
                    os.remove(query_result[0]) #delete original downloaded file

                    ind_timespickup_MDI = np.where(np.array(all_size_sieved_times_MDI) == hole_time_val_MDI)[0][0]
                    print('ind_timespickup_MDI:', ind_timespickup_MDI)
                    zoomed_time_range = TimeRange(str(hole_time_val_MDI),timedelta(hours=time_window))
                    #the zeroth entry didn't have it so that's why plus 1 in the brackets

                    fetch_inds_to_try_list = [] #gets reset for each new item
                    for time_val in all_size_sieved_times_MDI[ind_timespickup_MDI+1: ind_timespickup_MDI + look_ahead_MDI]:
                        #ind_local_times_to_try_list = []
                        #fetch_inds_to_try_list = [] #gets reset each time_val
                        print('time_val:', time_val)
                        if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                            ind_next_good_time = np.where(np.array(all_size_sieved_times_MDI) == time_val)[0]
                            print('ind_next_good_time:', ind_next_good_time)
                            #ind_local_times_to_try_list.append(ind_next_good_time)
                            #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                            fetch_indices_MDI_next_good = ind_4115_MDI[ind_next_good_time]
                            print('fetch_indices_MDI_next_good:', fetch_indices_MDI_next_good)

                            fetch_inds_to_try_list.append(fetch_indices_MDI_next_good)
                            print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))

                    for index in fetch_inds_to_try_list:
                        print('index:', index)
                        query_result = Fido.fetch(MDI_results[0,int(index)], path=f'{home_dir}MDI_96m', progress=False)
                        print('query_result:', query_result, len(query_result)) 
                        axis1_MDI_next_good,axis2_MDI_next_good,data_MDI_next_good = readfits(query_result[0])
                        print('data_MDI_next_good == None:', data_MDI_next_good.all() == None)

                        #######data_MDI_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!

                        if data_MDI_next_good is not None and axis1_MDI_next_good == axis2_MDI_next_good:

                            if not holes(query_result[0]): #so if not True; so no holes; can use image
                                reduced_MDI_data = data_reducer(data_MDI_next_good,flag,target_dimension,axis1_MDI_next_good)
                                print(reduced_MDI_data.shape)
                                time_data = MDI_results.get_response(0)[int(index)]['time']['start']
                                writefits(f'MDI_96m/SOHO_MDI_{time_data}_{target_dimension}', reduced_MDI_data)
                                reduced_MDI.append(reduced_MDI_data) #start building up cube for .h5py

                                all_2hr_sieved_times_MDI_times.append(time_val) #unsorted time location
                                print('len(all_2hr_sieved_times_MDI_times):', len(all_2hr_sieved_times_MDI_times))
                                all_2hr_sieved_times_MDI_times_inds_list.append(index)
                                print('len(all_2hr_sieved_times_MDI_times_inds_list):', len(all_2hr_sieved_times_MDI_times_inds_list))
                                all_fileids_MDI.append(MDI_results.get_response(0)[int(index)]['fileid'])
                                print('len(all_fileids_MDI):', len(all_fileids_MDI))
                                result_query_list_MDI.append(query_result)
                                print('len(result_query_list_MDI):', len(result_query_list_MDI))

                                os.remove(query_result[0]) #delete original downloaded file
                                break #don't need to continue, just exiting the for loop then

                            elif holes(query_result[0]): #so if True, if there are holes
                                holes_MDI.append(MDI_results.get_response(0)[int(index)]['fileid'])
                                os.remove(query_result[0])
                                continue 


                        elif data_MDI_next_good is None:
                            unreadable_file_ids_MDI.append(MDI_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue #continue the for loop


            elif data_MDI is None:
                unreadable_file_ids_MDI.append(MDI_results.get_response(0)[int(item)]['fileid'])

                #remove the time and associated index of bad_time_val_MDI!
                bad_time_val_MDI = MDI_results.get_response(0)[int(item)]['time']['start']
                print('bad_time_val_MDI:', bad_time_val_MDI)
                all_2hr_sieved_times_MDI_times.remove(bad_time_val_MDI)
                print('len(all_2hr_sieved_times_MDI_times):', len(all_2hr_sieved_times_MDI_times))
                ind_bad_time_val_MDI = np.where(np.array(all_size_sieved_times_MDI) == bad_time_val_MDI)[0][0]
                print('ind_bad_time_val_MDI:', ind_bad_time_val_MDI)
                all_2hr_sieved_times_MDI_times_inds_list.remove(ind_bad_time_val_MDI)
                print('len(all_2hr_sieved_times_MDI_times_inds_list):', len(all_2hr_sieved_times_MDI_times_inds_list))
                
                os.remove(query_result[0]) #delete original downloaded file

                ind_timespickup_MDI = np.where(np.array(all_size_sieved_times_MDI) == bad_time_val_MDI)[0][0]
                print('ind_timespickup_MDI:', ind_timespickup_MDI)
                zoomed_time_range = TimeRange(str(bad_time_val_MDI),timedelta(hours=time_window))
                #the zeroth entry didn't have it so that's why plus 1 in the brackets

                #ind_local_times_to_try_list = []
                fetch_inds_to_try_list = [] #gets reset for each new item
                for time_val in all_size_sieved_times_MDI[ind_timespickup_MDI+1: ind_timespickup_MDI + look_ahead_MDI]:
                    #ind_local_times_to_try_list = []
                    #fetch_inds_to_try_list = [] #gets reset each time_val
                    print('time_val:', time_val)
                    if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                        ind_next_good_time = np.where(np.array(all_size_sieved_times_MDI) == time_val)[0]
                        print('ind_next_good_time:', ind_next_good_time)
                        #ind_local_times_to_try_list.append(ind_next_good_time)
                        #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                        fetch_indices_MDI_next_good = ind_4115_MDI[ind_next_good_time]
                        print('fetch_indices_MDI_next_good:', fetch_indices_MDI_next_good)

                        fetch_inds_to_try_list.append(fetch_indices_MDI_next_good)
                        print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))

                for index in fetch_inds_to_try_list:
                    print('index:', index)
                    query_result = Fido.fetch(MDI_results[0,int(index)], path=f'{home_dir}MDI_96m', progress=False)
                    print('query_result:', query_result, len(query_result)) 
                    axis1_MDI_next_good,axis2_MDI_next_good,data_MDI_next_good = readfits(query_result[0])
                    print('data_MDI_next_good == None:', data_MDI_next_good.all() == None)

                    #data_MDI_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!

                    if data_MDI_next_good is not None and axis1_MDI_next_good == axis2_MDI_next_good:

                        if not holes(query_result[0]): #so if not True; so no holes; can use image
                            reduced_MDI_data = data_reducer(data_MDI_next_good,flag,target_dimension,axis1_MDI_next_good)
                            print(reduced_MDI_data.shape)
                            time_data = MDI_results.get_response(0)[int(index)]['time']['start']
                            writefits(f'MDI_96m/SOHO_MDI_{time_data}_{target_dimension}', reduced_MDI_data)
                            reduced_MDI.append(reduced_MDI_data) #start building up cube for .h5py

                            all_2hr_sieved_times_MDI_times.append(time_val) #unsorted time location
                            print('len(all_2hr_sieved_times_MDI_times):', len(all_2hr_sieved_times_MDI_times))
                            all_2hr_sieved_times_MDI_times_inds_list.append(index)
                            print('len(all_2hr_sieved_times_MDI_times_inds_list):', len(all_2hr_sieved_times_MDI_times_inds_list))
                            all_fileids_MDI.append(MDI_results.get_response(0)[int(index)]['fileid'])
                            print('len(all_fileids_MDI):', len(all_fileids_MDI))
                            result_query_list_MDI.append(query_result)
                            print('len(result_query_list_MDI):', len(result_query_list_MDI))

                            os.remove(query_result[0]) #delete original downloaded file
                            break #don't need to continue, just exiting the for loop then

                        elif holes(query_result[0]): #so if True, if there are holes
                            holes_MDI.append(MDI_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue 

                    elif data_MDI_next_good is None:
                        unreadable_file_ids_MDI.append(MDI_results.get_response(0)[int(index)]['fileid'])
                        os.remove(query_result[0])
                        continue #continue the for loop

    print('all_size_sieved_times_MDI:', all_size_sieved_times_MDI, len(all_size_sieved_times_MDI))
    print('all_2hr_sieved_times_MDI_times:', all_2hr_sieved_times_MDI_times, len(all_2hr_sieved_times_MDI_times))

    print('all_fileids_MDI:', all_fileids_MDI, len(all_fileids_MDI))
    print('error_list_MDI:', error_list_MDI)
    print('unreadable_file_ids_MDI:', unreadable_file_ids_MDI)

    print('reduced_MDI:', len(reduced_MDI)) #reduced_MDI
    print('holes_MDI:', holes_MDI, len(holes_MDI))
    
    #'''
    #MDI
    
    #LASCO_C2
    #'''
    date_time_pre_LASCO_C2 = '20000501-0000' #month in the middle
    date_time_start = parser.parse(date_time_pre_LASCO_C2)
    time_range = TimeRange(date_time_start, timedelta(days=14))
    print(time_range)
    
    LASCO_C2_results = Fido.search(avso.Time(time_range,date_time_start),avso.Provider('SDAC'),avso.Source('SOHO'),avso.Instrument('LASCO'),avso.Detector('C2'))
    print(LASCO_C2_results)
    
    LASCO_C2_results_number = LASCO_C2_results.file_num
    if LASCO_C2_results_number != 0:
        size_list_LASCO_C2 = [int(np.ceil(elem['size'] / 100.0))*100 for elem in LASCO_C2_results.get_response(0)[:]]
        print(np.unique(size_list_LASCO_C2), len(size_list_LASCO_C2))
    
    
        #rounded for LASCO C2 to 2100 since 2056 in earlier years and 2059 in later years are both proper file sizes. Only such sizes in this range! 4106 discontinued in later years!
        ind_2100_LASCO_C2 = np.where(np.array(size_list_LASCO_C2) == 2100.0)[0] 
        print(len(ind_2100_LASCO_C2))

        all_size_sieved_times_LASCO_C2 = [] #local list to populate at each loop
        all_2hr_sieved_times_LASCO_C2 = [] #local list to populate at each loop
        for value in ind_2100_LASCO_C2:
            all_size_sieved_times_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(value)]['time']['start'])
        print(all_size_sieved_times_LASCO_C2, len(all_size_sieved_times_LASCO_C2)) 
        all_size_sieved_times_LASCO_C2_copy = all_size_sieved_times_LASCO_C2.copy() #gets reset each loop too!

        for i,time_value in enumerate(all_size_sieved_times_LASCO_C2_copy):
            local_time_range = TimeRange(str(time_value),timedelta(hours=time_window))

            local_list = []
            for k,time_val in enumerate(all_size_sieved_times_LASCO_C2_copy[i:i+look_ahead_LASCO_C2]):
                if time_val in local_time_range:
                    local_list.append(time_val)
                print(local_list)
            if len(local_list) >= 1:
                for entry in local_list[1:]:
                    all_size_sieved_times_LASCO_C2_copy.remove(entry)
                all_2hr_sieved_times_LASCO_C2.append(local_list[0])

        all_2hr_sieved_times_LASCO_C2_times = np.unique(all_2hr_sieved_times_LASCO_C2) #np.unique() does np.array() and np.sort()
        print('all_2hr_sieved_times_LASCO_C2_times:', all_2hr_sieved_times_LASCO_C2_times, len(all_2hr_sieved_times_LASCO_C2_times))

        all_2hr_sieved_times_LASCO_C2_times_inds_list = np.hstack([np.where(np.array(all_size_sieved_times_LASCO_C2) == item)[0] for item in all_2hr_sieved_times_LASCO_C2_times])
        print(all_2hr_sieved_times_LASCO_C2_times_inds_list, len(all_2hr_sieved_times_LASCO_C2_times_inds_list))

        fetch_indices_LASCO_C2 = ind_2100_LASCO_C2[all_2hr_sieved_times_LASCO_C2_times_inds_list] #THESE ARE THE INDICES FOR FIDO TO FETCH!
        print(fetch_indices_LASCO_C2, len(fetch_indices_LASCO_C2))

        all_fileids_LASCO_C2 = []
        result_query_list_LASCO_C2 = []

        for item in fetch_indices_LASCO_C2[0:10]: ##### REMOVE [0:10] [2:4]!!!!
            query_result = Fido.fetch(LASCO_C2_results[0,int(item)], path=f'{home_dir}LASCO_C2', progress=False)
            print(query_result, len(query_result)) 
            if len(query_result) > 1: #means that errors have occured #should be something like this!; perhaps just pick up offending index and refetch again!
                print(item) #the offending index
                query_result = Fido.fetch(LASCO_C2_results[0,int(item)], path=f'{home_dir}LASCO_C2', progress=False) #redo the fetch
                if len(query_result) > 1:
                    print('Error file not found') #or raise error?
                    error_list_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(item)]['fileid'])
            
            result_query_list_LASCO_C2.append(query_result) #so all Fido error free results; prior to fits open check and hole check
            print('result_query_list_LASCO_C2:', result_query_list_LASCO_C2, len(result_query_list_LASCO_C2))
            all_fileids_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(item)]['fileid'])
            print('all_fileids_LASCO_C2:', all_fileids_LASCO_C2, len(all_fileids_LASCO_C2))
            
            axis1_LASCO_C2,axis2_LASCO_C2,data_LASCO_C2 = readfits(query_result[0])
            print('axis1_LASCO_C2,axis2_LASCO_C2,np.shape(data_LASCO_C2):', axis1_LASCO_C2,axis2_LASCO_C2,np.shape(data_LASCO_C2))
            
            #data_LASCO_C2 = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
            
            
            if data_LASCO_C2 is not None and axis1_LASCO_C2 == axis2_LASCO_C2:

                if not holes(query_result[0]): #so if not True; so no holes; can use image
                    reduced_LASCO_C2_data = data_reducer(data_LASCO_C2,flag,target_dimension,axis1_LASCO_C2)
                    print(reduced_LASCO_C2_data.shape)
                    time_data = LASCO_C2_results.get_response(0)[int(item)]['time']['start']
                    writefits(f'LASCO_C2/SOHO_LASCO_C2_{time_data}_{target_dimension}', reduced_LASCO_C2_data)
                    reduced_LASCO_C2.append(reduced_LASCO_C2_data) #start building up cube for .h5py
                    os.remove(query_result[0]) #delete original downloaded file

                elif holes(query_result[0]): #so if True, if there are holes
                    holes_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(item)]['fileid'])
                    hole_time_val_LASCO_C2 = LASCO_C2_results.get_response(0)[int(item)]['time']['start']
                    print('hole_time_val_LASCO_C2:', hole_time_val_LASCO_C2)
                    all_2hr_sieved_times_LASCO_C2_times.remove(hole_time_val_LASCO_C2)
                    print('len(all_2hr_sieved_times_LASCO_C2_times):', len(all_2hr_sieved_times_LASCO_C2_times))
                    ind_hole_time_val_LASCO_C2 = np.where(np.array(all_size_sieved_times_LASCO_C2) == hole_time_val_LASCO_C2)[0][0]
                    print('ind_hole_time_val_LASCO_C2:', ind_hole_time_val_LASCO_C2)
                    all_2hr_sieved_times_LASCO_C2_times_inds_list.remove(ind_hole_time_val_LASCO_C2)
                    print('len(all_2hr_sieved_times_LASCO_C2_times_inds_list):', len(all_2hr_sieved_times_LASCO_C2_times_inds_list))
                    
                    os.remove(query_result[0]) #delete original downloaded file

                    ind_timespickup_LASCO_C2 = np.where(np.array(all_size_sieved_times_LASCO_C2) == hole_time_val_LASCO_C2)[0][0]
                    print('ind_timespickup_LASCO_C2:', ind_timespickup_LASCO_C2)
                    zoomed_time_range = TimeRange(str(hole_time_val_LASCO_C2),timedelta(hours=time_window))
                    #the zeroth entry didn't have it so that's why plus 1 in the brackets
            
                    fetch_inds_to_try_list = [] #gets reset for each new item
                    for time_val in all_size_sieved_times_LASCO_C2[ind_timespickup_LASCO_C2+1: ind_timespickup_LASCO_C2 + look_ahead_LASCO_C2]:
                        #ind_local_times_to_try_list = []
                        #fetch_inds_to_try_list = [] #gets reset each time_val
                        print('time_val:', time_val)
                        if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                            print('hello')
                            ind_next_good_time = np.where(np.array(all_size_sieved_times_LASCO_C2) == time_val)[0]
                            print('ind_next_good_time:', ind_next_good_time)
                            #ind_local_times_to_try_list.append(ind_next_good_time)
                            #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                            fetch_indices_LASCO_C2_next_good = ind_2059_523_LASCO_C2[ind_next_good_time]
                            print('fetch_indices_LASCO_C2_next_good:', fetch_indices_LASCO_C2_next_good)

                            fetch_inds_to_try_list.append(fetch_indices_LASCO_C2_next_good)
                            print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))

                    for index in fetch_inds_to_try_list:
                        print('index:', index)
                        query_result = Fido.fetch(LASCO_C2_results[0,int(index)], path=f'{home_dir}LASCO_C2', progress=False)
                        print('query_result:', query_result, len(query_result)) 
                        axis1_LASCO_C2_next_good,axis2_LASCO_C2_next_good,data_LASCO_C2_next_good = readfits(query_result[0])
                        print('data_LASCO_C2_next_good == None:', data_LASCO_C2_next_good.all() == None)

                        #######data_LASCO_C2_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
    
                        if data_LASCO_C2_next_good is not None and axis1_LASCO_C2_next_good == axis2_LASCO_C2_next_good:

                            if not holes(query_result[0]): #so if not True; so no holes; can use image
                                reduced_LASCO_C2_data = data_reducer(data_LASCO_C2_next_good,flag,target_dimension,axis1_LASCO_C2_next_good)
                                print(reduced_LASCO_C2_data.shape)
                                time_data = LASCO_C2_results.get_response(0)[int(index)]['time']['start']
                                writefits(f'LASCO_C2/SOHO_LASCO_C2_{time_data}_{target_dimension}', reduced_LASCO_C2_data)
                                reduced_LASCO_C2.append(reduced_LASCO_C2_data) #start building up cube for .h5py

                                all_2hr_sieved_times_LASCO_C2_times.append(time_val) #unsorted time location
                                print('len(all_2hr_sieved_times_LASCO_C2_times):', len(all_2hr_sieved_times_LASCO_C2_times))
                                all_2hr_sieved_times_LASCO_C2_times_inds_list.append(index)
                                print('len(all_2hr_sieved_times_LASCO_C2_times_inds_list):', len(all_2hr_sieved_times_LASCO_C2_times_inds_list))
                                all_fileids_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(index)]['fileid'])
                                print('len(all_fileids_LASCO_C2):', len(all_fileids_LASCO_C2))
                                result_query_list_LASCO_C2.append(query_result)
                                print('len(result_query_list_LASCO_C2):', len(result_query_list_LASCO_C2))

                                os.remove(query_result[0]) #delete original downloaded file
                                break #don't need to continue, just exiting the for loop then

                            elif holes(query_result[0]): #so if True, if there are holes
                                holes_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(index)]['fileid'])
                                os.remove(query_result[0])
                                continue 

                        elif data_LASCO_C2_next_good is None:
                            unreadable_file_ids_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue #continue the for loop

            elif data_LASCO_C2 is None:
                unreadable_file_ids_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(item)]['fileid'])

                #remove the time and associated index of bad_time_val_LASCO_C2!
                bad_time_val_LASCO_C2 = LASCO_C2_results.get_response(0)[int(item)]['time']['start']
                print('bad_time_val_LASCO_C2:', bad_time_val_LASCO_C2)
                all_2hr_sieved_times_LASCO_C2_times.remove(bad_time_val_LASCO_C2)
                print('len(all_2hr_sieved_times_LASCO_C2_times):', len(all_2hr_sieved_times_LASCO_C2_times))
                ind_bad_time_val_LASCO_C2 = np.where(np.array(all_size_sieved_times_LASCO_C2) == bad_time_val_LASCO_C2)[0][0]
                print('ind_bad_time_val_LASCO_C2:', ind_bad_time_val_LASCO_C2)
                all_2hr_sieved_times_LASCO_C2_times_inds_list.remove(ind_bad_time_val_LASCO_C2)
                print('len(all_2hr_sieved_times_LASCO_C2_times_inds_list):', len(all_2hr_sieved_times_LASCO_C2_times_inds_list))
                
                os.remove(query_result[0]) #delete original downloaded file

                ind_timespickup_LASCO_C2 = np.where(np.array(all_size_sieved_times_LASCO_C2) == bad_time_val_LASCO_C2)[0][0]
                print('ind_timespickup_LASCO_C2:', ind_timespickup_LASCO_C2)
                zoomed_time_range = TimeRange(str(bad_time_val_LASCO_C2),timedelta(hours=time_window))
                #the zeroth entry didn't have it so that's why plus 1 in the brackets

                #ind_local_times_to_try_list = []
                fetch_inds_to_try_list = [] #gets reset for each new item
                for time_val in all_size_sieved_times_LASCO_C2[ind_timespickup_LASCO_C2+1: ind_timespickup_LASCO_C2 + look_ahead_LASCO_C2]:
                    #ind_local_times_to_try_list = []
                    #fetch_inds_to_try_list = [] #gets reset each time_val
                    print('time_val:', time_val)
                    if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                        ind_next_good_time = np.where(np.array(all_size_sieved_times_LASCO_C2) == time_val)[0]
                        print('ind_next_good_time:', ind_next_good_time)
                        #ind_local_times_to_try_list.append(ind_next_good_time)
                        #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                        fetch_indices_LASCO_C2_next_good = ind_2059_523_LASCO_C2[ind_next_good_time]
                        print('fetch_indices_LASCO_C2_next_good:', fetch_indices_LASCO_C2_next_good)

                        fetch_inds_to_try_list.append(fetch_indices_LASCO_C2_next_good)
                        print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))
  
                for index in fetch_inds_to_try_list:
                    print('index:', index)
                    query_result = Fido.fetch(LASCO_C2_results[0,int(index)], path=f'{home_dir}LASCO_C2', progress=False)
                    print('query_result:', query_result, len(query_result)) 
                    axis1_LASCO_C2_next_good,axis2_LASCO_C2_next_good,data_LASCO_C2_next_good = readfits(query_result[0])
                    print('data_LASCO_C2_next_good == None:', data_LASCO_C2_next_good.all() == None)

                    #data_LASCO_C2_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
  
                    if data_LASCO_C2_next_good is not None and axis1_LASCO_C2_next_good == axis2_LASCO_C2_next_good:
 
                        if not holes(query_result[0]): #so if not True; so no holes; can use image
                            reduced_LASCO_C2_data = data_reducer(data_LASCO_C2_next_good,flag,target_dimension,axis1_LASCO_C2_next_good)
                            print(reduced_LASCO_C2_data.shape)
                            time_data = LASCO_C2_results.get_response(0)[int(index)]['time']['start']
                            writefits(f'LASCO_C2/SOHO_LASCO_C2_{time_data}_{target_dimension}', reduced_LASCO_C2_data)
                            reduced_LASCO_C2.append(reduced_LASCO_C2_data) #start building up cube for .h5py

                            all_2hr_sieved_times_LASCO_C2_times.append(time_val) #unsorted time location
                            print('len(all_2hr_sieved_times_LASCO_C2_times):', len(all_2hr_sieved_times_LASCO_C2_times))
                            all_2hr_sieved_times_LASCO_C2_times_inds_list.append(index)
                            print('len(all_2hr_sieved_times_LASCO_C2_times_inds_list):', len(all_2hr_sieved_times_LASCO_C2_times_inds_list))
                            all_fileids_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(index)]['fileid'])
                            print('len(all_fileids_LASCO_C2):', len(all_fileids_LASCO_C2))
                            result_query_list_LASCO_C2.append(query_result)
                            print('len(result_query_list_LASCO_C2):', len(result_query_list_LASCO_C2))

                            os.remove(query_result[0]) #delete original downloaded file
                            break #don't need to continue, just exiting the for loop then

                        elif holes(query_result[0]): #so if True, if there are holes
                            holes_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue 

                    elif data_LASCO_C2_next_good is None:
                        unreadable_file_ids_LASCO_C2.append(LASCO_C2_results.get_response(0)[int(index)]['fileid'])
                        os.remove(query_result[0])
                        continue #continue the for loop
                
    print('all_size_sieved_times_LASCO_C2:', all_size_sieved_times_LASCO_C2, len(all_size_sieved_times_LASCO_C2))
    print('all_2hr_sieved_times_LASCO_C2_times:', all_2hr_sieved_times_LASCO_C2_times, len(all_2hr_sieved_times_LASCO_C2_times))

    print('all_fileids_LASCO_C2:', all_fileids_LASCO_C2, len(all_fileids_LASCO_C2))
    print('error_list_LASCO_C2:', error_list_LASCO_C2)
    print('unreadable_file_ids_LASCO_C2:', unreadable_file_ids_LASCO_C2)

    print('reduced_LASCO_C2:', len(reduced_LASCO_C2)) #reduced_LASCO_C2
    print('holes_LASCO_C2:', holes_LASCO_C2, len(holes_LASCO_C2))
    
    #LASCO_C2
    #''' 
    
    #LASCO_C3
    #'''
    date_time_pre_LASCO_C3 = '20040724-0000' #month in the middle
    date_time_start = parser.parse(date_time_pre_LASCO_C3)
    time_range = TimeRange(date_time_start, timedelta(days=14))
    print(time_range)
    
    LASCO_C3_results = Fido.search(avso.Time(time_range,date_time_start),avso.Provider('SDAC'),avso.Source('SOHO'),avso.Instrument('LASCO'),avso.Detector('C3'))
    print(LASCO_C3_results)
    
    LASCO_C3_results_number = LASCO_C3_results.file_num
    if LASCO_C3_results_number != 0:
        size_list_LASCO_C3 = [int(np.ceil(elem['size'] / 100.0))*100 for elem in LASCO_C3_results.get_response(0)[:]]
        print(np.unique(size_list_LASCO_C3), len(size_list_LASCO_C3))
    
    
        #rounded for LASCO C2 to 2100 since 2056 in earlier years and 2059 in later years are both proper file sizes. Only such sizes in this range! 4106 discontinued in later years!
        ind_2100_LASCO_C3 = np.where(np.array(size_list_LASCO_C3) == 2100.0)[0] 
        print(len(ind_2100_LASCO_C3))

        all_size_sieved_times_LASCO_C3 = [] #local list to populate at each loop
        all_2hr_sieved_times_LASCO_C3 = [] #local list to populate at each loop
        for value in ind_2100_LASCO_C3:
            all_size_sieved_times_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(value)]['time']['start'])
        print(all_size_sieved_times_LASCO_C3, len(all_size_sieved_times_LASCO_C3)) 
        all_size_sieved_times_LASCO_C3_copy = all_size_sieved_times_LASCO_C3.copy() #gets reset each loop too!

        for i,time_value in enumerate(all_size_sieved_times_LASCO_C3_copy):
            local_time_range = TimeRange(str(time_value),timedelta(hours=time_window))

            local_list = []
            for k,time_val in enumerate(all_size_sieved_times_LASCO_C3_copy[i:i+look_ahead_LASCO_C3]):
                if time_val in local_time_range:
                    local_list.append(time_val)
                print(local_list)
            if len(local_list) >= 1:
                for entry in local_list[1:]:
                    all_size_sieved_times_LASCO_C3_copy.remove(entry)
                all_2hr_sieved_times_LASCO_C3.append(local_list[0])

        all_2hr_sieved_times_LASCO_C3_times = np.unique(all_2hr_sieved_times_LASCO_C3) #np.unique() does np.array() and np.sort()
        print('all_2hr_sieved_times_LASCO_C3_times:', all_2hr_sieved_times_LASCO_C3_times, len(all_2hr_sieved_times_LASCO_C3_times))

        all_2hr_sieved_times_LASCO_C3_times_inds_list = np.hstack([np.where(np.array(all_size_sieved_times_LASCO_C3) == item)[0] for item in all_2hr_sieved_times_LASCO_C3_times])
        print(all_2hr_sieved_times_LASCO_C3_times_inds_list, len(all_2hr_sieved_times_LASCO_C3_times_inds_list))

        fetch_indices_LASCO_C3 = ind_2100_LASCO_C3[all_2hr_sieved_times_LASCO_C3_times_inds_list] #THESE ARE THE INDICES FOR FIDO TO FETCH!
        print(fetch_indices_LASCO_C3, len(fetch_indices_LASCO_C3))

        all_fileids_LASCO_C3 = []
        result_query_list_LASCO_C3 = []

        for item in fetch_indices_LASCO_C3[0:4]: ##### REMOVE [0:10] [2:4]!!!!
            query_result = Fido.fetch(LASCO_C3_results[0,int(item)], path=f'{home_dir}LASCO_C3', progress=False)
            print(query_result, len(query_result)) 
            if len(query_result) > 1: #means that errors have occured #should be something like this!; perhaps just pick up offending index and refetch again!
                print(item) #the offending index
                query_result = Fido.fetch(LASCO_C3_results[0,int(item)], path=f'{home_dir}LASCO_C3', progress=False) #redo the fetch
                if len(query_result) > 1:
                    print('Error file not found') #or raise error?
                    error_list_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(item)]['fileid'])
            
            result_query_list_LASCO_C3.append(query_result) #so all Fido error free results; prior to fits open check and hole check
            print('result_query_list_LASCO_C3:', result_query_list_LASCO_C3, len(result_query_list_LASCO_C3))
            all_fileids_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(item)]['fileid'])
            print('all_fileids_LASCO_C3:', all_fileids_LASCO_C3, len(all_fileids_LASCO_C3))
            
            axis1_LASCO_C3,axis2_LASCO_C3,data_LASCO_C3 = readfits(query_result[0])
            print('axis1_LASCO_C3,axis2_LASCO_C3,np.shape(data_LASCO_C3):', axis1_LASCO_C3,axis2_LASCO_C3,np.shape(data_LASCO_C3))
            
            #data_LASCO_C3 = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
                       
            if data_LASCO_C3 is not None and axis1_LASCO_C3 == axis2_LASCO_C3:

                if not holes(query_result[0]): #so if not True; so no holes; can use image
                    reduced_LASCO_C3_data = data_reducer(data_LASCO_C3,flag,target_dimension,axis1_LASCO_C3)
                    print(reduced_LASCO_C3_data.shape)
                    time_data = LASCO_C3_results.get_response(0)[int(item)]['time']['start']
                    writefits(f'LASCO_C3/SOHO_LASCO_C3_{time_data}_{target_dimension}', reduced_LASCO_C3_data)
                    reduced_LASCO_C3.append(reduced_LASCO_C3_data) #start building up cube for .h5py
                    os.remove(query_result[0]) #delete original downloaded file

                elif holes(query_result[0]): #so if True, if there are holes
                    holes_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(item)]['fileid'])
                    hole_time_val_LASCO_C3 = LASCO_C3_results.get_response(0)[int(item)]['time']['start']
                    print('hole_time_val_LASCO_C3:', hole_time_val_LASCO_C3)
                    all_2hr_sieved_times_LASCO_C3_times.remove(hole_time_val_LASCO_C3)
                    print('len(all_2hr_sieved_times_LASCO_C3_times):', len(all_2hr_sieved_times_LASCO_C3_times))
                    ind_hole_time_val_LASCO_C3 = np.where(np.array(all_size_sieved_times_LASCO_C3) == hole_time_val_LASCO_C3)[0][0]
                    print('ind_hole_time_val_LASCO_C3:', ind_hole_time_val_LASCO_C3)
                    all_2hr_sieved_times_LASCO_C3_times_inds_list.remove(ind_hole_time_val_LASCO_C3)
                    print('len(all_2hr_sieved_times_LASCO_C3_times_inds_list):', len(all_2hr_sieved_times_LASCO_C3_times_inds_list))
                    
                    os.remove(query_result[0]) #delete original downloaded file

                    ind_timespickup_LASCO_C3 = np.where(np.array(all_size_sieved_times_LASCO_C3) == hole_time_val_LASCO_C3)[0][0]
                    print('ind_timespickup_LASCO_C3:', ind_timespickup_LASCO_C3)
                    zoomed_time_range = TimeRange(str(hole_time_val_LASCO_C3),timedelta(hours=time_window))
                    #the zeroth entry didn't have it so that's why plus 1 in the brackets
            
                    fetch_inds_to_try_list = [] #gets reset for each new item
                    for time_val in all_size_sieved_times_LASCO_C3[ind_timespickup_LASCO_C3+1: ind_timespickup_LASCO_C3 + look_ahead_LASCO_C3]:
                        #ind_local_times_to_try_list = []
                        #fetch_inds_to_try_list = [] #gets reset each time_val
                        print('time_val:', time_val)
                        if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                            print('hello')
                            ind_next_good_time = np.where(np.array(all_size_sieved_times_LASCO_C3) == time_val)[0]
                            print('ind_next_good_time:', ind_next_good_time)
                            #ind_local_times_to_try_list.append(ind_next_good_time)
                            #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                            fetch_indices_LASCO_C3_next_good = ind_2059_523_LASCO_C3[ind_next_good_time]
                            print('fetch_indices_LASCO_C3_next_good:', fetch_indices_LASCO_C3_next_good)

                            fetch_inds_to_try_list.append(fetch_indices_LASCO_C3_next_good)
                            print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))

                    for index in fetch_inds_to_try_list:
                        print('index:', index)
                        query_result = Fido.fetch(LASCO_C3_results[0,int(index)], path=f'{home_dir}LASCO_C3', progress=False)
                        print('query_result:', query_result, len(query_result)) 
                        axis1_LASCO_C3_next_good,axis2_LASCO_C3_next_good,data_LASCO_C3_next_good = readfits(query_result[0])
                        print('data_LASCO_C3_next_good == None:', data_LASCO_C3_next_good.all() == None)

                        #######data_LASCO_C3_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
    
                        if data_LASCO_C3_next_good is not None and axis1_LASCO_C3_next_good == axis2_LASCO_C3_next_good:

                            if not holes(query_result[0]): #so if not True; so no holes; can use image
                                reduced_LASCO_C3_data = data_reducer(data_LASCO_C3_next_good,flag,target_dimension,axis1_LASCO_C3_next_good)
                                print(reduced_LASCO_C3_data.shape)
                                time_data = LASCO_C3_results.get_response(0)[int(index)]['time']['start']
                                writefits(f'LASCO_C3/SOHO_LASCO_C3_{time_data}_{target_dimension}', reduced_LASCO_C3_data)
                                reduced_LASCO_C3.append(reduced_LASCO_C3_data) #start building up cube for .h5py

                                all_2hr_sieved_times_LASCO_C3_times.append(time_val) #unsorted time location
                                print('len(all_2hr_sieved_times_LASCO_C3_times):', len(all_2hr_sieved_times_LASCO_C3_times))
                                all_2hr_sieved_times_LASCO_C3_times_inds_list.append(index)
                                print('len(all_2hr_sieved_times_LASCO_C3_times_inds_list):', len(all_2hr_sieved_times_LASCO_C3_times_inds_list))
                                all_fileids_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(index)]['fileid'])
                                print('len(all_fileids_LASCO_C3):', len(all_fileids_LASCO_C3))
                                result_query_list_LASCO_C3.append(query_result)
                                print('len(result_query_list_LASCO_C3):', len(result_query_list_LASCO_C3))

                                os.remove(query_result[0]) #delete original downloaded file
                                break #don't need to continue, just exiting the for loop then

                            elif holes(query_result[0]): #so if True, if there are holes
                                holes_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(index)]['fileid'])
                                os.remove(query_result[0])
                                continue 

                        elif data_LASCO_C3_next_good is None:
                            unreadable_file_ids_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue #continue the for loop

            elif data_LASCO_C3 is None:
                unreadable_file_ids_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(item)]['fileid'])

                #remove the time and associated index of bad_time_val_LASCO_C3!
                bad_time_val_LASCO_C3 = LASCO_C3_results.get_response(0)[int(item)]['time']['start']
                print('bad_time_val_LASCO_C3:', bad_time_val_LASCO_C3)
                all_2hr_sieved_times_LASCO_C3_times.remove(bad_time_val_LASCO_C3)
                print('len(all_2hr_sieved_times_LASCO_C3_times):', len(all_2hr_sieved_times_LASCO_C3_times))
                ind_bad_time_val_LASCO_C3 = np.where(np.array(all_size_sieved_times_LASCO_C3) == bad_time_val_LASCO_C3)[0][0]
                print('ind_bad_time_val_LASCO_C3:', ind_bad_time_val_LASCO_C3)
                all_2hr_sieved_times_LASCO_C3_times_inds_list.remove(ind_bad_time_val_LASCO_C3)
                print('len(all_2hr_sieved_times_LASCO_C3_times_inds_list):', len(all_2hr_sieved_times_LASCO_C3_times_inds_list))
                
                os.remove(query_result[0]) #delete original downloaded file

                ind_timespickup_LASCO_C3 = np.where(np.array(all_size_sieved_times_LASCO_C3) == bad_time_val_LASCO_C3)[0][0]
                print('ind_timespickup_LASCO_C3:', ind_timespickup_LASCO_C3)
                zoomed_time_range = TimeRange(str(bad_time_val_LASCO_C3),timedelta(hours=time_window))
                #the zeroth entry didn't have it so that's why plus 1 in the brackets

                #ind_local_times_to_try_list = []
                fetch_inds_to_try_list = [] #gets reset for each new item
                for time_val in all_size_sieved_times_LASCO_C3[ind_timespickup_LASCO_C3+1: ind_timespickup_LASCO_C3 + look_ahead_LASCO_C3]:
                    #ind_local_times_to_try_list = []
                    #fetch_inds_to_try_list = [] #gets reset each time_val
                    print('time_val:', time_val)
                    if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                        ind_next_good_time = np.where(np.array(all_size_sieved_times_LASCO_C3) == time_val)[0]
                        print('ind_next_good_time:', ind_next_good_time)
                        #ind_local_times_to_try_list.append(ind_next_good_time)
                        #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                        fetch_indices_LASCO_C3_next_good = ind_2059_523_LASCO_C3[ind_next_good_time]
                        print('fetch_indices_LASCO_C3_next_good:', fetch_indices_LASCO_C3_next_good)

                        fetch_inds_to_try_list.append(fetch_indices_LASCO_C3_next_good)
                        print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))
  
                for index in fetch_inds_to_try_list:
                    print('index:', index)
                    query_result = Fido.fetch(LASCO_C3_results[0,int(index)], path=f'{home_dir}LASCO_C3', progress=False)
                    print('query_result:', query_result, len(query_result)) 
                    axis1_LASCO_C3_next_good,axis2_LASCO_C3_next_good,data_LASCO_C3_next_good = readfits(query_result[0])
                    print('data_LASCO_C3_next_good == None:', data_LASCO_C3_next_good.all() == None)

                    #data_LASCO_C3_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
  
                    if data_LASCO_C3_next_good is not None and axis1_LASCO_C3_next_good == axis2_LASCO_C3_next_good:
 
                        if not holes(query_result[0]): #so if not True; so no holes; can use image
                            reduced_LASCO_C3_data = data_reducer(data_LASCO_C3_next_good,flag,target_dimension,axis1_LASCO_C3_next_good)
                            print(reduced_LASCO_C3_data.shape)
                            time_data = LASCO_C3_results.get_response(0)[int(index)]['time']['start']
                            writefits(f'LASCO_C3/SOHO_LASCO_C3_{time_data}_{target_dimension}', reduced_LASCO_C3_data)
                            reduced_LASCO_C3.append(reduced_LASCO_C3_data) #start building up cube for .h5py

                            all_2hr_sieved_times_LASCO_C3_times.append(time_val) #unsorted time location
                            print('len(all_2hr_sieved_times_LASCO_C3_times):', len(all_2hr_sieved_times_LASCO_C3_times))
                            all_2hr_sieved_times_LASCO_C3_times_inds_list.append(index)
                            print('len(all_2hr_sieved_times_LASCO_C3_times_inds_list):', len(all_2hr_sieved_times_LASCO_C3_times_inds_list))
                            all_fileids_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(index)]['fileid'])
                            print('len(all_fileids_LASCO_C3):', len(all_fileids_LASCO_C3))
                            result_query_list_LASCO_C3.append(query_result)
                            print('len(result_query_list_LASCO_C3):', len(result_query_list_LASCO_C3))

                            os.remove(query_result[0]) #delete original downloaded file
                            break #don't need to continue, just exiting the for loop then

                        elif holes(query_result[0]): #so if True, if there are holes
                            holes_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue 

                    elif data_LASCO_C3_next_good is None:
                        unreadable_file_ids_LASCO_C3.append(LASCO_C3_results.get_response(0)[int(index)]['fileid'])
                        os.remove(query_result[0])
                        continue #continue the for loop
                
    print('all_size_sieved_times_LASCO_C3:', all_size_sieved_times_LASCO_C3, len(all_size_sieved_times_LASCO_C3))
    print('all_2hr_sieved_times_LASCO_C3_times:', all_2hr_sieved_times_LASCO_C3_times, len(all_2hr_sieved_times_LASCO_C3_times))

    print('all_fileids_LASCO_C3:', all_fileids_LASCO_C3, len(all_fileids_LASCO_C3))
    print('error_list_LASCO_C3:', error_list_LASCO_C3)
    print('unreadable_file_ids_LASCO_C3:', unreadable_file_ids_LASCO_C3)

    print('reduced_LASCO_C3:', len(reduced_LASCO_C3)) #reduced_LASCO_C3
    print('holes_LASCO_C3:', holes_LASCO_C3, len(holes_LASCO_C3))
    
    #LASCO_C3
    #'''   
    
    #EIT171
    #'''
    
    EIT171_results = Fido.search(avso.Time(time_range,date_time_start),avso.Source('SOHO'),avso.Instrument('EIT'),avso.Provider('SDAC'),avso.Wavelength(171 * avso.u.Angstrom,171 * avso.u.Angstrom))
    print(EIT171_results)

    EIT171_results_number = EIT171_results.file_num #need to check that file number is non-zero!
    if EIT171_results_number != 0:
        size_list_EIT171 = [elem['size'] for elem in EIT171_results.get_response(0)[:]]
        print(np.unique(size_list_EIT171), len(size_list_EIT171))

        ind_2059_EIT171 = np.where(np.array(size_list_EIT171) == 2059)[0]
        print(len(ind_2059_EIT171))

        all_size_sieved_times_EIT171 = [] #local list to populate at each loop
        all_2hr_sieved_times_EIT171 = [] #local list to populate at each loop

        for value in ind_2059_EIT171:
            all_size_sieved_times_EIT171.append(EIT171_results.get_response(0)[int(value)]['time']['start'])
        print(all_size_sieved_times_EIT171, len(all_size_sieved_times_EIT171)) 
        all_size_sieved_times_EIT171_copy = all_size_sieved_times_EIT171.copy() #gets reset each two-month period.

        for i,time_value in enumerate(all_size_sieved_times_EIT171_copy):
            local_time_range = TimeRange(str(time_value),timedelta(hours=time_window))

            local_list = []
            for k,time_val in enumerate(all_size_sieved_times_EIT171_copy[i:i+look_ahead_EIT171]):
                if time_val in local_time_range:
                    local_list.append(time_val)
                print(local_list)
            if len(local_list) >= 1:
                for entry in local_list[1:]:
                    all_size_sieved_times_EIT171_copy.remove(entry)
                all_2hr_sieved_times_EIT171.append(local_list[0])
				
        all_2hr_sieved_times_EIT171_times = list(np.unique(all_2hr_sieved_times_EIT171)) #np.unique() does np.array() and np.sort()
        print('all_2hr_sieved_times_EIT171_times:', all_2hr_sieved_times_EIT171_times, len(all_2hr_sieved_times_EIT171_times))

        all_2hr_sieved_times_EIT171_times_inds_list = list(np.hstack([np.where(np.array(all_size_sieved_times_EIT171) == item)[0] for item in all_2hr_sieved_times_EIT171_times]))
        print(all_2hr_sieved_times_EIT171_times_inds_list, len(all_2hr_sieved_times_EIT171_times_inds_list))

        fetch_indices_EIT171 = ind_2059_EIT171[all_2hr_sieved_times_EIT171_times_inds_list] #THESE ARE THE INDICES FOR FIDO TO FETCH!
        print(fetch_indices_EIT171, len(fetch_indices_EIT171))

        all_fileids_EIT171 = []
        result_query_list_EIT171 = []

		
        for item in fetch_indices_EIT171[0:10]: ##### REMOVE [0:10], [0:2], [2:4]!!!!
            query_result = Fido.fetch(EIT171_results[0,int(item)], path=f'{home_dir}EIT171', progress=False) #once fetched its downloaded!
            print('query_result:', query_result, len(query_result)) 
            if len(query_result) > 1: #means that errors have occured #should be something like this!; perhaps just pick up offending index and refetch again!
                print(item) #the offending index
                query_result = Fido.fetch(EIT171_results[0,int(item)], path=f'{home_dir}EIT171', progress=False) #re-fetch
                if len(query_result) > 1:
                    print('Error file not found') #or raise error?
                    error_list_EIT171.append(EIT171_results.get_response(0)[int(item)]['fileid'])

            result_query_list_EIT171.append(query_result) #so all Fido error free results; prior to fits open check and hole check
            print('result_query_list_EIT171:', result_query_list_EIT171, len(result_query_list_EIT171))
            all_fileids_EIT171.append(EIT171_results.get_response(0)[int(item)]['fileid'])
            print('all_fileids_EIT171:', all_fileids_EIT171, len(all_fileids_EIT171))

            axis1_EIT171,axis2_EIT171,data_EIT171 = readfits(query_result[0]) #queries, the way they're constructed here, is that they return one item at a time.
            print('axis1_EIT171,axis2_EIT171,np.shape(data_EIT171):', axis1_EIT171,axis2_EIT171,np.shape(data_EIT171))

            #data_EIT171 = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!

            if data_EIT171 is not None and axis1_EIT171 == axis2_EIT171:

                if not holes(query_result[0]): #so if not True; so no holes; can use image
                    reduced_EIT171_data = data_reducer(data_EIT171,flag,target_dimension,axis1_EIT171)
                    print(reduced_EIT171_data.shape)
                    time_data = EIT171_results.get_response(0)[int(item)]['time']['start']
                    writefits(f'EIT171/SOHO_EIT171_{time_data}_{target_dimension}', reduced_EIT171_data)
                    reduced_EIT171.append(reduced_EIT171_data) #start building up cube for .h5py
                    os.remove(query_result[0]) #delete original downloaded file

                elif holes(query_result[0]): #so if True, if there are holes
                    holes_EIT171.append(EIT171_results.get_response(0)[int(item)]['fileid'])
                    hole_time_val_EIT171 = EIT171_results.get_response(0)[int(item)]['time']['start']
                    print('hole_time_val_EIT171:', hole_time_val_EIT171)
                    all_2hr_sieved_times_EIT171_times.remove(hole_time_val_EIT171)
                    print('len(all_2hr_sieved_times_EIT171_times):', len(all_2hr_sieved_times_EIT171_times))
                    ind_hole_time_val_EIT171 = np.where(np.array(all_size_sieved_times_EIT171) == hole_time_val_EIT171)[0][0]
                    print('ind_hole_time_val_EIT171:', ind_hole_time_val_EIT171)
                    all_2hr_sieved_times_EIT171_times_inds_list.remove(ind_hole_time_val_EIT171)
                    print('len(all_2hr_sieved_times_EIT171_times_inds_list):', len(all_2hr_sieved_times_EIT171_times_inds_list))
                    
                    os.remove(query_result[0]) #delete original downloaded file

                    ind_timespickup_EIT171 = np.where(np.array(all_size_sieved_times_EIT171) == hole_time_val_EIT171)[0][0]
                    print('ind_timespickup_EIT171:', ind_timespickup_EIT171)
                    zoomed_time_range = TimeRange(str(hole_time_val_EIT171),timedelta(hours=time_window))
                    #the zeroth entry didn't have it so that's why plus 1 in the brackets
            
                    fetch_inds_to_try_list = [] #gets reset for each new item
                    for time_val in all_size_sieved_times_EIT171[ind_timespickup_EIT171+1: ind_timespickup_EIT171 + look_ahead_EIT171]:
                        #ind_local_times_to_try_list = []
                        #fetch_inds_to_try_list = [] #gets reset each time_val
                        print('time_val:', time_val)
                        if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                            ind_next_good_time = np.where(np.array(all_size_sieved_times_EIT171) == time_val)[0]
                            print('ind_next_good_time:', ind_next_good_time)
                            #ind_local_times_to_try_list.append(ind_next_good_time)
                            #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                            fetch_indices_EIT171_next_good = ind_2059_EIT171[ind_next_good_time]
                            print('fetch_indices_EIT171_next_good:', fetch_indices_EIT171_next_good)

                            fetch_inds_to_try_list.append(fetch_indices_EIT171_next_good)
                            print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))

                    for index in fetch_inds_to_try_list:
                        print('index:', index)
                        query_result = Fido.fetch(EIT171_results[0,int(index)], path=f'{home_dir}EIT171', progress=False)
                        print('query_result:', query_result, len(query_result)) 
                        axis1_EIT171_next_good,axis2_EIT171_next_good,data_EIT171_next_good = readfits(query_result[0])
                        print('data_EIT171_next_good == None:', data_EIT171_next_good.all() == None)

                        #######data_EIT171_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
    
                        if data_EIT171_next_good is not None and axis1_EIT171_next_good == axis2_EIT171_next_good:

                            if not holes(query_result[0]): #so if not True; so no holes; can use image
                                reduced_EIT171_data = data_reducer(data_EIT171_next_good,flag,target_dimension,axis1_EIT171_next_good)
                                print(reduced_EIT171_data.shape)
                                time_data = EIT171_results.get_response(0)[int(index)]['time']['start']
                                writefits(f'EIT171/SOHO_EIT171_{time_data}_{target_dimension}', reduced_EIT171_data)
                                reduced_EIT171.append(reduced_EIT171_data) #start building up cube for .h5py

                                all_2hr_sieved_times_EIT171_times.append(time_val) #unsorted time location
                                print('len(all_2hr_sieved_times_EIT171_times):', len(all_2hr_sieved_times_EIT171_times))
                                all_2hr_sieved_times_EIT171_times_inds_list.append(index)
                                print('len(all_2hr_sieved_times_EIT171_times_inds_list):', len(all_2hr_sieved_times_EIT171_times_inds_list))
                                all_fileids_EIT171.append(EIT171_results.get_response(0)[int(index)]['fileid'])
                                print('len(all_fileids_EIT171):', len(all_fileids_EIT171))
                                result_query_list_EIT171.append(query_result)
                                print('len(result_query_list_EIT171):', len(result_query_list_EIT171))

                                os.remove(query_result[0]) #delete original downloaded file
                                break #don't need to continue, just exiting the for loop then

                            elif holes(query_result[0]): #so if True, if there are holes
                                holes_EIT171.append(EIT171_results.get_response(0)[int(index)]['fileid'])
                                os.remove(query_result[0])
                                continue 


                        elif data_EIT171_next_good is None:
                            unreadable_file_ids_EIT171.append(EIT171_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue #continue the for loop


            elif data_EIT171 is None:
                unreadable_file_ids_EIT171.append(EIT171_results.get_response(0)[int(item)]['fileid'])

                #remove the time and associated index of bad_time_val_EIT171!
                bad_time_val_EIT171 = EIT171_results.get_response(0)[int(item)]['time']['start']
                print('bad_time_val_EIT171:', bad_time_val_EIT171)
                all_2hr_sieved_times_EIT171_times.remove(bad_time_val_EIT171)
                print('len(all_2hr_sieved_times_EIT171_times):', len(all_2hr_sieved_times_EIT171_times))
                ind_bad_time_val_EIT171 = np.where(np.array(all_size_sieved_times_EIT171) == bad_time_val_EIT171)[0][0]
                print('ind_bad_time_val_EIT171:', ind_bad_time_val_EIT171)
                all_2hr_sieved_times_EIT171_times_inds_list.remove(ind_bad_time_val_EIT171)
                print('len(all_2hr_sieved_times_EIT171_times_inds_list):', len(all_2hr_sieved_times_EIT171_times_inds_list))
                
                os.remove(query_result[0]) #delete original downloaded file

                ind_timespickup_EIT171 = np.where(np.array(all_size_sieved_times_EIT171) == bad_time_val_EIT171)[0][0]
                print('ind_timespickup_EIT171:', ind_timespickup_EIT171)
                zoomed_time_range = TimeRange(str(bad_time_val_EIT171),timedelta(hours=time_window))
                #the zeroth entry didn't have it so that's why plus 1 in the brackets

                #ind_local_times_to_try_list = []
                fetch_inds_to_try_list = [] #gets reset for each new item
                for time_val in all_size_sieved_times_EIT171[ind_timespickup_EIT171+1: ind_timespickup_EIT171 + look_ahead_EIT171]:
                    #ind_local_times_to_try_list = []
                    #fetch_inds_to_try_list = [] #gets reset each time_val
                    print('time_val:', time_val)
                    if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                        ind_next_good_time = np.where(np.array(all_size_sieved_times_EIT171) == time_val)[0]
                        print('ind_next_good_time:', ind_next_good_time)
                        #ind_local_times_to_try_list.append(ind_next_good_time)
                        #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                        fetch_indices_EIT171_next_good = ind_2059_EIT171[ind_next_good_time]
                        print('fetch_indices_EIT171_next_good:', fetch_indices_EIT171_next_good)

                        fetch_inds_to_try_list.append(fetch_indices_EIT171_next_good)
                        print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))
  
                for index in fetch_inds_to_try_list:
                    print('index:', index)
                    query_result = Fido.fetch(EIT171_results[0,int(index)], path=f'{home_dir}EIT171', progress=False)
                    print('query_result:', query_result, len(query_result)) 
                    axis1_EIT171_next_good,axis2_EIT171_next_good,data_EIT171_next_good = readfits(query_result[0])
                    print('data_EIT171_next_good == None:', data_EIT171_next_good.all() == None)

                    #data_EIT171_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
  
                    if data_EIT171_next_good is not None and axis1_EIT171_next_good == axis2_EIT171_next_good:
 
                        if not holes(query_result[0]): #so if not True; so no holes; can use image
                            reduced_EIT171_data = data_reducer(data_EIT171_next_good,flag,target_dimension,axis1_EIT171_next_good)
                            print(reduced_EIT171_data.shape)
                            time_data = EIT171_results.get_response(0)[int(index)]['time']['start']
                            writefits(f'EIT171/SOHO_EIT171_{time_data}_{target_dimension}', reduced_EIT171_data)
                            reduced_EIT171.append(reduced_EIT171_data) #start building up cube for .h5py

                            all_2hr_sieved_times_EIT171_times.append(time_val) #unsorted time location
                            print('len(all_2hr_sieved_times_EIT171_times):', len(all_2hr_sieved_times_EIT171_times))
                            all_2hr_sieved_times_EIT171_times_inds_list.append(index)
                            print('len(all_2hr_sieved_times_EIT171_times_inds_list):', len(all_2hr_sieved_times_EIT171_times_inds_list))
                            all_fileids_EIT171.append(EIT171_results.get_response(0)[int(index)]['fileid'])
                            print('len(all_fileids_EIT171):', len(all_fileids_EIT171))
                            result_query_list_EIT171.append(query_result)
                            print('len(result_query_list_EIT171):', len(result_query_list_EIT171))

                            os.remove(query_result[0]) #delete original downloaded file
                            break #don't need to continue, just exiting the for loop then

                        elif holes(query_result[0]): #so if True, if there are holes
                            holes_EIT171.append(EIT171_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue 

                    elif data_EIT171_next_good is None:
                        unreadable_file_ids_EIT171.append(EIT171_results.get_response(0)[int(index)]['fileid'])
                        os.remove(query_result[0])
                        continue #continue the for loop
                
    print('all_size_sieved_times_EIT171:', all_size_sieved_times_EIT171, len(all_size_sieved_times_EIT171))
    print('all_2hr_sieved_times_EIT171_times:', all_2hr_sieved_times_EIT171_times, len(all_2hr_sieved_times_EIT171_times))

    print('all_fileids_EIT171:', all_fileids_EIT171, len(all_fileids_EIT171))
    print('error_list_EIT171:', error_list_EIT171)
    print('unreadable_file_ids_EIT171:', unreadable_file_ids_EIT171)

    print('reduced_EIT171:', len(reduced_EIT171)) #reduced_EIT171
    print('holes_EIT171:', holes_EIT171, len(holes_EIT171))
        
    #EIT171
    #'''
    
    
    
    #EIT304
    #'''
    EIT304_results = Fido.search(avso.Time(time_range,date_time_start),avso.Source('SOHO'),avso.Instrument('EIT'),avso.Provider('SDAC'),avso.Wavelength(304 * avso.u.Angstrom,304 * avso.u.Angstrom))
    print(EIT304_results)

    EIT304_results_number = EIT304_results.file_num #need to check that file number is non-zero!
    if EIT304_results_number != 0:
        size_list_EIT304 = [elem['size'] for elem in EIT304_results.get_response(0)[:]]
        print(np.unique(size_list_EIT304), len(size_list_EIT304))

        ind_2059_EIT304 = np.where(np.array(size_list_EIT304) == 2059)[0]
        print(len(ind_2059_EIT304))

        all_size_sieved_times_EIT304 = [] #local list to populate at each loop
        all_2hr_sieved_times_EIT304 = [] #local list to populate at each loop

        for value in ind_2059_EIT304:
            all_size_sieved_times_EIT304.append(EIT304_results.get_response(0)[int(value)]['time']['start'])
        print(all_size_sieved_times_EIT304, len(all_size_sieved_times_EIT304)) 
        all_size_sieved_times_EIT304_copy = all_size_sieved_times_EIT304.copy() #gets reset each two-month period.

        for i,time_value in enumerate(all_size_sieved_times_EIT304_copy):
            local_time_range = TimeRange(str(time_value),timedelta(hours=time_window))

            local_list = []
            for k,time_val in enumerate(all_size_sieved_times_EIT304_copy[i:i+look_ahead_EIT304]):
                if time_val in local_time_range:
                    local_list.append(time_val)
                print(local_list)
            if len(local_list) >= 1:
                for entry in local_list[1:]:
                    all_size_sieved_times_EIT304_copy.remove(entry)
                all_2hr_sieved_times_EIT304.append(local_list[0])
				
        all_2hr_sieved_times_EIT304_times = list(np.unique(all_2hr_sieved_times_EIT304)) #np.unique() does np.array() and np.sort()
        print('all_2hr_sieved_times_EIT304_times:', all_2hr_sieved_times_EIT304_times, len(all_2hr_sieved_times_EIT304_times))

        all_2hr_sieved_times_EIT304_times_inds_list = list(np.hstack([np.where(np.array(all_size_sieved_times_EIT304) == item)[0] for item in all_2hr_sieved_times_EIT304_times]))
        print(all_2hr_sieved_times_EIT304_times_inds_list, len(all_2hr_sieved_times_EIT304_times_inds_list))

        fetch_indices_EIT304 = ind_2059_EIT304[all_2hr_sieved_times_EIT304_times_inds_list] #THESE ARE THE INDICES FOR FIDO TO FETCH!
        print(fetch_indices_EIT304, len(fetch_indices_EIT304))

        all_fileids_EIT304 = []
        result_query_list_EIT304 = []

		
        for item in fetch_indices_EIT304[0:10]: ##### REMOVE [0:10], [0:2], [2:4]!!!!
            query_result = Fido.fetch(EIT304_results[0,int(item)], path=f'{home_dir}EIT304', progress=False) #once fetched its downloaded!
            print('query_result:', query_result, len(query_result)) 
            if len(query_result) > 1: #means that errors have occured #should be something like this!; perhaps just pick up offending index and refetch again!
                print(item) #the offending index
                query_result = Fido.fetch(EIT304_results[0,int(item)], path=f'{home_dir}EIT304', progress=False) #re-fetch
                if len(query_result) > 1:
                    print('Error file not found') #or raise error?
                    error_list_EIT304.append(EIT304_results.get_response(0)[int(item)]['fileid'])

            result_query_list_EIT304.append(query_result) #so all Fido error free results; prior to fits open check and hole check
            print('result_query_list_EIT304:', result_query_list_EIT304, len(result_query_list_EIT304))
            all_fileids_EIT304.append(EIT304_results.get_response(0)[int(item)]['fileid'])
            print('all_fileids_EIT304:', all_fileids_EIT304, len(all_fileids_EIT304))

            axis1_EIT304,axis2_EIT304,data_EIT304 = readfits(query_result[0]) #queries, the way they're constructed here, is that they return one item at a time.
            print('axis1_EIT304,axis2_EIT304,np.shape(data_EIT304):', axis1_EIT304,axis2_EIT304,np.shape(data_EIT304))

            #data_EIT304 = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!

            if data_EIT304 is not None and axis1_EIT304 == axis2_EIT304:

                if not holes(query_result[0]): #so if not True; so no holes; can use image
                    reduced_EIT304_data = data_reducer(data_EIT304,flag,target_dimension,axis1_EIT304)
                    print(reduced_EIT304_data.shape)
                    time_data = EIT304_results.get_response(0)[int(item)]['time']['start']
                    writefits(f'EIT304/SOHO_EIT304_{time_data}_{target_dimension}', reduced_EIT304_data)
                    reduced_EIT304.append(reduced_EIT304_data) #start building up cube for .h5py
                    os.remove(query_result[0]) #delete original downloaded file

                elif holes(query_result[0]): #so if True, if there are holes
                    holes_EIT304.append(EIT304_results.get_response(0)[int(item)]['fileid'])
                    hole_time_val_EIT304 = EIT304_results.get_response(0)[int(item)]['time']['start']
                    print('hole_time_val_EIT304:', hole_time_val_EIT304)
                    all_2hr_sieved_times_EIT304_times.remove(hole_time_val_EIT304)
                    print('len(all_2hr_sieved_times_EIT304_times):', len(all_2hr_sieved_times_EIT304_times))
                    ind_hole_time_val_EIT304 = np.where(np.array(all_size_sieved_times_EIT304) == hole_time_val_EIT304)[0][0]
                    print('ind_hole_time_val_EIT304:', ind_hole_time_val_EIT304)
                    all_2hr_sieved_times_EIT304_times_inds_list.remove(ind_hole_time_val_EIT304)
                    print('len(all_2hr_sieved_times_EIT304_times_inds_list):', len(all_2hr_sieved_times_EIT304_times_inds_list))
                    
                    os.remove(query_result[0]) #delete original downloaded file

                    ind_timespickup_EIT304 = np.where(np.array(all_size_sieved_times_EIT304) == hole_time_val_EIT304)[0][0]
                    print('ind_timespickup_EIT304:', ind_timespickup_EIT304)
                    zoomed_time_range = TimeRange(str(hole_time_val_EIT304),timedelta(hours=time_window))
                    #the zeroth entry didn't have it so that's why plus 1 in the brackets
            
                    fetch_inds_to_try_list = [] #gets reset for each new item
                    for time_val in all_size_sieved_times_EIT304[ind_timespickup_EIT304+1: ind_timespickup_EIT304 + look_ahead_EIT304]:
                        #ind_local_times_to_try_list = []
                        #fetch_inds_to_try_list = [] #gets reset each time_val
                        print('time_val:', time_val)
                        if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                            ind_next_good_time = np.where(np.array(all_size_sieved_times_EIT304) == time_val)[0]
                            print('ind_next_good_time:', ind_next_good_time)
                            #ind_local_times_to_try_list.append(ind_next_good_time)
                            #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                            fetch_indices_EIT304_next_good = ind_2059_EIT304[ind_next_good_time]
                            print('fetch_indices_EIT304_next_good:', fetch_indices_EIT304_next_good)

                            fetch_inds_to_try_list.append(fetch_indices_EIT304_next_good)
                            print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))

                    for index in fetch_inds_to_try_list:
                        print('index:', index)
                        query_result = Fido.fetch(EIT304_results[0,int(index)], path=f'{home_dir}EIT304', progress=False)
                        print('query_result:', query_result, len(query_result)) 
                        axis1_EIT304_next_good,axis2_EIT304_next_good,data_EIT304_next_good = readfits(query_result[0])
                        print('data_EIT304_next_good == None:', data_EIT304_next_good.all() == None)

                        #######data_EIT304_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
    
                        if data_EIT304_next_good is not None and axis1_EIT304_next_good == axis2_EIT304_next_good:

                            if not holes(query_result[0]): #so if not True; so no holes; can use image
                                reduced_EIT304_data = data_reducer(data_EIT304_next_good,flag,target_dimension,axis1_EIT304_next_good)
                                print(reduced_EIT304_data.shape)
                                time_data = EIT304_results.get_response(0)[int(index)]['time']['start']
                                writefits(f'EIT304/SOHO_EIT304_{time_data}_{target_dimension}', reduced_EIT304_data)
                                reduced_EIT304.append(reduced_EIT304_data) #start building up cube for .h5py

                                all_2hr_sieved_times_EIT304_times.append(time_val) #unsorted time location
                                print('len(all_2hr_sieved_times_EIT304_times):', len(all_2hr_sieved_times_EIT304_times))
                                all_2hr_sieved_times_EIT304_times_inds_list.append(index)
                                print('len(all_2hr_sieved_times_EIT304_times_inds_list):', len(all_2hr_sieved_times_EIT304_times_inds_list))
                                all_fileids_EIT304.append(EIT304_results.get_response(0)[int(index)]['fileid'])
                                print('len(all_fileids_EIT304):', len(all_fileids_EIT304))
                                result_query_list_EIT304.append(query_result)
                                print('len(result_query_list_EIT304):', len(result_query_list_EIT304))

                                os.remove(query_result[0]) #delete original downloaded file
                                break #don't need to continue, just exiting the for loop then

                            elif holes(query_result[0]): #so if True, if there are holes
                                holes_EIT304.append(EIT304_results.get_response(0)[int(index)]['fileid'])
                                os.remove(query_result[0])
                                continue 


                        elif data_EIT304_next_good is None:
                            unreadable_file_ids_EIT304.append(EIT304_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue #continue the for loop


            elif data_EIT304 is None:
                unreadable_file_ids_EIT304.append(EIT304_results.get_response(0)[int(item)]['fileid'])

                #remove the time and associated index of bad_time_val_EIT304!
                bad_time_val_EIT304 = EIT304_results.get_response(0)[int(item)]['time']['start']
                print('bad_time_val_EIT304:', bad_time_val_EIT304)
                all_2hr_sieved_times_EIT304_times.remove(bad_time_val_EIT304)
                print('len(all_2hr_sieved_times_EIT304_times):', len(all_2hr_sieved_times_EIT304_times))
                ind_bad_time_val_EIT304 = np.where(np.array(all_size_sieved_times_EIT304) == bad_time_val_EIT304)[0][0]
                print('ind_bad_time_val_EIT304:', ind_bad_time_val_EIT304)
                all_2hr_sieved_times_EIT304_times_inds_list.remove(ind_bad_time_val_EIT304)
                print('len(all_2hr_sieved_times_EIT304_times_inds_list):', len(all_2hr_sieved_times_EIT304_times_inds_list))
                
                os.remove(query_result[0]) #delete original downloaded file

                ind_timespickup_EIT304 = np.where(np.array(all_size_sieved_times_EIT304) == bad_time_val_EIT304)[0][0]
                print('ind_timespickup_EIT304:', ind_timespickup_EIT304)
                zoomed_time_range = TimeRange(str(bad_time_val_EIT304),timedelta(hours=time_window))
                #the zeroth entry didn't have it so that's why plus 1 in the brackets

                #ind_local_times_to_try_list = []
                fetch_inds_to_try_list = [] #gets reset for each new item
                for time_val in all_size_sieved_times_EIT304[ind_timespickup_EIT304+1: ind_timespickup_EIT304 + look_ahead_EIT304]:
                    #ind_local_times_to_try_list = []
                    #fetch_inds_to_try_list = [] #gets reset each time_val
                    print('time_val:', time_val)
                    if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                        ind_next_good_time = np.where(np.array(all_size_sieved_times_EIT304) == time_val)[0]
                        print('ind_next_good_time:', ind_next_good_time)
                        #ind_local_times_to_try_list.append(ind_next_good_time)
                        #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                        fetch_indices_EIT304_next_good = ind_2059_EIT304[ind_next_good_time]
                        print('fetch_indices_EIT304_next_good:', fetch_indices_EIT304_next_good)

                        fetch_inds_to_try_list.append(fetch_indices_EIT304_next_good)
                        print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))
  
                for index in fetch_inds_to_try_list:
                    print('index:', index)
                    query_result = Fido.fetch(EIT304_results[0,int(index)], path=f'{home_dir}EIT304', progress=False)
                    print('query_result:', query_result, len(query_result)) 
                    axis1_EIT304_next_good,axis2_EIT304_next_good,data_EIT304_next_good = readfits(query_result[0])
                    print('data_EIT304_next_good == None:', data_EIT304_next_good.all() == None)

                    #data_EIT304_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
  
                    if data_EIT304_next_good is not None and axis1_EIT304_next_good == axis2_EIT304_next_good:
 
                        if not holes(query_result[0]): #so if not True; so no holes; can use image
                            reduced_EIT304_data = data_reducer(data_EIT304_next_good,flag,target_dimension,axis1_EIT304_next_good)
                            print(reduced_EIT304_data.shape)
                            time_data = EIT304_results.get_response(0)[int(index)]['time']['start']
                            writefits(f'EIT304/SOHO_EIT304_{time_data}_{target_dimension}', reduced_EIT304_data)
                            reduced_EIT304.append(reduced_EIT304_data) #start building up cube for .h5py

                            all_2hr_sieved_times_EIT304_times.append(time_val) #unsorted time location
                            print('len(all_2hr_sieved_times_EIT304_times):', len(all_2hr_sieved_times_EIT304_times))
                            all_2hr_sieved_times_EIT304_times_inds_list.append(index)
                            print('len(all_2hr_sieved_times_EIT304_times_inds_list):', len(all_2hr_sieved_times_EIT304_times_inds_list))
                            all_fileids_EIT304.append(EIT304_results.get_response(0)[int(index)]['fileid'])
                            print('len(all_fileids_EIT304):', len(all_fileids_EIT304))
                            result_query_list_EIT304.append(query_result)
                            print('len(result_query_list_EIT304):', len(result_query_list_EIT304))

                            os.remove(query_result[0]) #delete original downloaded file
                            break #don't need to continue, just exiting the for loop then

                        elif holes(query_result[0]): #so if True, if there are holes
                            holes_EIT304.append(EIT304_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue 

                    elif data_EIT304_next_good is None:
                        unreadable_file_ids_EIT304.append(EIT304_results.get_response(0)[int(index)]['fileid'])
                        os.remove(query_result[0])
                        continue #continue the for loop
                
    print('all_size_sieved_times_EIT304:', all_size_sieved_times_EIT304, len(all_size_sieved_times_EIT304))
    print('all_2hr_sieved_times_EIT304_times:', all_2hr_sieved_times_EIT304_times, len(all_2hr_sieved_times_EIT304_times))

    print('all_fileids_EIT304:', all_fileids_EIT304, len(all_fileids_EIT304))
    print('error_list_EIT304:', error_list_EIT304)
    print('unreadable_file_ids_EIT304:', unreadable_file_ids_EIT304)

    print('reduced_EIT304:', len(reduced_EIT304)) #reduced_EIT304
    print('holes_EIT304:', holes_EIT304, len(holes_EIT304))    
    
    #EIT304
    #'''
    
    
    
    #EIT284
    #'''
    EIT284_results = Fido.search(avso.Time(time_range,date_time_start),avso.Source('SOHO'),avso.Instrument('EIT'),avso.Provider('SDAC'),avso.Wavelength(284 * avso.u.Angstrom,284 * avso.u.Angstrom))
    print(EIT284_results)

    EIT284_results_number = EIT284_results.file_num #need to check that file number is non-zero!
    if EIT284_results_number != 0:
        size_list_EIT284 = [elem['size'] for elem in EIT284_results.get_response(0)[:]]
        print(np.unique(size_list_EIT284), len(size_list_EIT284))

        ind_2059_EIT284 = np.where(np.array(size_list_EIT284) == 2059)[0]
        print(len(ind_2059_EIT284))

        all_size_sieved_times_EIT284 = [] #local list to populate at each loop
        all_2hr_sieved_times_EIT284 = [] #local list to populate at each loop

        for value in ind_2059_EIT284:
            all_size_sieved_times_EIT284.append(EIT284_results.get_response(0)[int(value)]['time']['start'])
        print(all_size_sieved_times_EIT284, len(all_size_sieved_times_EIT284)) 
        all_size_sieved_times_EIT284_copy = all_size_sieved_times_EIT284.copy() #gets reset each two-month period.

        for i,time_value in enumerate(all_size_sieved_times_EIT284_copy):
            local_time_range = TimeRange(str(time_value),timedelta(hours=time_window))

            local_list = []
            for k,time_val in enumerate(all_size_sieved_times_EIT284_copy[i:i+look_ahead_EIT284]):
                if time_val in local_time_range:
                    local_list.append(time_val)
                print(local_list)
            if len(local_list) >= 1:
                for entry in local_list[1:]:
                    all_size_sieved_times_EIT284_copy.remove(entry)
                all_2hr_sieved_times_EIT284.append(local_list[0])
				
        all_2hr_sieved_times_EIT284_times = list(np.unique(all_2hr_sieved_times_EIT284)) #np.unique() does np.array() and np.sort()
        print('all_2hr_sieved_times_EIT284_times:', all_2hr_sieved_times_EIT284_times, len(all_2hr_sieved_times_EIT284_times))

        all_2hr_sieved_times_EIT284_times_inds_list = list(np.hstack([np.where(np.array(all_size_sieved_times_EIT284) == item)[0] for item in all_2hr_sieved_times_EIT284_times]))
        print(all_2hr_sieved_times_EIT284_times_inds_list, len(all_2hr_sieved_times_EIT284_times_inds_list))

        fetch_indices_EIT284 = ind_2059_EIT284[all_2hr_sieved_times_EIT284_times_inds_list] #THESE ARE THE INDICES FOR FIDO TO FETCH!
        print(fetch_indices_EIT284, len(fetch_indices_EIT284))

        all_fileids_EIT284 = []
        result_query_list_EIT284 = []

		
        for item in fetch_indices_EIT284[0:10]: ##### REMOVE [0:10], [0:2], [2:4]!!!!
            query_result = Fido.fetch(EIT284_results[0,int(item)], path=f'{home_dir}EIT284', progress=False) #once fetched its downloaded!
            print('query_result:', query_result, len(query_result)) 
            if len(query_result) > 1: #means that errors have occured #should be something like this!; perhaps just pick up offending index and refetch again!
                print(item) #the offending index
                query_result = Fido.fetch(EIT284_results[0,int(item)], path=f'{home_dir}EIT284', progress=False) #re-fetch
                if len(query_result) > 1:
                    print('Error file not found') #or raise error?
                    error_list_EIT284.append(EIT284_results.get_response(0)[int(item)]['fileid'])

            result_query_list_EIT284.append(query_result) #so all Fido error free results; prior to fits open check and hole check
            print('result_query_list_EIT284:', result_query_list_EIT284, len(result_query_list_EIT284))
            all_fileids_EIT284.append(EIT284_results.get_response(0)[int(item)]['fileid'])
            print('all_fileids_EIT284:', all_fileids_EIT284, len(all_fileids_EIT284))

            axis1_EIT284,axis2_EIT284,data_EIT284 = readfits(query_result[0]) #queries, the way they're constructed here, is that they return one item at a time.
            print('axis1_EIT284,axis2_EIT284,np.shape(data_EIT284):', axis1_EIT284,axis2_EIT284,np.shape(data_EIT284))

            #data_EIT284 = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!

            if data_EIT284 is not None and axis1_EIT284 == axis2_EIT284:

                if not holes(query_result[0]): #so if not True; so no holes; can use image
                    reduced_EIT284_data = data_reducer(data_EIT284,flag,target_dimension,axis1_EIT284)
                    print(reduced_EIT284_data.shape)
                    time_data = EIT284_results.get_response(0)[int(item)]['time']['start']
                    writefits(f'EIT284/SOHO_EIT284_{time_data}_{target_dimension}', reduced_EIT284_data)
                    reduced_EIT284.append(reduced_EIT284_data) #start building up cube for .h5py
                    os.remove(query_result[0]) #delete original downloaded file

                elif holes(query_result[0]): #so if True, if there are holes
                    holes_EIT284.append(EIT284_results.get_response(0)[int(item)]['fileid'])
                    hole_time_val_EIT284 = EIT284_results.get_response(0)[int(item)]['time']['start']
                    print('hole_time_val_EIT284:', hole_time_val_EIT284)
                    all_2hr_sieved_times_EIT284_times.remove(hole_time_val_EIT284)
                    print('len(all_2hr_sieved_times_EIT284_times):', len(all_2hr_sieved_times_EIT284_times))
                    ind_hole_time_val_EIT284 = np.where(np.array(all_size_sieved_times_EIT284) == hole_time_val_EIT284)[0][0]
                    print('ind_hole_time_val_EIT284:', ind_hole_time_val_EIT284)
                    all_2hr_sieved_times_EIT284_times_inds_list.remove(ind_hole_time_val_EIT284)
                    print('len(all_2hr_sieved_times_EIT284_times_inds_list):', len(all_2hr_sieved_times_EIT284_times_inds_list))
                    
                    os.remove(query_result[0]) #delete original downloaded file

                    ind_timespickup_EIT284 = np.where(np.array(all_size_sieved_times_EIT284) == hole_time_val_EIT284)[0][0]
                    print('ind_timespickup_EIT284:', ind_timespickup_EIT284)
                    zoomed_time_range = TimeRange(str(hole_time_val_EIT284),timedelta(hours=time_window))
                    #the zeroth entry didn't have it so that's why plus 1 in the brackets
            
                    fetch_inds_to_try_list = [] #gets reset for each new item
                    for time_val in all_size_sieved_times_EIT284[ind_timespickup_EIT284+1: ind_timespickup_EIT284 + look_ahead_EIT284]:
                        #ind_local_times_to_try_list = []
                        #fetch_inds_to_try_list = [] #gets reset each time_val
                        print('time_val:', time_val)
                        if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                            ind_next_good_time = np.where(np.array(all_size_sieved_times_EIT284) == time_val)[0]
                            print('ind_next_good_time:', ind_next_good_time)
                            #ind_local_times_to_try_list.append(ind_next_good_time)
                            #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                            fetch_indices_EIT284_next_good = ind_2059_EIT284[ind_next_good_time]
                            print('fetch_indices_EIT284_next_good:', fetch_indices_EIT284_next_good)

                            fetch_inds_to_try_list.append(fetch_indices_EIT284_next_good)
                            print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))

                    for index in fetch_inds_to_try_list:
                        print('index:', index)
                        query_result = Fido.fetch(EIT284_results[0,int(index)], path=f'{home_dir}EIT284', progress=False)
                        print('query_result:', query_result, len(query_result)) 
                        axis1_EIT284_next_good,axis2_EIT284_next_good,data_EIT284_next_good = readfits(query_result[0])
                        print('data_EIT284_next_good == None:', data_EIT284_next_good.all() == None)

                        #######data_EIT284_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
    
                        if data_EIT284_next_good is not None and axis1_EIT284_next_good == axis2_EIT284_next_good:

                            if not holes(query_result[0]): #so if not True; so no holes; can use image
                                reduced_EIT284_data = data_reducer(data_EIT284_next_good,flag,target_dimension,axis1_EIT284_next_good)
                                print(reduced_EIT284_data.shape)
                                time_data = EIT284_results.get_response(0)[int(index)]['time']['start']
                                writefits(f'EIT284/SOHO_EIT284_{time_data}_{target_dimension}', reduced_EIT284_data)
                                reduced_EIT284.append(reduced_EIT284_data) #start building up cube for .h5py

                                all_2hr_sieved_times_EIT284_times.append(time_val) #unsorted time location
                                print('len(all_2hr_sieved_times_EIT284_times):', len(all_2hr_sieved_times_EIT284_times))
                                all_2hr_sieved_times_EIT284_times_inds_list.append(index)
                                print('len(all_2hr_sieved_times_EIT284_times_inds_list):', len(all_2hr_sieved_times_EIT284_times_inds_list))
                                all_fileids_EIT284.append(EIT284_results.get_response(0)[int(index)]['fileid'])
                                print('len(all_fileids_EIT284):', len(all_fileids_EIT284))
                                result_query_list_EIT284.append(query_result)
                                print('len(result_query_list_EIT284):', len(result_query_list_EIT284))

                                os.remove(query_result[0]) #delete original downloaded file
                                break #don't need to continue, just exiting the for loop then

                            elif holes(query_result[0]): #so if True, if there are holes
                                holes_EIT284.append(EIT284_results.get_response(0)[int(index)]['fileid'])
                                os.remove(query_result[0])
                                continue 


                        elif data_EIT284_next_good is None:
                            unreadable_file_ids_EIT284.append(EIT284_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue #continue the for loop


            elif data_EIT284 is None:
                unreadable_file_ids_EIT284.append(EIT284_results.get_response(0)[int(item)]['fileid'])

                #remove the time and associated index of bad_time_val_EIT284!
                bad_time_val_EIT284 = EIT284_results.get_response(0)[int(item)]['time']['start']
                print('bad_time_val_EIT284:', bad_time_val_EIT284)
                all_2hr_sieved_times_EIT284_times.remove(bad_time_val_EIT284)
                print('len(all_2hr_sieved_times_EIT284_times):', len(all_2hr_sieved_times_EIT284_times))
                ind_bad_time_val_EIT284 = np.where(np.array(all_size_sieved_times_EIT284) == bad_time_val_EIT284)[0][0]
                print('ind_bad_time_val_EIT284:', ind_bad_time_val_EIT284)
                all_2hr_sieved_times_EIT284_times_inds_list.remove(ind_bad_time_val_EIT284)
                print('len(all_2hr_sieved_times_EIT284_times_inds_list):', len(all_2hr_sieved_times_EIT284_times_inds_list))
                
                os.remove(query_result[0]) #delete original downloaded file

                ind_timespickup_EIT284 = np.where(np.array(all_size_sieved_times_EIT284) == bad_time_val_EIT284)[0][0]
                print('ind_timespickup_EIT284:', ind_timespickup_EIT284)
                zoomed_time_range = TimeRange(str(bad_time_val_EIT284),timedelta(hours=time_window))
                #the zeroth entry didn't have it so that's why plus 1 in the brackets

                #ind_local_times_to_try_list = []
                fetch_inds_to_try_list = [] #gets reset for each new item
                for time_val in all_size_sieved_times_EIT284[ind_timespickup_EIT284+1: ind_timespickup_EIT284 + look_ahead_EIT284]:
                    #ind_local_times_to_try_list = []
                    #fetch_inds_to_try_list = [] #gets reset each time_val
                    print('time_val:', time_val)
                    if time_val in zoomed_time_range: #this is the next fitting time in the list, slightly less than 2hrs seperated theoretically
                        ind_next_good_time = np.where(np.array(all_size_sieved_times_EIT284) == time_val)[0]
                        print('ind_next_good_time:', ind_next_good_time)
                        #ind_local_times_to_try_list.append(ind_next_good_time)
                        #print('ind_local_times_to_try_list:', ind_local_times_to_try_list)
                        fetch_indices_EIT284_next_good = ind_2059_EIT284[ind_next_good_time]
                        print('fetch_indices_EIT284_next_good:', fetch_indices_EIT284_next_good)

                        fetch_inds_to_try_list.append(fetch_indices_EIT284_next_good)
                        print('fetch_inds_to_try_list:', fetch_inds_to_try_list, len(fetch_inds_to_try_list))
  
                for index in fetch_inds_to_try_list:
                    print('index:', index)
                    query_result = Fido.fetch(EIT284_results[0,int(index)], path=f'{home_dir}EIT284', progress=False)
                    print('query_result:', query_result, len(query_result)) 
                    axis1_EIT284_next_good,axis2_EIT284_next_good,data_EIT284_next_good = readfits(query_result[0])
                    print('data_EIT284_next_good == None:', data_EIT284_next_good.all() == None)

                    #data_EIT284_next_good = None ###### NEED TO CHANGE THIS: COMMENT OUT! Just for testing!
  
                    if data_EIT284_next_good is not None and axis1_EIT284_next_good == axis2_EIT284_next_good:
 
                        if not holes(query_result[0]): #so if not True; so no holes; can use image
                            reduced_EIT284_data = data_reducer(data_EIT284_next_good,flag,target_dimension,axis1_EIT284_next_good)
                            print(reduced_EIT284_data.shape)
                            time_data = EIT284_results.get_response(0)[int(index)]['time']['start']
                            writefits(f'EIT284/SOHO_EIT284_{time_data}_{target_dimension}', reduced_EIT284_data)
                            reduced_EIT284.append(reduced_EIT284_data) #start building up cube for .h5py

                            all_2hr_sieved_times_EIT284_times.append(time_val) #unsorted time location
                            print('len(all_2hr_sieved_times_EIT284_times):', len(all_2hr_sieved_times_EIT284_times))
                            all_2hr_sieved_times_EIT284_times_inds_list.append(index)
                            print('len(all_2hr_sieved_times_EIT284_times_inds_list):', len(all_2hr_sieved_times_EIT284_times_inds_list))
                            all_fileids_EIT284.append(EIT284_results.get_response(0)[int(index)]['fileid'])
                            print('len(all_fileids_EIT284):', len(all_fileids_EIT284))
                            result_query_list_EIT284.append(query_result)
                            print('len(result_query_list_EIT284):', len(result_query_list_EIT284))

                            os.remove(query_result[0]) #delete original downloaded file
                            break #don't need to continue, just exiting the for loop then

                        elif holes(query_result[0]): #so if True, if there are holes
                            holes_EIT284.append(EIT284_results.get_response(0)[int(index)]['fileid'])
                            os.remove(query_result[0])
                            continue 

                    elif data_EIT284_next_good is None:
                        unreadable_file_ids_EIT284.append(EIT284_results.get_response(0)[int(index)]['fileid'])
                        os.remove(query_result[0])
                        continue #continue the for loop
                
    print('all_size_sieved_times_EIT284:', all_size_sieved_times_EIT284, len(all_size_sieved_times_EIT284))
    print('all_2hr_sieved_times_EIT284_times:', all_2hr_sieved_times_EIT284_times, len(all_2hr_sieved_times_EIT284_times))

    print('all_fileids_EIT284:', all_fileids_EIT284, len(all_fileids_EIT284))
    print('error_list_EIT284:', error_list_EIT284)
    print('unreadable_file_ids_EIT284:', unreadable_file_ids_EIT284)

    print('reduced_EIT284:', len(reduced_EIT284)) #reduced_EIT284
    print('holes_EIT284:', holes_EIT284, len(holes_EIT284))
    
    #EIT284
    #''' 

    
    time_range.next() #Sunpy iterator to go for the next 2 months #also have time_range.previous() to go back. #### UNCOMMENT!
    print('time_range.next:', time_range)


			
		

	
	
	
	
	
#This needs to be included, especially the time_range.next() step
'''
all_size_sieved_times_EIT195_grand.append(all_size_sieved_times_EIT195) #list of lists, will need to np.hstack() this grand list!
all_2hr_sieved_times_EIT195_grand.append(all_2hr_sieved_times_EIT195) #list of lists, will need to np.hstack() this grand list!
all_fileids_EIT195_grand.append(all_fileids_EIT195) #list of lists, will need to np.hstack() this grand list!
result_query_list_EIT195_grand.append(result_query_list_EIT195)

all_size_sieved_times_MDI_grand.append(all_size_sieved_times_MDI) #list of lists, will need to np.hstack() this grand list!
all_2hr_sieved_times_MDI_grand.append(all_2hr_sieved_times_MDI) #list of lists, will need to np.hstack() this grand list!
all_fileids_MDI_grand.append(all_fileids_MDI) #list of lists, will need to np.hstack() this grand list!
result_query_list_MDI_grand.append(result_query_list_MDI)


time_range.next() #Sunpy iterator to go for the next 2 months #also have time_range.previous() to go back. #### UNCOMMENT!
print('time_range.next:', time_range) #### UNCOMMENT!
'''

print('len(all_size_sieved_times_EIT195_grand):',len(all_size_sieved_times_EIT195_grand))
print('all_2hr_sieved_times_EIT195_grand:', all_2hr_sieved_times_EIT195_grand, len(all_2hr_sieved_times_EIT195_grand))
print('len(all_fileids_EIT195_grand):', len(all_fileids_EIT195_grand))
print('len(result_query_list_EIT195_grand):', len(result_query_list_EIT195_grand))

print('len(all_size_sieved_times_MDI_grand):',len(all_size_sieved_times_MDI_grand))
print('len(all_2hr_sieved_times_MDI_grand):', len(all_2hr_sieved_times_MDI_grand))
print('len(all_fileids_MDI_grand):', len(all_fileids_MDI_grand))
print('len(result_query_list_MDI_grand):', len(result_query_list_MDI_grand))

print('unreadable_file_ids_EIT195:', unreadable_file_ids_EIT195, len(unreadable_file_ids_EIT195)) 

"""
HERE MAKE THE .H5PY FILES AND COMPRESS THEM!
Export all times that were downloaded -> txt file; sort -n; one column. '\n'
export all fileids that were downloaded -> txt or csv file
export file of all names of reduced products for each product
"""



#### SAVE AS COMPRESSED .H5 FILE
#### PROVIDE OMNI DATASET TRUTH LABEL !!!!
#### DOWNLOAD SCRIPT THAT CALCULATES FRAC DIM AS THIS HAS THE ROUTUINE FOR BLOCK POOLING MIN AND MAX
#### SUBSAMPLE FILE! [::4,8].T[::4,8].T
#### 7 STATS AND FRAC DIM OPTION TO CALCULATE
#### EIT195, MDI, LASCO C2, {then C3 then 284, then 171, 304} -> check good sizes of these!
#### save: 1. original time list, 2. reduced time list (from filtering by size), 3. 'fileid' -> the urls to all files? pre-pended by nascom?
#### COMPRESSED .H5 FILE!

#### time_range.next() forward time range by 2 months each time! to stay within 10k VSOclient query max limit! until reach 1 May 2011 or therabouts as HMI from SDO takes over already in 2010. 

#### hdul_eit = fits.open(downloaded_file_eit195_result) #because path to file is already there!
#### hdr_eit = hdul_eit[0].header

#Check that NAXIS=2 and not 3; otherwise have movie! #size requirement constraint should take care of this issue


#include OMNI2 data set from '96-'11 -> saved on seagate; include on Zenodo!
#follow notes in code to comment out certain lines or uncomment others; like [0:1] or [0:2] or certain times used as place holder currenty!
#also like change days_fixed from 31 -> 60

#package into h5py file and then compress it!
#see if want to have all 7 products inside of it plus accompanying lists of times and OMNI2 data

#remove un-needed print statements 




'''
#look in Andong's script for more switches, ideas of how to switch on options
def main(start_time, end_time, lr, mini_batch, epochs, img_size, home_dir, summary_dir): #img_size = 32x32

	###Clock program
	clock() #initialize clock
	start_pg = clock()

	print(f'start_time: {start_time}, end_time: {end_time}, lr: {lr}, mini_batch: {mini_batch}, epochs: {epochs}, img_size: {img_size}')
	print(f'home_dir: {home_dir}, summary_dir: {summary_dir}')

	experiment_dir = summary_dir+"soho_omni_regular_FULLMDI_minBzminusmean_pooled"+datetime.now().strftime('%Y-%m-%d-%H-%M') #ADDED MDI NAME, NEED TO DO 6 TIMES #fusion
	print(f'experiment_dir: {experiment_dir}')


	start_time = str(start_time)
	end_time = str(end_time)
	lr = float(lr)
	epochs = int(epochs)
	mini_batch = int(mini_batch)
	img_size = int(img_size)

	# Data pre collection
	starttime = ''.join(start_time.split('-'))
	endtime = ''.join(end_time.split('-'))

	home_dir = str(home_dir) #'/home/carl/Documents/' #'/export/scratch2/carl/'



	end_script = clock()
	print('script took {} seconds'.format(end_script))



if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='DL experiment parameters')
	parser.add_argument('--start_date', metavar='time', required=True, help='yyyy-mm-dd')
	parser.add_argument('--end_date', metavar='time', required=True, help='yyyy-mm-dd')
	parser.add_argument('--lr', metavar='learning rate', required=True, help='float val, 0.001')
	parser.add_argument('--miniBatch', metavar='miniBatch', required=True, help='int val, 10')
	parser.add_argument('--epochs', metavar='epochs', required=True, help='int val, 10')
	parser.add_argument('--img_size', metavar='image size', required=False, help='int val, 256', default = 256, type = int)
	parser.add_argument('--dir', metavar='home directory', required=False, help='str, default is /export/scratch2/carl/, need / in the end', default = '/export/scratch2/carl/', type = str )
	parser.add_argument('--sum_dir', metavar='model summaries directory', required=False, help='str, default is /export/scratch1/carl/', default = '/export/scratch1/carl/', type = str )

	args = parser.parse_args()
	main(
		start_time = args.start_date,
		end_time = args.end_date,
		lr = args.lr,
		epochs = args.epochs,
		img_size = args.img_size,
		home_dir = args.dir,
		summary_dir = args.sum_dir,
		mini_batch = args.miniBatch)
'''
