from tqdm import tqdm
from time import process_time 
from SOHO_utility.product_time_sync import *

def main(date_start, date_finish, time_step, home_dir, bases, option): #with date_start, date_finish can further regulate which fits files want to start/end with in case different from original prior choice!

    start_process_time = process_time() #initialize clock
    
    product_list = [] #list for all products #will need to hsplit this by base_list_len - won't result in equal division: need another approach
    slice_start_ind_list = [] #to keep track of the corresponding subset of cube slice indices if the original corresponding time range has been modified.
    slice_end_ind_list = []
    
    base_list = bases.split(',')
    print('base_list:', base_list)
    print('time_step in hours:', time_step)
    time_step_sec = time_step*3600
        
    for base in base_list:
        base = base.strip(' ')
        
        time_step_prev = time_step_prev_reader(home_dir, pattern = f'*{base}*[!sync].h5')
        print('time_step_prev:', time_step_prev)
        
        if (option == 'Y') or (option == 'y'): #so deal with .fits files since contain all info that need such as the time and dim. additionally have h5 cube and csv file too.
            data_raw_times = fits_times_reader(home_dir, base) #based on indiv base!
            product_list.append(times_actualizer(data_raw_times, date_start, date_finish)[0]) #so the first return of times_actualizer which now returns two objects
            slice_start_ind_list.append(times_actualizer(data_raw_times, date_start, date_finish)[2][0]) #first start entry
            slice_end_ind_list.append(times_actualizer(data_raw_times, date_start, date_finish)[3][-1]) #last finish entry
                            
            min_time_diff = min_time_step(times_actualizer(data_raw_times, date_start, date_finish)[1])
            print('min_time_diff in hours:', min_time_diff)
            if time_step_prev > time_step: 
                #so times_actualizer()[1] because want the data_times original and not the revised ones in order to have full range when take min.
                raise ValueError("selected time step has to be greater than or equal to the original time step used (e.g., can not choose < 6 hrs if original time step was 6 hrs)")
                
            elif not dimension_checker_from_fits(home_dir, base_list): #strictly the for loop is about a single base but here already raise exception for all bases because need len(base_list) here!
                raise ValueError("there is a mismatch in dimensionality among one or more products (e.g., 128x128 vs 256x256) in the directory")
                
            else:
                continue

        elif (option == 'N') or (option == 'n'): #so dealing with the downloaded data cubes and csv files but no fits files. so need times from csv and dim from h5 cube or csv. both h5 and csv used.
            csv_uniq_times = csv_times_reader(home_dir, pattern = f'*{base}*[!sync].csv')
            product_list.append(times_actualizer(csv_uniq_times, date_start, date_finish)[0])
            slice_start_ind_list.append(times_actualizer(csv_uniq_times, date_start, date_finish)[2][0])
            slice_end_ind_list.append(times_actualizer(csv_uniq_times, date_start, date_finish)[3][-1])
                            
            min_time_diff = min_time_step(times_actualizer(csv_uniq_times, date_start, date_finish)[1])
            print('min_time_diff in hours:', min_time_diff)
            if time_step_prev > time_step: 
                #so times_actualizer()[1] because want the data_times original and not the revised ones in order to have full range when take min.
                raise ValueError("selected time step has to be greater than or equal to the original time step used (e.g., can not choose < 6 hrs if original time step was 6 hrs)")
                                    
            elif not dimension_checker_from_h5cube_csv(home_dir, base_list): #strictly the for loop is about a single base but here already raise exception for all bases
                raise ValueError("there is a mismatch in dimensionality among one or more products reported in the h5 and csv files in the directory (e.g., 128x128 vs 256x256)")                
                
            else:
                continue

    ind_min_len = shortest_prod_list_index_finder(product_list)
        
    synch_time_inds_list, synch_time_list = sync_times_and_inds(product_list, ind_min_len, time_step, time_step_prev)
    synch_time_inds_list_mod, synch_time_list_mod = sync_times_and_inds_sort_by_product(synch_time_inds_list, synch_time_list)

    for i,base in tqdm(enumerate(base_list)):
        base = base.strip(' ')
        cube_data, cube_dim = cube_data_reader(home_dir, base, pattern = f'*{base}*[!sync].h5')
        cube_sync_maker(home_dir, base, cube_data, cube_dim, slice_start_ind_list[i], slice_end_ind_list[i], synch_time_inds_list_mod[i], date_start, date_finish, time_step_prev, time_step)
        csv_time_sync_writer(home_dir, base, date_start, date_finish, cube_dim, synch_time_list_mod[i], time_step_prev, time_step)

    end_process_time = process_time()
    time_of_process = end_process_time - start_process_time
    print(f'time to run in seconds:', time_of_process)
    

if __name__ == '__main__':
    import argparse
    parser_args = argparse.ArgumentParser(description='SOHO Products Synchronization')
    parser_args.add_argument('--date_start', metavar='time', required=True, help='yyyy-mm-dd, 1996-01-01 is earliest start', type = str)
    parser_args.add_argument('--date_finish', metavar='time', required=True, help='yyyy-mm-dd, 2011-05-01 is recommended latest finish', type = str)
    parser_args.add_argument('--time_step', metavar='time', required=True, help='int, time step in hours', type = int) #was previously time_window
    parser_args.add_argument('--home_dir', metavar='home directory', required=True, help='str, e.g., "/home/user/Documents/", need "/" in the end', type = str)
    parser_args.add_argument('--option', metavar='*.fits files present?', required=True, help='str, Y/N or y/n ', type = str)    
    parser_args.add_argument('--products', metavar='products to sync', required=True, help='str, Enter all the following or a subset thereof, in any order, seperated by commas: "EIT195, MDI_96m, LASCO_C2, LASCO_C3, EIT171, EIT304, EIT284"', type = str)   
    
    args = parser_args.parse_args()
    main(
        date_start = args.date_start,
        date_finish = args.date_finish,
        time_step = args.time_step,
        home_dir = args.home_dir,
        option = args.option,
        bases = args.products)
         
