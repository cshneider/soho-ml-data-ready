### Retroactively seeds metadata for MDI or HMI. Using the DRMS protocol='fits' is time consuming since a data export request is issued to JSOC 

import drms
import h5py
import numpy as np
import csv
from tqdm import tqdm
import warnings
from pandas.core.common import SettingWithCopyWarning
from Mission_utility.product_time_sync import csv_times_reader
from Mission_utility import downsample_header_local

def main(image_size_output, path_to_mag_cube, mag_cube_name, base, mission):

    warnings.simplefilter(action = "ignore", category = FutureWarning)
    warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

    ### This put all the FITS header keywords, for all the MDI or HMI data products available in a pandas dataframe called mag_keys. ###
    ### This part takes about 10 minutes ### 
    
    client = drms.Client()
    query_mag = 'mdi.fd_M_96m_lev182[]' ### or ('hmi.M_720s') if want SDO HMI instead of SOHO MDI
    mag_keys = client.query(query_mag, key=drms.const.all) 
    print('len(mag_keys):', len(mag_keys))
    
    print('image_size_output:', image_size_output)
    print('path_to_mag_cube:', path_to_mag_cube)
    print('mag_cube_name:', mag_cube_name)
    print('base:', base)
    print('mission:', mission)

    mag_keys_list = list(client.keys('mdi.fd_M_96m_lev182')) #or ('hmi.M_720s') if want SDO HMI instead of SOHO MDI
    print('mag_keys_list:', mag_keys_list)

    cube_orig = h5py.File(f'{path_to_mag_cube}{mag_cube_name}','r') 
    print(list(cube_orig.keys()))

    cube_orig_data = cube_orig[list(cube_orig.keys())[0]][:] #cube_orig[f'{base}_{mission}_{image_size_output}'][:]
    print('np.shape(cube_orig_data):', np.shape(cube_orig_data))

    times_list = csv_times_reader(path_to_mag_cube, pattern = f'*{base}*{mission}*[!sync].csv')

    print('times_list[0:10]:', times_list[0:10])
    print('times_list[-10:]:', times_list[-10:])
    print('np.shape(times_list):', np.shape(times_list))


    ### creat cube copy with data from original cube and add the metadata via attributes which can now write ###
    full_mag_cube_name = f'{path_to_mag_cube}{mag_cube_name}'
    mag_cube_name_new = full_mag_cube_name.split('.')[0] + '_retroactive_metadata.h5'
    print(mag_cube_name_new)

    data_cube_new = h5py.File(mag_cube_name_new,'w')
    data_cube_new.create_dataset(f'{base}_{mission}_{image_size_output}_metadata', data=cube_orig_data, compression="gzip")

    counter = 0
    meta_data_dict = {}
    
    for t_pre in tqdm(times_list): #[0:2] saftey check
        
        t_drms_split = str(drms.to_datetime(t_pre)).split(' ')
        t_tai = '_'.join((t_drms_split[0].replace('-','.'),t_drms_split[1]))+'_TAI'
        
        ### this original method line below never completes on JSOC as starts fast but after 800 files starts to slow down exponentially ###
        #query = client.query(f'mdi.fd_M_96m_lev182[{t_tai}]', key = client.keys('mdi.fd_M_96m_lev182')) with client = drms.Client(email,verbose=False)
        
        query_pre = mag_keys.loc[mag_keys['T_REC'] == t_tai]
        query = mag_keys.loc[query_pre.index[0]]
        
        query_metadata_update = downsample_header_local(mission, image_size_output, query, mag_keys)
        
        for j, key in enumerate(mag_keys):
            if (key == 'COMMENT') or (key == 'HISTORY'):
                key1 = f'{key}{counter}'
                ##########data_cube_new.attrs[f'{key1}_{counter}'] = query_metadata_update[key] #[0]    
                meta_data_dict[f'{key1}_{counter}'] = query_metadata_update[key]               
            else:
                ##########data_cube_new.attrs[f'{key}_{counter}'] = query_metadata_update[key] #[0]
                meta_data_dict[f'{key}_{counter}'] = query_metadata_update[key]
        
        #########data_cube_new.attrs[f'COMMENT_{counter}'] = f'Zeros outside solar disk for {base}'
        meta_data_dict[f'COMMENT_{counter}'] = f'Zeros outside solar disk for {base}'
        
        counter += 1
        
    data_cube_new.attrs.update(meta_data_dict)
    data_cube_new.close()


if __name__ == '__main__':
    import argparse
    parser_args = argparse.ArgumentParser(description='MDI and HMI retroactive metadata seeding')
    parser_args.add_argument('--image_size_output', metavar='image size', required=True, help='e.g., 128, same as entered when made cube', type = int)
    parser_args.add_argument('--cube_dir', metavar='MDI or HMI cube path', required=True, help='str, e.g., "/home/user/Documents/", need "/" in the end', type = str)
    parser_args.add_argument('--cube_name', metavar='name of .h5 cube', required=True, help='str, e.g., 1999-02-02-16:00:00_to_2011-01-01-00:00:00_MDI_96m_subsample_6_LASCOlev1-N_SOHO_128.h5', type = str)
    parser_args.add_argument('--product', metavar='MDI_96m or HMI_720s', required=True, help='str, Enter either MDI_96m or HMI_720s', type = str)
    parser_args.add_argument('--mission', metavar='SOHO or SDO', required=True, help='SOHO or SDO mission', type = str)    
    
    args = parser_args.parse_args()
    main(
        image_size_output = args.image_size_output,
        path_to_mag_cube = args.cube_dir,
        mag_cube_name = args.cube_name,
        base = args.product,
        mission = args.mission) 


