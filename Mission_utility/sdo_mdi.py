import numpy as np
from sunpy.time import TimeRange
from datetime import time, timedelta

"""
Several functions are dependent on the base product and/or mission, and this  class
is meant to organize and simplify those functions. This class is for SDO products,
in addition to MDI. All base-specific properties and functions are defined in this
class file. 
"""
class SDO_MDI:
    
    #static attributes
    class_type = 'SDO_MDI'
    euv_matches = ['94', '131', '171', '193', '211', '304', '335'] 
    uv_matches = ['1600', '1700']
    matches_MDI_HMI_AIA = ['mdi', '96m', 'hmi', '720s', 'aia']
    matches_eit = ['171', '304', '284']
    holes_method_num = 1
    

    
    
    def __init__(self, base_full, lev1_lasco, fits_headers, time_window):
        self.base_full = base_full
        self.lev1_lasco = lev1_lasco
        self.fits_headers = fits_headers
        self.time_window = time_window
        
    
    #dictionary for various scenarios
    
    """
    Uses base/product input information to assign any necessary properties that will
    be used by future functions.
    """
    def set_base_dictionary(self):
        if 'MDI' in self.base_full:
            self.base = 'MDI'
            self.client_export = 'mdi.fd_M_96m_lev182'
            self.mission = 'SOHO'
            self.orig_img_size = 1024
        elif 'HMI' in self.base_full:
            self.base = 'HMI'
            self.client_export = 'hmi.M_720s'
            self.mission = 'SDO'
            self.orig_img_size = 4096
        elif 'AIA' in self.base_full:
            self.base = 'AIA'
            self.wavelen = int(self.base_full[3:])
            self.mission = 'SDO'
            self.orig_img_size = 4096
        else:
            print('Not a valid base name')
            

        
        
    #class methods
    
    """
    Base-specific helper function for holes(), returns True or False depending on if file has
    missing pixel data.
    """
    def are_holes(self, data, x_coord, filename, rsquared, blank_val, missing_vals):
        if any([x in filename for x in self.matches_MDI_HMI_AIA]):
            rad1 = float(x_coord)
            rad2 = 0.6*float(x_coord)
            indices_rad1 = np.where(rsquared.flatten() < rad1**2)[0]
            indices_rad2 = np.where(rsquared.flatten() < rad2**2)[0]
            zeros_ind = np.where(data.flatten()[indices_rad1] == 0.)[0]
            nan_ind = np.where(data.flatten()[indices_rad2] != data.flatten()[indices_rad2])[0]
            blank_ind = np.where(data.flatten()[indices_rad2] == blank_val)[0]
            
            if (len(nan_ind) > 100) or (missing_vals > 0) or (len(blank_ind) > 100): #only nan_ind for calibrated JSOC MDI images which can have many zeros #worked for uncalibrated images: zeros_nan_ind_len > 100:
                return True #so image not useable as there are holes
            else:
                return False #can use this image
    
    """
    Used SunPy's Fido to search for the user specified products
    """
    def product_search(self, time_range,client):
        
        ts = '_'.join(str(time_range.start).split(' '))+'_TAI'
        tf = '_'.join(str(time_range.end).split(' '))+'_TAI'
        
        
        if self.base == 'AIA':
            
            if any([x in self.base_full for x in self.euv_matches]):
                product_results = client.export(f'aia.lev1_euv_12s[{ts}-{tf}@{self.time_window}h][{self.wavelen}]{{image}}')
            
            elif any([x in self.base_full for x in self.uv_matches]):
                product_results = client.export(f'aia.lev1_uv_24s[{ts}-{tf}@{self.time_window}h][{self.wavelen}]{{image}}')   
            
            else:
                product_results = client.export(f'aia.lev1_vis_1h[{ts}-{tf}][{self.wavelen}]{{image}}')
                #product_results = client.export(f'aia.lev1_vis_1h[{ts}-{tf}][{4500}]{{image}}')
        else:
            export_str = self.client_export + f'[{ts}-{tf}]'
            #product_results = client.export(f'mdi.fd_M_96m_lev182[{ts}-{tf}]') 
            product_results = client.export(export_str)                                             
                    

        
        return product_results, client
    
    
    
    
    
    
    """
    Using Fido's returned query object that has now been sieved for proper data size in the case of EIT and LASCO products and user specified time window, the fileid is extracted from the accompanying Fido dictionary and wget is used to retrieve the product. Fido.fetch() method is used to obtain calibrated MDI products.
    """
    def product_retriever(self, product_results,indiv_ind,url_prefix,home_dir,email,size_sieved_df, client): 
        out_dir = f'{home_dir}{self.base_full}_{self.mission}/'
        time_DRMS = size_sieved_df.loc[size_sieved_df['orig_ind']==indiv_ind]['time_at_ind'].values[0]
        
        time_range_for_DRMS = TimeRange(time_DRMS, timedelta(minutes = 1)) #minutes hardcoded to 1 since just fetching one JSOC image at a time and want small time around it.
        ts_DRMS = '_'.join(str(time_range_for_DRMS.start).split(' '))+'_TAI'
        tf_DRMS = '_'.join(str(time_range_for_DRMS.end).split(' '))+'_TAI'
        
        if (self.fits_headers == 'Y') or (self.fits_headers == 'y'):
            
            if self.base == 'AIA':
                
                if any([x in self.base_full for x in self.euv_matches]):
                    client_export_drms = client.export(f'aia.lev1_euv_12s[{ts_DRMS}-{tf_DRMS}][{self.wavelen}]{{image}}', method='url', protocol='fits')
            
                elif any([x in self.base_full for x in self.uv_matches]):
                    client_export_drms = client.export(f'aia.lev1_uv_24s[{ts_DRMS}-{tf_DRMS}][{self.wavelen}]{{image}}', method='url', protocol='fits')
                else:
                    client_export_drms = client.export(f'aia.lev1_vis_1h[{ts_DRMS}-{tf_DRMS}][{self.wavelen}]{{image}}')
            
            else:
                export_str = self.client_export + f'[{ts_DRMS}-{tf_DRMS}]'
                client_export_drms = client.export(export_str, method='url', protocol='fits')
        
        elif (self.fits_headers == 'N') or (self.fits_headers == 'n'):                 
            
            if self.base == 'AIA':
                
                if any([x in self.base_full for x in self.euv_matches]):
                    client_export_drms = client.export(f'aia.lev1_euv_12s[{ts_DRMS}-{tf_DRMS}][{self.wavelen}]{{image}}')
            
                elif any([x in self.base_full for x in self.uv_matches]):
                    client_export_drms = client.export(f'aia.lev1_uv_24s[{ts_DRMS}-{tf_DRMS}][{self.wavelen}]{{image}}')           
            
                else:
                    client_export_drms = client.export(f'aia.lev1_vis_1h[{ts_DRMS}-{tf_DRMS}][{self.wavelen}]{{image}}')
            else:
                export_str = self.client_export + f'[{ts_DRMS}-{tf_DRMS}]'
                client_export_drms = client.export(export_str)

        if client_export_drms.status == 0:
            query_result_pre = client_export_drms.download(out_dir,0) #always the first elemement since just downloading one MDI image at a time
            query_result = list(query_result_pre.download)
        
        elif client_export_drms.status != 0:
            query_result_pre = []
            query_result = []
                            
        
        if (list(size_sieved_df['time_at_ind']) != []) and (query_result == []):
            print('sleep for 15 minutes and then retry DRMS download')
            time.sleep(900) 
            query_result_pre = client_export_drms.download(out_dir,0)
            query_result = list(query_result_pre.download)
    
        return query_result
    
    
    
    
    """
    Returns the indices corresponding to the proper sizes of each of the data products to ensure that they are single, two dimensional images.
    """
    def index_of_sizes(self, product_results):

        calib_file_num = product_results.data.count()['record'] #product_results.file_num works for Fido queries
        
        if calib_file_num == 0:
            ind = []
        else:
            ind = np.arange(calib_file_num)
               
    
        return ind   
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            