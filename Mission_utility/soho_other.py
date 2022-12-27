from skimage.transform import probabilistic_hough_line
from skimage.feature import canny, blob_log
import numpy as np
import drms
from sunpy.net import Fido, attrs as a #from sunpy.net.vso import attrs as avso
from sunpy.time import TimeRange
import astropy.units as u
from astropy.io import fits
from astropy.io.fits import Header
from datetime import datetime, date, time, timedelta
import shlex, subprocess

class SOHO_no_MDI:
    
    #static attributes
    class_type = 'SOHO_other'
    orig_img_size = 1024
    matches_eit = ['171', '304', '284']
    mission = 'SOHO'
    holes_method_num = 0
    
    
    
    def __init__(self, base_full, lev1_lasco, fits_headers):
        self.base_full = base_full
        self.lev1_lasco = lev1_lasco
        self.fits_headers = fits_headers
        
    #dictionaries for various scenarios
    def set_base_dictionary(self):
        if 'EIT' in self.base_full:
            self.base = 'EIT'
            self.wavelen = int(self.base_full[3:6])
            self.fido_parms = ''
        elif 'LASCO' in self.base_full:
            self.base = 'LASCO'
            self.detector = self.base_full.split('_')[1]

            

        
        
        
    #class methods
    
    def are_holes(self, data, x_coord, filename, rsquared, blank_val, missing_vals):
        if 'LASCO_C3' in filename: 
            rad = 0.8*x_coord
            indices = np.where(rsquared.flatten() < rad**2)[0]
            zeros_ind = np.where(data.flatten()[indices] == 0.)[0]
            zeros_ind_len = len(zeros_ind)  
            nan_ind = np.where(data.flatten()[indices] != data.flatten()[indices])[0]
            nan_ind_len = len(nan_ind)
                    
            if (zeros_ind_len > 100) or (nan_ind_len > 100): ###### For the case if holes are nans in lev1 images. Holes are 0 in lev0 images. ######
                return True #so image not useable as there are holes
            else:
                return False #can use this image   
    
        
        elif 'LASCO_C2' in filename:
            rad1 = 160 #this seems good
            #print('rad1:', rad1)
            rad2 = int(x_coord)
            indices = np.where((rad2**2 > rsquared.flatten()) & (rsquared.flatten() > rad1**2))[0]
            zeros_ind = np.where(data.flatten()[indices] == 0.)[0]
            zeros_ind_len = len(zeros_ind)
            nan_ind = np.where(data.flatten()[indices] != data.flatten()[indices])[0]
            nan_ind_len = len(nan_ind)
         
            if (zeros_ind_len > 100) or (nan_ind_len > 100): ###### For the case if holes are nans in lev1 images. Holes are 0 in lev0 images. ######
                return True #so image not useable as there are holes
            else:
                return False #can use this image
        
   
    
    """
    Used SunPy's Fido to search for the user specified products
    """
    def product_search(self, time_range,client):
        
        
        if self.base == 'EIT':
            wavelen = self.wavelen
            product_results = Fido.search(a.Time(time_range.start,time_range.end), a.Source('SOHO'), a.Instrument('EIT'), a.Provider('SDAC'), a.Wavelength(wavelen * u.Angstrom, wavelen * u.Angstrom))
            #product_results = Fido.search(a.Time(time_range.start,time_range.end), a.Source('SOHO'), a.Instrument('EIT'), a.Provider('SDAC'), a.Wavelength(195 * u.Angstrom, 195 * u.Angstrom))
            
        elif self.base == 'LASCO': 
            detector = self.detector
            print(detector)
            #product_results = Fido.search(a.Time(time_range.start,time_range.end), a.Provider('SDAC'), a.Source('SOHO'), a.Instrument('LASCO'), a.Detector(detector))
            product_results = Fido.search(a.Time(time_range.start,time_range.end), a.Provider('SDAC'), a.Source('SOHO'), a.Instrument('LASCO'), a.Detector('C3'))

        return product_results, client
    
    
    """
    Using Fido's returned query object that has now been sieved for proper data size in the case of EIT and LASCO products and user specified time window, the fileid is extracted from the accompanying Fido dictionary and wget is used to retrieve the product. Fido.fetch() method is used to obtain calibrated MDI products.
    """
    def product_retriever(self, product_results,indiv_ind,url_prefix,home_dir,email,all_size_sieved_times_pre, client): 

        fileid = product_results[0,:][int(indiv_ind)]['fileid']
        item_wget =  url_prefix + fileid
        cmd = 'wget' + ' ' + '--retry-connrefused' + ' ' + '--waitretry=1' + ' ' + '--read-timeout=20' + ' ' + '--timeout=15' + ' ' + '-t' + ' ' + '0' + ' ' + '--continue' + ' ' + item_wget + ' ' + '-P' + ' ' + f'{home_dir}{self.base}_{self.mission}'     
        args = shlex.split(cmd)    
    
        try: 
            wget_output = subprocess.check_output(args, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as err:
            print('Error output:\n', err.output) 
            print('sleep for 15 minutes and then retry command')
            time.sleep(900)
            wget_output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    
        downloaded_fileid = fileid.split('/')[-1]
        query_result = [f'{home_dir}{self.base}_{self.mission}/{downloaded_fileid}']
    
        return query_result
    
    
    
    """
    Returns the indices corresponding to the proper sizes of each of the data products to ensure that they are single, two dimensional images.
    Note: the original version of this code extracted data from VSO using a size unit of kibibyte (it appears that the website may as well). With the update,
    it appears that Fido lists file size in mebibytes instead. Hence, the new version uses mebibytes instead of kibibytes. The Kibibyte values are listed
    for reference, or if looking for files through VSO directly.
    """
    def index_of_sizes(self, product_results): 
        if self.base == 'EIT':
            if '195' in self.base_full:
                size_list = [elem['Size'] for elem in product_results[0,:]]
                size_list_np = [size_list[i].value for i in range(len(size_list))]
                ind_2059 = np.where(np.array(size_list_np) == 2.01074)[0] #this is for 1024 x 1024 pixels #OG: 2059 kibibyte
                ind_523 = np.where(np.array(size_list_np) == 0.51074)[0] ##this is for 512 x 512 pixels #OG: 523 kibibyte
                ind = np.sort(list(ind_2059) + list(ind_523)) #important to sort here since combining two lists!
            #elif (any([x in self.base_full for x in self.matches_EIT])): 
            else:
                size_list = [elem['Size'].value for elem in product_results[0,:]]
                ind = np.where(np.array(size_list) == 2.01074)[0]
        
        elif self.base == 'LASCO': 
            lev1_LASCO = self.lev1_lasco
        
            if (lev1_LASCO == 'Y') or (lev1_LASCO == 'y'):
                size_list = [elem['Size'].value for elem in product_results[0,:]]
                ind = np.where(np.array(size_list) == 4.00977)[0] #4106 is for calibrated level-1.0 #2100.0 was for uncalibrated level-0.5 
            
            elif (lev1_LASCO == 'N') or (lev1_LASCO == 'n'):
                #NUMBERS NOT CORRECT, THIS CANNOT BE USED AS OF 12/12/22
                size_list = [int(np.ceil(elem['Size'].value / 102400.0))/100 for elem in product_results[0,:]]
                ind = np.where(np.array(size_list) == 2100.0)[0] #4106 is for calibrated level-1.0 #2100.0 was for uncalibrated level-0.5         
            
            #[int(np.ceil(elem['size'] / 100.0))*100 for elem in product_results.get_response(0)[:]] #this was for uncalibrated level-0.5
      
    
        return ind
    
   
    
   
    """
    Planet and comet transients filter. Uses Probabilistic Hough Transform on Canny edge processed images to detect straight lines. 
    The lines object returned by probabilistic_hough_line() is a list.
    Lines == False means can use image.
    """
    def planet_comet_transient_filter(self, data_product): #only C2, C3, provide URLS: for C3 lev 1: use theta explicitely for only the horizontal angles since pile on is very well apparent!!! 
    
       
        data_product_log10 = np.log10(data_product)
        edges_data_product_log10 = canny(data_product_log10,1)
        
        if self.detector == 'C3':
            lines = probabilistic_hough_line(edges_data_product_log10, threshold=80, line_length=5, line_gap=0, seed=1, theta=np.array([-np.pi/2, np.pi/2]))
        
                        
        elif self.detector == 'C2':
            lines = probabilistic_hough_line(edges_data_product_log10, threshold=80, line_length=5, line_gap=0, seed=1, theta=np.array([-np.pi/2, np.pi/2]))
        
        if len(lines) > 5:
            return True #can't use this image as straight lines detected
        
        elif not lines:
            return False #can use this image as no straight lines detected
        
        
        
     

    

