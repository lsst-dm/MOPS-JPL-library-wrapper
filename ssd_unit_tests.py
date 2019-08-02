#!/bin/python

######################################################################## 
# This file contains unit tests for f2py wrapped ssd module 
# kindly made available to LSST by the NASA Jet Propulsion Laboratory, 
# California Institute of Technology   
#
# written by S. Eggl 20190720
#######################################################################

import ssd
import numpy as np
import sys

###################################################################    
# DEFINE UNIT TESTS
################################################################### 

def test_ssd_compute_n_ephemerides(obs_code,Tref,g_param,h_v,epoch_mjd,elements,src,cov_int,n_eph,ut0):
	#precision and accuracy tests for SSD ephemeris generation
	err=0
    	#nominal values for right ascension for ut0 array from JPL Horizons [deg]
	ra_horizons=[260.67263,297.64340,305.5466,297.41694]
	dec_horizons=[-25.86976,-26.52211,-30.45119,-30.54095]
    
  	#nominal values from previous runs [rad]
	ra_nom= [4.54959563, 5.19485734, 5.33279937, 5.19102786]
	dec_nom=[-0.45151257, -0.46289822, -0.53147261, -0.53302494]

    	#call f2py wrapped SSD ephemeris function 
	try:
		[ra, dec, err_ra, err_dec, smaa, smia,pa, v_minus_h, errcod]=ssd.ssd.ssd_compute_n_ephemerides(epoch_mjd,elements, src, obs_code, g_param, ut0, cov_int)
	except:
		print('Error in ssd_compute_n_ephemerides!')
	  
    	#test accuracy        	
	dra=abs(np.array(ra_horizons)-np.rad2deg(np.array(ra)))
	ddec=abs(np.array(dec_horizons)-np.rad2deg(np.array(dec)))

	try:
		if(np.any(dra>1.e-2)):
			err=1
			raise ValueError("Accuracy test failed!")
	except ValueError as ve:
		print(ve)
    
    	#test precision
	dra=abs(np.array(ra_nom)-np.array(ra))
	ddec=abs(np.array(dec_nom)-np.array(dec))

	try:   
		if(np.any(dra>1.e-8)):
    	#precision test failed    
			err=2
			raise ValueError("Precision test failed!")
	except ValueError as ve:
		print(ve)
 
	return err

def test_ssd_propagate_elements(epoch_mjd,elements,ut0):
    	#precision test for SSD orbit element propagation
	err=0
    
    	#cometary orbital elements for four different ut0 = np.array([58468.0,58568.0,58668.0,58768.0])
    	# e,q[au],tp[MJD],node[deg],peri[deg],inc[deg]
	nom_elems=np.array([[7.54302726e-03, 3.22537457e+00, 5.89903144e+04, 5.79572750e+01,2.92000246e+02, 8.94402109e+00],[7.25538845e-03, 3.22458726e+00, 5.90130437e+04, 5.79500333e+01,2.95873173e+02, 8.94218982e+00], [7.02123213e-03, 3.22375590e+00, 5.90336215e+04, 5.79415992e+01,2.99391009e+02, 8.94102085e+00],[6.82773298e-03, 3.22323407e+00, 5.90481155e+04, 5.79337711e+01, 3.01872078e+02, 8.94045418e+00]])
    
	j=0
    	#propagate orbital elements to different epochs:
	for t in ut0:
		try:
			[out_elems,errcod]=ssd.ssd.ssd_propagate_elements(epoch_mjd,elements, t)
        
			if(errcod != 0):
        # orbit element propagation test failed due to internal error   
				err=2
				raise ValueError('Orbit element propagation test failed due to internal error!')
  
			delem=np.array(abs(out_elems-nom_elems[j]))/nom_elems[j]
        
			if(np.any(delem>2.e-9)):
        # orbit element propagation accuracy test failed 
				err=1
				raise ValueError('Orbit element propagation accuracy test failed!')
		except ValueError as ve:
        		print(ve)

		except:
			print('Error in ssd_propagate_elements')
		j=j+1
	return err

###################################################################    
# START UNIT TESTS
###################################################################   

def ssd_unit_test():
###################################################################    
# DEFINE CONSTANTS
################################################################### 
    obs_code='807'
    ut0 = np.array([58468.0,58568.0,58668.0,58768.0])
    cov_int=1
    n_eph=len(ut0)
    Tref=2400000.5

###################################################################    
# DEFINE TARGET ASTEROID
###################################################################    
        
     #TARGET ASTEROID
    tname='Siegfried'

    #spectral slope parameter G
    g_param=0.15

    #absolute magnitude H
    h_v=12.8

    #epoch of orbital elements at Tref [MJD]
    epoch_mjd=58600.0
        
    #cometary orbital elements: e,q[au],tp[MJD],node[deg],peri[deg],inc[deg]       
    elements=[0.007179709583604317,
          3.224300394965053,
          59019.99074555654,
          57.94736169556705,
          297.060245083253,
          8.94174142498718]
    
    #square root covariance matrix for cometary orbital elements
    src=['-7.30758692296977E-9',
         '4.382705220637333E-8',
         '-1.50567315673694E-7',
         '-1.580393656215238E-8',
         '4.990217319865802E-8',
         '-2.144863132279895E-5',
         '-1.522924868782715E-9',
         '4.310160088087949E-9',
         '-.0001765466699452076',
         '-5.220128811753229E-7',
         '-1.227256440856608E-8',
         '3.471999146585647E-8',
         '-.001280090122935523',
         '7.088113260704621E-8',
         '-3.807627393885189E-6',
         '-2.723681462245841E-10',
         '7.394772596575109E-10',
         '1.160018565127471E-5',
         '8.240009079701517E-8',
         '-4.736490292424078E-8',
         '-8.370061313329297E-8']
    

    #initialize ssd
    try:
    	ssd.ssd.ssd_init()
    except:
        print('ssd_init failed.')

    #test ephemeris generation
    try:
        err=test_ssd_compute_n_ephemerides(obs_code,Tref,g_param,h_v,epoch_mjd,elements,src,cov_int,n_eph,ut0)
        if(err!=0):   	
            raise ValueError("ssd_compute_n_ephemerides failed with error code %d" % err) 	 
    except ValueError as ve:
        print(ve)
        sys.exit()
    except:
        print("ssd_compute_n_ephemerides failed with error code %d" % err)
        sys.exit()
            
    #test orbit element propagation
    try:
        err=test_ssd_propagate_elements(epoch_mjd,elements,ut0)
        if(err!=0):   	
            raise ValueError("ssd_propagate_elements failed with error code %d" % err) 	 
    except ValueError as ve:
        print(ve)
        sys.exit()
    except:
        print("ssd_compute_n_ephemerides failed with error code %d" % err)
        sys.exit()
            
    print('All unit tests passed.')
            
    return err

###################################################################    
# PERFORM UNIT TESTS
###################################################################   

err=ssd_unit_test()

