import pyfits as pf
import h5py
import numpy as np
import kmeans_radec
import matplotlib.pyplot as plt
import treecorr
import seaborn as sns
from kmeans_radec import KMeans, kmeans_sample
import csv
import pandas as pd
import matplotlib
from astropy import constants as const, units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import LambdaCDM
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
import multiprocessing
from multiprocessing import Pool
from astropy import constants as const, units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import LambdaCDM
import dat

# Import constants
G = const.G.to('pc3/Msun s2')
c = const.c.to('pc/s')
inf = np.inf
h, O_matter, O_lambda = [0.67, 0.319, 1]
cosmo = LambdaCDM(H0=h*100, Om0=O_matter, Ode0=O_lambda)

##### loading lens & source ######

def lens_source(zl1,zl2,zs1,zs2,zpeak):
        
	lens_z = lname[:,4]
	source_z = sname["zphot"]
	source_zp2 = sname["zp2"] #secondary photo-z peak
        lens_mask = (lens_z>zl1)&(lens_z<zl2)
	
	if zpeak == False:

	   source_mask = (source_z>zs1)&(source_z<zs2)
        
        elif zpeak == True:
	    
	    source_mask = (source_z>zs1)&(source_z<zs2)&(source_zp2<0)
        

	lens = lname[lens_mask]
	source = sname[source_mask]
	   
        return lens, source

def dsigma(lens, source):

    source_ra = source["ra"]
    source_dec = source["dec"]
    source_z = source["zphot"]
    source_e1 = -1.* source["gamma1"]
    source_e2 = source["gamma2"]
    source_w = source["weight"]

    lens_z = lens[:,4]
    lens_ra = lens[:,3]
    lens_dec = lens[:,2]
    
    lens_Dc = cosmo.comoving_distance(np.median(lens_z))
    source_Dc = cosmo.comoving_distance(source_z)
    lens_Da = cosmo.angular_diameter_distance(np.median(lens_z))    
    
    source_Dc =  source_Dc*1e6
    lens_Dc =  lens_Dc*1e6
    lens_Da = lens_Da*1e6
    
    DlsoDs = (source_Dc - lens_Dc)/source_Dc
    Sigma_crit = (c.value**2)/(4*np.pi*G.value) * 1/(lens_Da * DlsoDs)
    
    lens_cat = treecorr.Catalog(x=lens_ra, y=lens_dec, x_units='degree', y_units='degree')
    source_cat=treecorr.Catalog(x=source_ra,y=source_dec,
                                g1=source_e1*Sigma_crit, g2=source_e2*Sigma_crit,
				w=source_w*Sigma_crit**-2., x_units='degree',y_units='degree')


    ng = treecorr.NGCorrelation(nbins = 5, min_sep=0.2, max_sep=6, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    r, xi_t , xi_x , w , npairs = ng.meanr, ng.xi, ng.xi_im, ng.weight, ng.npairs
		   
    print "done with lens source cross-correlation"

    lens_Dc = cosmo.comoving_distance(np.median(rand_z))
    lens_Da = cosmo.angular_diameter_distance(np.median(rand_z))    
    
    source_Dc =  source_Dc*1e6
    lens_Dc =  lens_Dc*1e6
    lens_Da = lens_Da*1e6
    
    DlsoDs = (source_Dc - lens_Dc)/source_Dc
    Sigma_crit = (c.value**2)/(4*np.pi*G.value) * 1/(lens_Da * DlsoDs)

    random_cat = treecorr.Catalog(x=rname[:,2], y=rname[:,1], 
	                          x_units='degree', y_units='degree')
    source_cat=treecorr.Catalog(x=source_ra,y=source_dec,
                                g1=source_e1*Sigma_crit, g2=source_e2*Sigma_crit,
				w=source_w*Sigma_crit**-2., x_units='degree',y_units='degree')
    ng = treecorr.NGCorrelation(nbins = 5, min_sep=0.2, max_sep=6, sep_units='arcmin', verbose=1)
    ng.process(random_cat, source_cat)
    rr, xi_tr , xi_xr,  wr = ng.meanr, ng.xi, ng.xi_im, ng.weight

    return r, xi_t, xi_x, w, xi_tr, xi_xr, wr, npairs, np.mean(Sigma_crit)

if __name__ == '__main__':

   halfpi, pi, twopi = [f*np.pi for f in [0.5, 1, 2]]
   degs, rads = 180/pi, pi/180

   sname  = dat.load_source()
   lname  = dat.load_lens()
   rname  = dat.load_random()
   nrands = len(rname)

   zl1, zl2 = 0.5, 0.6
   zs1, zs2 = 0.6, 3.0
   zpeak = False

   lens, source = lens_source(zl1, zl2, zs1, zs2, zpeak)
     
   for i in [0, 1, 2, 3]:

       zl1 = 0.4 + i*0.1
       zl2 = 0.4 + (i+1)*0.1
       zs1 , zs2 = zl2 + 0.1 , 4.0

       lens, source = lens_source(zl1, zl2, zs1, zs2, zpeak)
       lens_z = lens[:,4]
   
       print "randoms to lens ratio" , nrands * 1.0/ len(lens)
       

       print "assigning redshifts to randoms"
       hist, bins = np.histogram(lens_z, bins = 50)
       bin_midpoints = bins[:-1] + np.diff(bins)/2
       cdf = np.cumsum(hist)
       cdf = cdf / float(cdf[-1])
       nrands = len(rname)
       values = np.random.rand(nrands)
       value_bins = np.searchsorted(cdf, values)
       rand_z = bin_midpoints[value_bins]
       
       
       R, St, Sx , W, Str, Sxr, Wr, Npairs ,Sc = dsigma(lens, source)
       print "done with dsigma_i"

       Rc = cosmo.kpc_comoving_per_arcmin(np.median(lens_z))*R

       dsigma_file = h5py.File("pau_deltasigma_bkg_z_"+str(zl1)+"_"+str(zl2)+"_"+str(zpeak)+".h5", "w")
       dsigma_file.create_dataset("r", (5,), data=Rc)
       dsigma_file.create_dataset("dsigma_lens", (5,),data = St)
       dsigma_file.create_dataset("dsigma_random", (5,),data = Str)
       dsigma_file.create_dataset("npairs", (5,),data = Npairs)
       dsigma_file.create_dataset("Sc", (1,),data = Sc)
       dsigma_file.close()

       print "done with writing the file for the lens bin" , i
