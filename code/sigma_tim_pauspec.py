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
plt.switch_backend("Agg")

# Import constants
G = const.G.to('pc3/Msun s2')
c = const.c.to('pc/s')
inf = np.inf
h, O_matter, O_lambda = [0.67, 0.319, 1]
cosmo = LambdaCDM(H0=h*100, Om0=O_matter, Ode0=O_lambda)

##### loading lens & source ######

def lens_source(zl1,zl2,zs1,zs2,zpeak,qcut,lens_spec):
        '''  
        lname = np.vstack([zb, zs, qz, zb_mean, chi2, n_band, best_run, ebv, i_auto, ra, dec]).T
	'''
	lens_z = lname[:,0]
	qz = lname[:,2]
        i_auto = lname[:,-3]

        lens_mask = (lens_z>zl1)&(lens_z<zl2)&(qz<np.percentile(qz, qcut))&(i_auto<22.5)
	source_z = sname[:,-1]
	source_zp2 = sname[:,6] #secondary photo-z peak
	
	if zpeak == False:

	   source_mask = (source_z>zs1)&(source_z<zs2)
        
        elif zpeak == True:
	    
	   source_mask = (source_z>zs1)&(source_z<zs2)&(source_zp2<0)

	lens = lname[lens_mask]
	source = sname[source_mask]
        
	   
	lens_z = lens[:,0]
        source_z = source[:,-1]

        plt.figure(figsize = (10,10))
	sns.distplot(lens_z, norm_hist = True, kde = False, hist = True, label = "PAU "+str(qcut)+"%", hist_kws={"histtype": "step", "linewidth": 1.5, "alpha": 1})
	sns.distplot(source_z, norm_hist = True, kde = False, hist = True, label = "COSMOS source redshift distribution", hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1})
	plt.xlabel("$z$", fontsize = 40)
	plt.ylabel("$n(z)$", fontsize = 40)
	plt.xlim([0,6])
	plt.legend(loc = 'best', fontsize = 15)
        plt.savefig("/home/vakili/public_html/zdist_zlmin_"+str(zl1)+"_zlmax_"+str(zl2)+"_qcut_"+str(qcut)+".png")
        plt.close() 
	return lens, source

def dsigma(lens, source, lens_spec):

    source_ra = source[:,0]
    source_dec = source[:,1]
    source_z = source[:,-1]
    source_e1 = -1.* source[:,2]
    source_e2 = source[:,3]
    source_w = source[:,4]
    
    if lens_spec == False:
	
       lens_z = lens[:,0]

    elif lens_spec == True:
	   
       lens_z = lens[:,1]
    
    lens_ra = lens[:,-2]
    lens_dec = lens[:,-1]
    
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

   sname  = dat.load_tim()
   lname  = dat.load_lens_with_qz_specz()
   rname  = dat.load_random()
   nrands = len(rname)

   zl1, zl2 = 0.5, 0.6
   zs1, zs2 = 0.6, 3.0
   zpeak = True
   qcut = 75
   lens_spec = True

   lens, source = lens_source(zl1, zl2, zs1, zs2, zpeak, qcut, lens_spec)
     
   for i in [10,11]:
    for qcut in [25,50,75,100]:
     for lens_spec in [True, False]:

       zl1 = 0.0 + i*0.1
       zl2 = 0.0 + (i+1)*0.1
       zs1 , zs2 = zl2 + 0.2 , 10.0
       print "lens redshifts", zl1, zl2
       
       lens, source = lens_source(zl1, zl2, zs1, zs2, zpeak, qcut, lens_spec)
       if lens_spec == False:
	
         lens_z = lens[:,0]
	
       elif lens_spec == True:
	   
         lens_z = lens[:,1]
   
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
       
       
       R, St, Sx , W, Str, Sxr, Wr, Npairs ,Sc = dsigma(lens, source, lens_spec)
       print "done with dsigma_i"

       Rc = cosmo.kpc_comoving_per_arcmin(np.median(lens_z))*R

       dsigma_file = h5py.File("output/pau_tim_deltasigma_bkg_z_"+str(zl1)+"_"+str(zl2)+"_zpeak_"+str(zpeak)+"_qcut_"+str(qcut)+"_lens_spec_"+str(lens_spec)+".h5", "w")
       dsigma_file.create_dataset("r", (5,), data=Rc)
       dsigma_file.create_dataset("dsigma_lens", (5,),data = St)
       dsigma_file.create_dataset("dsigma_random", (5,),data = Str)
       dsigma_file.create_dataset("npairs", (5,),data = Npairs)
       dsigma_file.create_dataset("Sc", (1,),data = Sc)
       dsigma_file.close()

       print "done with writing the file for the lens bin" , i
