import pyfits as pf
import numpy as np

def load_folder():

    return "/disks/shear14/mj/cosmos/"

def load_lens_with_qz_specz():
    '''column names:
    #
    zb, odds, zb_mean, chi2, n_band, best_run, ebv, i_auto, ra, dec
    best_run are taken from table 1 of eriksen et al 2018:

    1 False None Ell1, Ell2, Ell3, Ell4, Ell5, Ell6
    2 False None Ell6, Ell7, S0, Sa, Sb, Sc
    3 True None Sc, Sd, Sdm, SB0, SB1, SB2
    4 True None SB2, SB3, SB4, SB5, SB6, SB7, SB8, SB9, SB10, SB11
    5 False None BC03(0.008, 0.509), BC03(0.008, 8.0), BC03(0.02, 0.509), BC03(0.02, 2.1),
    BC03(0.02, 2.6), BC03(0.02, 3.75)
    6-15 True Calzetti SB4, SB5, SB6, SB7, SB8, SB9, SB10, SB11
    16-25 True Calzetti+Bump 1 SB4, SB5, SB6, SB7, SB8, SB9, SB10, SB11
    26-35 True Calzetti+Bump 2 SB4, SB5, SB6, SB7, SB8, SB9, SB10, SB11

    '''
    folder_path = load_folder()
    filename = "pauxcosmos_zs.fits"
    y = pf.open(folder_path+filename)[1].data
    zb, zs, qz, zb_mean, chi2, n_band, best_run, ebv, i_auto, ra, dec = y['zb'], y['zspec'], y['qz'], y['zb_mean'], y['chi2'], y['n_band'], y['best_run'], y['ebv'], y['i_auto'], y['ra'], y['dec']

    return np.vstack([zb, zs, qz, zb_mean, chi2, n_band, best_run, ebv, i_auto, ra, dec]).T

def load_lens_with_qz_pz():
    '''column names:
    #
    zb, odds, zb_mean, chi2, n_band, best_run, ebv, i_auto, ra, dec
    best_run are taken from table 1 of eriksen et al 2018:

    1 False None Ell1, Ell2, Ell3, Ell4, Ell5, Ell6
    2 False None Ell6, Ell7, S0, Sa, Sb, Sc
    3 True None Sc, Sd, Sdm, SB0, SB1, SB2
    4 True None SB2, SB3, SB4, SB5, SB6, SB7, SB8, SB9, SB10, SB11
    5 False None BC03(0.008, 0.509), BC03(0.008, 8.0), BC03(0.02, 0.509), BC03(0.02, 2.1),
    BC03(0.02, 2.6), BC03(0.02, 3.75)
    6-15 True Calzetti SB4, SB5, SB6, SB7, SB8, SB9, SB10, SB11
    16-25 True Calzetti+Bump 1 SB4, SB5, SB6, SB7, SB8, SB9, SB10, SB11
    26-35 True Calzetti+Bump 2 SB4, SB5, SB6, SB7, SB8, SB9, SB10, SB11

    '''
    folder_path = load_folder()
    filename = "pauxcosmos_zp.fits"
    y = pf.open(folder_path+filename)[1].data
    zb, qz, zb_mean, chi2, n_band, best_run, ebv, i_auto, radius, sersic, ra, dec = y['zb'], y['qz'], y['zb_mean'], y['chi2'], y['n_band'], y['best_run'], y['ebv'], y['i_auto'], np.sqrt(y['acs_a_image']*y['acs_b_image']), y['sersic_n_gim2d'], y['ra'], y['dec']

    return np.vstack([zb, qz, zb_mean, chi2, n_band, best_run, ebv, i_auto, radius, sersic, ra, dec]).T

def load_lens_with_mstar():

    folder_path = load_folder()
    filename = "pauxcosmos_zp.fits"
    pau = pf.open(folder_path+filename)[1].data
    # _paudm_id   mass_med   mass_med_min68   mass_med_max68   mass_best   sfr_med      sfr_med_min68   sfr_med_max68   sfr_best      NB_mass                 NB_mass_err             NB_sfr(10Myrs)         NB_sfr(10Myrs)_err     BB_mass                 BB_mass_err             BB_sfr(10Myrs)        BB_sfr(10Myrs)_err  

    starname = "st_mass.dat"
    
    
    with open(folder_path+starname) as f:
        lines = f.readlines()
    #print lines
    
    y = []
    i = 0
    for line in lines:
        i = i + 1 
	if ('""' not in line.split()):
	    #print i
	    #print line.split()

            y.append(np.array(line.split()).astype(float))
	   
    st = np.array(y)
    st_id = st[:,0]
    pau_id = pau['paudm_id']

    mask_one = np.where(np.in1d(pau_id, st_id))[0]
    mask_two = np.where(np.in1d(st_id, pau_id[mask_one]))[0]
    arg_one = np.argsort(pau_id[mask_one])
    arg_two = np.argsort(st_id[mask_two])

    print pau_id[mask_one][arg_one]
    print st_id[mask_two][arg_two]

    pau = pau[mask_one][arg_one]
    st = st[mask_two][arg_two]
    
    y = pau

    zb, qz, zb_mean, chi2, n_band, best_run, ebv, i_auto, radius, sersic, ra, dec = y['zb'], y['qz'], y['zb_mean'], y['chi2'], y['n_band'], y['best_run'], y['ebv'], y['i_auto'], np.sqrt(y['acs_a_image']*y['acs_b_image']), y['sersic_n_gim2d'], y['ra'], y['dec']

    return np.vstack([zb, np.log10(st[:,9]), np.log10(st[:,10]), st[:,11], st[:,12], qz, zb_mean, chi2, n_band, best_run, ebv, i_auto, radius, sersic, ra, dec]).T

def load_lens():
    '''coloumn names:
       #zb jk  dec   ra       z      Iauto   K         Z       R       G        V      B        U
       objects of interest according to EG: zb=2 corresponding to 0.4<z<0.8
    '''
    folder_path = load_folder()
    filename = "zCOSMOS-in-wq6jk25Is1510.004-19.025.019.025.0zCOSMOS-Bext4.dat"
    
    return np.loadtxt(folder_path+filename)

def load_lens_cosmos():


    folder_path = load_folder()
    filename = "pauxcosmos_zp.fits"
    y = pf.open(folder_path+filename)[1].data
        
    header = pf.open(folder_path+filename)[1].header
    
    return y, header

def load_random():
    '''column names: JK-regions, DEC, RA

    approximate randoms created by EG:
    The mask in my randoms is only approximate, but I have tested that it
    gives good results for auto and cross correlation of galaxy counts.
    The way I built this mask is by taking all objects flagged as masked
    in the  Ilbert et al.,( ApJ, 2008, astro-ph/0809.2101)
    sample to Iauto<25 and then built a mask with radius of 15arcsec
    around these objects. This is conservative as it creates large halos
    around bright stars
    '''

    folder_path = load_folder()
    filename = "zCOSMOS-rndwq6jk25Is1510.004-19.025.019.025.0zCOSMOS-Bext4.dat"
    
    return np.loadtxt(folder_path+filename)
    

def load_source():
    '''
    cosmos shear data from alexie.
    columns:
    `unique_id`, `ra`, `dec`, `zphot`, `zp2`, `photoz`, `gamma1`, `gamma2`, `weight`
    
    shape_noise = 0.27 weight = 1.0/(var_e1+(shape_noise)^2)
    Field zp2 is the redshift of a secondary peak if it exists. If you want more robust photoz's, you can limit sources to those with zp2<0
    When they exist, photoz-s have been replaced with speczs.
    L'Aigle et al 2015 photoz catalog
    Leauthaud 2007 shape catalog
    '''
    folder_path = load_folder()
    filename = folder_path+"source.fits"
    source = pf.open(filename)
    
    return source[1].data

def load_tim():
    '''
    cosmos shear data from tim.
    columns:
    `ra`, `dec`, `gamma1`, `gamma2`, `weight`, `photoz`, `zp2`, `zphot`
    
    shape_noise = 0.25 weight = 1.0/(var_e1+(shape_noise)^2)
    Field zp2 is the redshift of a secondary peak if it exists. If you want more robust photoz's, you can limit sources to those with zp2<0
    When they exist, photoz-s have been replaced with speczs.
    L'Aigle et al 2015 photoz catalog
    Leauthaud 2007 shape catalog
    '''
    
    return np.loadtxt('new_tim.txt')

if __name__ == '__main__':


    load_lens_with_mstar()
