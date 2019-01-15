import h5py
import matplotlib.pyplot as plt
plt.switch_backend("Agg")
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
import numpy as np
import dat


def plot(zl1, zl2, zpeak, qcut):

   plt.figure(figsize = (12,12))
  
   lens_spec = False
   wr = h5py.File("output/pau_tim_deltasigma_bkg_z_"+str(zl1)+"_"+str(zl2)+"_zpeak_"+str(zpeak)+"_qcut_"+str(qcut)+"_lens_spec_"+str(lens_spec)+".h5")
   logr = wr["r"][:]
   lens = wr["dsigma_lens"][:]
   random = wr["dsigma_random"][:]
   Sc = wr["Sc"][0]
   err = 0.25*Sc/wr["npairs"][:]**.5 
   print -1.*lens + 1.*random
   plt.errorbar(logr, -1*lens+1.*random, err, fmt='o', label = str(zl1)+"$<z_{l}<$"+str(zl2)+", PAU "+str(qcut)+"%"+", PAU photo-z")
   #plt.fill_between(logr, -1*lens+1.*random-err, -1*lens+1.*random+ err, alpha = 0.3, label = str(zl1)+"$<z_{l}<$"+str(zl2)+", PAU "+str(qcut)+"%"+", PAU photo-z")
   #plt.plot(logr, -1.*random)
   #plt.plot(logr, lens - random)
  
   lens_spec = True
   wr = h5py.File("output/pau_tim_deltasigma_bkg_z_"+str(zl1)+"_"+str(zl2)+"_zpeak_"+str(zpeak)+"_qcut_"+str(qcut)+"_lens_spec_"+str(lens_spec)+".h5")
   logr = wr["r"][:]
   lens = wr["dsigma_lens"][:]
   random = wr["dsigma_random"][:]
   Sc = wr["Sc"][0]
   err = 0.25*Sc/wr["npairs"][:]**.5 
   print -1.*lens + 1.* random
   plt.errorbar(logr+0.1*logr, -1*lens+1.*random, err, fmt='o',label = str(zl1)+"$<z_{l}<$"+str(zl2)+", PAU "+str(qcut)+"%"+", COSMOS spec-z")
   #plt.fill_between(logr, -1*lens+1.*random-err, -1*lens+1.*random+ err, alpha = 0.3, label = str(zl1)+"$<z_{l}<$"+str(zl2)+", PAU "+str(qcut)+"%"+", COSMOS spec-z")
   #plt.plot(logr, -1.*random)
   #plt.plot(logr, lens - random)
  
   #plt.yscale("log")
   plt.xscale("log")
   plt.xlim([90, 5000])
   #plt.ylim([.5, 31]) 
   if zpeak == False:    
      title = "including sources with a secondary photo-z peak"
   elif zpeak == True:
      title = "excluding sources with a secondary photo-z peak"
   plt.title(title, fontsize = 30)
   plt.xticks(fontsize = 30)
   plt.yticks(fontsize = 30)
   plt.legend(fontsize = 28, edgecolor = None)
   plt.xlabel(r"$R \; [\mathrm{kpc}/h]$", fontsize = 40)
   plt.ylabel(r"$\Delta \Sigma \; [M_{\odot}/\mathrm{pc}^{2}]$", fontsize = 40)
   
   plt.savefig("fig/sigma_pau_tim_deltasigma_bkg_z_"+str(zl1)+"_"+str(zl2)+"_zpeak_"+str(zpeak)+"_qcut_"+str(qcut)+".png")
   plt.close()
   
   return None

if __name__ == '__main__':

   zpeak  = True
   
   for i in [6,7,8,9]:
    for qcut in [25,50,75,100]:

       zl1 = 0.0 + i*0.1
       zl2 = 0.0 + (i+1)*0.1
   
   #for i in [0, 1, 2]:
   # for qcut in [25,50,75,100]:

   #    zl1 = 0.0 + i*0.3
   #    zl2 = 0.0 + (i+1)*0.3
       zs1 , zs2 = zl2 + 0.2 , 10.0
       plot(zl1, zl2, zpeak, qcut)
