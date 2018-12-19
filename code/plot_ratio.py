import h5py
import matplotlib.pyplot as plt
plt.switch_backend("Agg")
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
import numpy as np
import dat

if __name__ == '__main__':


   plt.figure(figsize = (18,10))
   zpeak = True 
   for i in range(0,4):

       zl1 = 0.4 + i*0.1
       zl2 = 0.4 + (i+1)*0.1
        
       wr = h5py.File("pau_tim_deltasigma_bkg_z_"+str(zl1)+"_"+str(zl2)+"_"+str(zpeak)+".h5", "r")
       logr = wr["r"][:]
       lens = wr["dsigma_lens"][:]
       random = wr["dsigma_random"][:]
       Sc = wr["Sc"][0]
       err = 0.25*Sc/wr["npairs"][:]**.5 
       
       wr = h5py.File("pau_deltasigma_bkg_z_"+str(zl1)+"_"+str(zl2)+"_"+str(zpeak)+".h5", "r")
       logr1 = wr["r"][:]
       lens1 = wr["dsigma_lens"][:]
       random1 = wr["dsigma_random"][:]
       Sc1 = wr["Sc"][0]
       err1 = 0.25*Sc1/wr["npairs"][:]**.5 
       
       plt.fill_between(logr, (-1*lens+1.*random-err)/(-1*lens1+random1), (-1*lens+1.*random+ err)/(-1*lens1+random1), alpha = 0.3, label = str(zl1)+"$<z_{l}<$"+str(zl2))
       
       #plt.yscale("log")
       plt.xscale("log")
       plt.xlim([99, 5000])
       #plt.ylim([.5, 31]) 
   if zpeak == False:    
      title = "including sources with a secondary photo-z peak"
   elif zpeak == True:
      title = "excluding sources with a secondary photo-z peak"
   plt.title(title, fontsize = 30)
   plt.xticks(fontsize = 30)
   plt.yticks(fontsize = 30)
   plt.legend(fontsize = 28, edgecolor = None)
   plt.xlabel(r"$R \; [\mathrm{kpc}/h]$", fontsize = 30)
   plt.ylabel(r"$\Delta \Sigma_{\rm cosmos, 2010}/\Delta \Sigma_{\rm cosmos, 2007}$", fontsize = 30)
   plt.savefig("/home/vakili/public_html/sheartheta_pau_ratio"+str(zpeak)+".png")
   plt.close()
