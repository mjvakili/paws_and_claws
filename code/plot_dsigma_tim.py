import h5py
import matplotlib.pyplot as plt
plt.switch_backend("Agg")
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
import numpy as np
import dat

if __name__ == '__main__':


   plt.figure(figsize = (12,12))
   zpeak = False 
   for i in range(0,4):

       zl1 = 0.4 + i*0.1
       zl2 = 0.4 + (i+1)*0.1
        
       wr = h5py.File("pau_tim_deltasigma_bkg_z_"+str(zl1)+"_"+str(zl2)+"_"+str(zpeak)+".h5", "r")
       logr = wr["r"][:]
       lens = wr["dsigma_lens"][:]
       random = wr["dsigma_random"][:]
       Sc = wr["Sc"][0]
       err = 0.25*Sc/wr["npairs"][:]**.5 
       print err
       #plt.errorbar(logr, -1*lens+1.*random, err, fmt='o', label = str(zl1)+"$<z_{l}<$"+str(zl2))
       plt.fill_between(logr, -1*lens+1.*random-err, -1*lens+1.*random+ err, alpha = 0.3, label = str(zl1)+"$<z_{l}<$"+str(zl2))
       #plt.plot(logr, -1.*random)
       #plt.plot(logr, lens - random)
       plt.yscale("log")
       plt.xscale("log")
       plt.xlim([99, 5000])
       plt.ylim([.5, 31]) 
   if zpeak == False:    
      title = "Schrabback+2010, including sources with a secondary photo-z peak"
   elif zpeak == True:
      title = "Schrabback+2010, excluding sources with a secondary photo-z peak"
   plt.title(title, fontsize = 30)
   plt.xticks(fontsize = 30)
   plt.yticks(fontsize = 30)
   plt.legend(fontsize = 28, edgecolor = None)
   plt.xlabel(r"$R \; [\mathrm{kpc}/h]$", fontsize = 40)
   plt.ylabel(r"$\Delta \Sigma \; [M_{\odot}/\mathrm{pc}^{2}]$", fontsize = 40)
   plt.savefig("/home/vakili/public_html/sheartheta_pau_tim"+str(zpeak)+".png")
   plt.close()


   plt.figure(figsize = (12,12))
   import seaborn as sns 
   source = dat.load_source()
   zs = source['zphot']
   zp2 = source['zp2']
   for i in range(0,4):

       zl1 = 0.4 + i*0.1
       zl2 = 0.4 + (i+1)*0.1
       
       if zpeak == True:
          zsr = zs[(zs>zl2+0.1)&(zp2<0)]
       elif zpeak == False:
          zsr = zs[(zs>zl2+0.1)]
       sns.distplot(zsr, hist = False, kde = True)	  
       #plt.yscale("log")
       #plt.xscale("log")
       #plt.xlim([99, 5000])
       #plt.ylim([.5, 31]) 
   if zpeak == False:    
      title = "Schrabback+2010, including sources with a secondary photo-z peak"
   elif zpeak == True:
      title = "Schrabback+2010, excluding sources with a secondary photo-z peak"
   plt.title(title, fontsize = 30)
   plt.xticks(fontsize = 30)
   plt.yticks(fontsize = 30)
   plt.legend(fontsize = 28, edgecolor = None)
   plt.xlabel(r"$z$", fontsize = 40)
   plt.ylabel(r"$n(z)$", fontsize = 40)
   plt.savefig("/home/vakili/public_html/redshift_pau_tim"+str(zpeak)+".png")
   plt.close()

