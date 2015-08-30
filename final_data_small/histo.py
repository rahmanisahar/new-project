import numpy as np
from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits

m31=Table.read('m31_data_f1.fits')
m101=Table.read('m101_data_f.fits')
tt=fits.open('m31_data_f1.fits')
t=fits.open('m101_data_f.fits')
cols=tt[1].columns
cols1=t[1].columns
for n in range(23,43):
    plt.hist(m31[cols.names[n]],bins=20, histtype='stepfilled', normed=True, color='r', alpha=0.5, label='M31')
    plt.rc('xtick', labelsize=8) 
    plt.hist(m101[cols.names[n]],bins=20, histtype='stepfilled', normed=True, color='g', alpha=0.5, label='m101')
    plt.title(cols.names[n])
    plt.legend()
    #plt.show()
    plt.savefig('histograms\plot'+str(n)+'.pdf')
