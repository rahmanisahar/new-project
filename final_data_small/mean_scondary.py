import numpy as np
from astropy.table import Table
from scipy.stats import nanmean


t=Table.read('secondary.fits')
t['12plogOH']=t['12plogOH']/nanmean(t['12plogOH'])

t['PAH12.7eqw']=t['PAH12.7eqw']/nanmean(t['PAH12.7eqw'])
t['PAH11.3eqw']=t['PAH11.3eqw']/nanmean(t['PAH11.3eqw'])
t['PAH7.7eqw']=t['PAH7.7eqw']/nanmean(t['PAH7.7eqw'])
t['PAH17.0eqw']=t['PAH17.0eqw']/nanmean(t['PAH17.0eqw'])

t['PAH12.7flx']=t['PAH12.7flx']/nanmean(t['PAH12.7flx'])
t['PAH11.3flx']=t['PAH11.3flx']/nanmean(t['PAH11.3flx'])
t['PAH7.7flx']=t['PAH7.7flx']/nanmean(t['PAH7.7flx'])
t['PAH17.0flx']=t['PAH17.0flx']/nanmean(t['PAH17.0flx'])

t['SFR_FUV_PER_PC_SQ']=t['SFR_FUV_PER_PC_SQ']/nanmean(t['SFR_FUV_PER_PC_SQ'])
t['TIR_LUM_PER_PC_SQ']=t['TIR_LUM_PER_PC_SQ']/nanmean(t['TIR_LUM_PER_PC_SQ'])
t['STELLAR_MASS_PER_PC_SQ']=t['STELLAR_MASS_PER_PC_SQ']/nanmean(t['STELLAR_MASS_PER_PC_SQ'])
t['METALLICITY_PER_PC_SQ']=t['METALLICITY_PER_PC_SQ']/nanmean(t['METALLICITY_PER_PC_SQ'])
t['TOTAL_GAS_PER_MASS_OF_SUN_PER_PC_SQ']=t['TOTAL_GAS_PER_MASS_OF_SUN_PER_PC_SQ']/nanmean(t['TOTAL_GAS_PER_MASS_OF_SUN_PER_PC_SQ'])

t.write('secondary_mean.fits')
