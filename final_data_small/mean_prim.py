import numpy as np
from astropy.table import Table
from scipy.stats import nanmean


t=Table.read('primary.fits')
t['PAH5.7eqw']=t['PAH5.7eqw']/nanmean(t['PAH5.7eqw'])

t['PAH6.2eqw']=t['PAH6.2eqw']/nanmean(t['PAH6.2eqw'])

t['PAH7.4eqw']=t['PAH7.4eqw']/nanmean(t['PAH7.4eqw'])

t['PAH7.6eqw']=t['PAH7.6eqw']/nanmean(t['PAH7.6eqw'])

t['PAH7.9eqw']=t['PAH7.9eqw']/nanmean(t['PAH7.9eqw'])

t['PAH8.3eqw']=t['PAH8.3eqw']/nanmean(t['PAH8.3eqw'])

t['PAH8.6eqw']=t['PAH8.6eqw']/nanmean(t['PAH8.6eqw'])

t['PAH10.7eqw']=t['PAH10.7eqw']/nanmean(t['PAH10.7eqw'])

t['PAH11.23eqw']=t['PAH11.23eqw']/nanmean(t['PAH11.23eqw'])

t['PAH11.33eqw']=t['PAH11.33eqw']/nanmean(t['PAH11.33eqw'])

t['PAH12.0eqw']=t['PAH12.0eqw']/nanmean(t['PAH12.0eqw'])

t['PAH12.62eqw']=t['PAH12.62eqw']/nanmean(t['PAH12.62eqw'])

t['PAH12.69eqw']=t['PAH12.69eqw']/nanmean(t['PAH12.69eqw'])

#t['PAH14.0eqw']=t['PAH14.0eqw']/nanmean(t['PAH14.0eqw'])

t['PAH16.45eqw']=t['PAH16.45eqw']/nanmean(t['PAH16.45eqw'])

t['PAH17.04eqw']=t['PAH17.04eqw']/nanmean(t['PAH17.04eqw'])

t['PAH17.39eqw']=t['PAH17.39eqw']/nanmean(t['PAH17.39eqw'])

#t['PAH17.87eqw']=t['PAH17.87eqw']/nanmean(t['PAH17.87eqw'])

t['PAH5.7flx']=t['PAH5.7flx']/nanmean(t['PAH5.7flx'])

t['PAH6.2flx']=t['PAH6.2flx']/nanmean(t['PAH6.2flx'])

t['PAH7.4flx']=t['PAH7.4flx']/nanmean(t['PAH7.4flx'])

t['PAH7.6flx']=t['PAH7.6flx']/nanmean(t['PAH7.6flx'])

t['PAH7.9flx']=t['PAH7.9flx']/nanmean(t['PAH7.9flx'])

t['PAH8.3flx']=t['PAH8.3flx']/nanmean(t['PAH8.3flx'])

t['PAH8.6flx']=t['PAH8.6flx']/nanmean(t['PAH8.6flx'])

t['PAH10.7flx']=t['PAH10.7flx']/nanmean(t['PAH10.7flx'])

t['PAH11.23flx']=t['PAH11.23flx']/nanmean(t['PAH11.23flx'])

t['PAH11.33flx']=t['PAH11.33flx']/nanmean(t['PAH11.33flx'])

t['PAH12.0flx']=t['PAH12.0flx']/nanmean(t['PAH12.0flx'])

t['PAH12.62flx']=t['PAH12.62flx']/nanmean(t['PAH12.62flx'])

t['PAH12.69flx']=t['PAH12.69flx']/nanmean(t['PAH12.69flx'])

#t['PAH14.0flx']=t['PAH14.0flx']/nanmean(t['PAH14.0flx'])

t['PAH16.45flx']=t['PAH16.45flx']/nanmean(t['PAH16.45flx'])

t['PAH17.04flx']=t['PAH17.04flx']/nanmean(t['PAH17.04flx'])

t['PAH17.39flx']=t['PAH17.39flx']/nanmean(t['PAH17.39flx'])

#t['PAH17.87flx']=t['PAH17.87flx']/nanmean(t['PAH17.87flx'])

t['ArII']=t['ArII']/nanmean(t['ArII'])

t['ArIII']=t['ArIII']/nanmean(t['ArIII'])

t['SIV']=t['SIV']/nanmean(t['SIV'])

t['NeII']=t['NeII']/nanmean(t['NeII'])

t['NeIII']=t['NeIII']/nanmean(t['NeIII'])

t['SIII']=t['SIII']/nanmean(t['SIII'])

t['IRAC1_LUM_PER_PC_SQ']=t['IRAC1_LUM_PER_PC_SQ']/nanmean(t['IRAC1_LUM_PER_PC_SQ'])

t['IRAC2_LUM_PER_PC_SQ']=t['IRAC2_LUM_PER_PC_SQ']/nanmean(t['IRAC2_LUM_PER_PC_SQ'])

t['IRAC3_LUM_PER_PC_SQ']=t['IRAC3_LUM_PER_PC_SQ']/nanmean(t['IRAC3_LUM_PER_PC_SQ'])

t['IRAC4_LUM_PER_PC_SQ']=t['IRAC4_LUM_PER_PC_SQ']/nanmean(t['IRAC4_LUM_PER_PC_SQ'])

t['MIPS24_LUM_PER_PC_SQ']=t['MIPS24_LUM_PER_PC_SQ']/nanmean(t['MIPS24_LUM_PER_PC_SQ'])

t['MIPS70_LUM_PER_PC_SQ']=t['MIPS70_LUM_PER_PC_SQ']/nanmean(t['MIPS70_LUM_PER_PC_SQ'])


t['PACS100_LUM_PER_PC_SQ']=t['PACS100_LUM_PER_PC_SQ']/nanmean(t['PACS100_LUM_PER_PC_SQ'])

t['PACS160_LUM_PER_PC_SQ']=t['PACS160_LUM_PER_PC_SQ']/nanmean(t['PACS160_LUM_PER_PC_SQ'])

t['SPIRE250_LUM_PER_PC_SQ']=t['SPIRE250_LUM_PER_PC_SQ']/nanmean(t['SPIRE250_LUM_PER_PC_SQ'])

t['SPIRE350_LUM_PER_PC_SQ']=t['SPIRE350_LUM_PER_PC_SQ']/nanmean(t['SPIRE350_LUM_PER_PC_SQ'])

t['SPIRE500_LUM_PER_PC_SQ']=t['SPIRE500_LUM_PER_PC_SQ']/nanmean(t['SPIRE500_LUM_PER_PC_SQ'])

t.write('primary_mean.fits')


