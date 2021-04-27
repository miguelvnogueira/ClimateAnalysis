import numpy as np
from netCDF4 import Dataset
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import calendar
import os

sys.path.append("../Funcoes")

from Load_Models import Cenario_Props
from Load_Models import Ler_Modelo


from Compute_Extremes import Simple_Average


Cenario='rcp26ec'

driveCORDEX,driveSAVE,Ano0i,Ano0f,Anoi,Anof,Cenarioleg=Cenario_Props(Cenario)

GCMs=np.load(driveSAVE+'/Lista_GCMs.npy')
RCMs=np.load(driveSAVE+'/Lista_RCMs.npy')
EXPs=np.load(driveSAVE+'/Lista_EXPs.npy')
VERs=np.load(driveSAVE+'/Lista_VERs.npy')
ndias=np.load(driveSAVE+'/Lista_ndias.npy')
nMod=len(GCMs)


                
# CRIAR PASTA DO CENARIO SE NAO EXISTIR                
dirsave=driveSAVE+'/'+Cenario
if os.path.isdir(dirsave)==False:                        
	os.mkdir(dirsave)

for cmod in [35]:
#for cmod in range(0,nMod):

	print(GCMs[cmod]+'/'+RCMs[cmod])

      	tasmax,tasmax_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmax',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
	tasmax=tasmax-273.15

        tasmin,tasmin_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmin',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        tasmin=tasmin-273.15

	dtr=tasmax-tasmin
	tasmea = (tasmax+tasmin)/2.

	del tasmax
	del tasmin

	if tasmax_flag==1 and tasmin_flag==1:
		Lt,Ly,Lx=tasmea.shape
	
		#++++++++++++++++++++++++++++
		# Tmea simple avearge
		Tm_MON,Tm_SAZ,Tm_CLI=Simple_Average(tasmea,ndias[cmod],Ano0i,Anoi,Anof)	

                # save TMea simple average
                dirsave=driveSAVE+'/'+Cenario+'/Tm'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_Tm_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('Tm_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=Tm_MON[:,:,:]
                TMmonw.units='deg C'
                TMmonw.long_name='Climatological average daily mean 2m air temperature per month'

                TMsazw = f.createVariable('Tm_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=Tm_SAZ[:,:,:]
                TMsazw.units='deg C'
                TMsazw.long_name='Climatological average daily mean 2m air temperature per season'

                TMcliw = f.createVariable('Tm_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=Tm_CLI[:,:]
                TMcliw.units='deg C'
                TMcliw.long_name='Climatological average daily mean 2m air temperature'

                f.close()
		#+++++++++++++++++++
	

                #++++++++++++++++++++++++++++
                # DTR simple avearge
                DTR_MON,DTR_SAZ,DTR_CLI=Simple_Average(dtr,ndias[cmod],Ano0i,Anoi,Anof)

                # save DTR simple average
                dirsave=driveSAVE+'/'+Cenario+'/DTR'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_DTR_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('DTR_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=DTR_MON[:,:,:]
                TMmonw.units='deg C'
                TMmonw.long_name='Climatological average daily 2m air temperature range per month'

                TMsazw = f.createVariable('DTR_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=DTR_SAZ[:,:,:]
                TMsazw.units='deg C'
                TMsazw.long_name='Climatological average daily mean 2m air temperature per season'

                TMcliw = f.createVariable('DTR_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=DTR_CLI[:,:]
                TMcliw.units='deg C'
                TMcliw.long_name='Climatological average daily mean 2m air temperature'

                f.close()
                #+++++++++++++++++++

		
