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
from Compute_Variables import Compute_RH

Cenario='rcp26mc'

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


print(Cenario)

#for cmod in range(0,nMod):
for cmod in [33]:
#for cmod in [7]:


      	tasmax,tasmax_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmax',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        tasmin,tasmin_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmin',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        ps,ps_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'ps',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        huss,huss_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'huss',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)

	if tasmax_flag==1 and tasmin_flag==1 and ps_flag==1 and huss_flag==1:
		tasmea = (tasmax+tasmin)/2.
		del tasmax
		del tasmin
		RH,RH_flag=Compute_RH(tasmea,ps,huss)
		del ps
		del huss
	else:
		RH_flag=0
	
	if RH_flag==1:
		Lt,Ly,Lx=RH.shape
	
		#++++++++++++++++++++++++++++
		# RH simple avearge
		RH_MON,RH_SAZ,RH_CLI=Simple_Average(RH,ndias[cmod],Ano0i,Anoi,Anof)	

                # save RH simple average
                dirsave=driveSAVE+'/'+Cenario+'/rh'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_rh_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('rh_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=RH_MON[:,:,:]
                TMmonw.units='%'
                TMmonw.long_name='Climatological averaged daily average relative humidity at 2-m by month'

                TMsazw = f.createVariable('rh_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=RH_SAZ[:,:,:]
                TMsazw.units='%'
                TMsazw.long_name='Climatological averaged daily average relative humidity at 2-m by season'

                TMcliw = f.createVariable('rh_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=RH_CLI[:,:]
                TMcliw.units='%'
                TMcliw.long_name='Climatological averaged daily average relative humidity at 2-m'

                f.close()
		#+++++++++++++++++++
	
	
