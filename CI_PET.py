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



print(Cenario)

for cmod in range(0,nMod):


      	pet,pet_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'pet',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
	if pet_flag==1:

		Lt,Ly,Lx=pet.shape
	
		#++++++++++++++++++++++++++++
		# PET simple avearge
		pet_MON,pet_SAZ,pet_CLI=Simple_Average(pet,ndias[cmod],Ano0i,Anoi,Anof)	

                # save TMAX simple average
                dirsave=driveSAVE+'/'+Cenario+'/pet'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_pet_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('pet_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=pet_MON[:,:,:]
                TMmonw.units='mm/day'
                TMmonw.long_name='Climatological averaged daily potential evapotranspiration by month (Penman-Monteith) '

                TMsazw = f.createVariable('pet_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=pet_SAZ[:,:,:]
                TMsazw.units='mm/day'
                TMsazw.long_name='Climatological averaged daily potential evapotranspiration by season (Penman-Monteith)'

                TMcliw = f.createVariable('pet_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=pet_CLI[:,:]
                TMcliw.units='mm/day'
                TMcliw.long_name='Climatological averaged daily potential evapotranspiration (Penman-Monteith)'

                f.close()
		#+++++++++++++++++++
	
	
