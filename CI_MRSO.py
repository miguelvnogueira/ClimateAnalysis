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
from Compute_Extremes import GT_numdays
from Compute_Extremes import LT_numdays
from Compute_Extremes import Simple_Maximum

Cenario='rcp26bc'

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
for cmod in [11]:

      	mrso,mrso_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'mrso',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
	

	if mrso_flag==1:
		Lt,Ly,Lx=mrso.shape
	
		#++++++++++++++++++++++++++++
		# vh10 simple avearge
		mrso_MON,mrso_SAZ,mrso_CLI=Simple_Average(mrso,ndias[cmod],Ano0i,Anoi,Anof)	

                # save mrso simple average
                dirsave=driveSAVE+'/'+Cenario+'/mrso'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_mrso_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('mrso_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=mrso_MON[:,:,:]
                TMmonw.units='mm'
                TMmonw.long_name='Climatological averaged Total Soil Moisture Content for each month'

                TMsazw = f.createVariable('mrso_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=mrso_SAZ[:,:,:]
                TMsazw.units='mm'
                TMsazw.long_name='Climatological averaged Total Soil Moisture Content for each season'

                TMcliw = f.createVariable('mrso_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=mrso_CLI[:,:]
                TMcliw.units='mm'
                TMcliw.long_name='Climatological averaged Total Soil Moisture Content'

                f.close()
		#+++++++++++++++++++
	
	
