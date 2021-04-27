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
from Compute_Extremes import HW_numdays
from Compute_Extremes import GT_numdays
from Compute_Extremes import CDGT_numdays

Cenario='rcp85ec'

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



for cmod in range(0,nMod):

        url=driveSAVE+'/historical/tasmax_P90/PI_tasmaxP90_EUR-11_'+GCMs[cmod]+'_historical_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(1971)+str(2000)+'.nc' 
        ncfile=Dataset(url,'r') 
	tasmax_p90=np.array(ncfile.variables['tasmax_p90']) 
        ncfile.close() 

      	tasmax,tasmax_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmax',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
	tasmax=tasmax-273.15

	if tasmax_flag==1:
		Lt,Ly,Lx=tasmax.shape
	

		#++++++++++++++++++++++++++++++++++++++++++++++++
		# MAXIMUM NUMBER OF CONSECUTIVE DAYS WITH Tx>35
                CDTxg35_MON,CDTxg35_SAZ,CDTxg35_CLI=CDGT_numdays(tasmax,35,ndias[cmod],Ano0i,Anoi,Anof)

                # save TMAX simple average
                dirsave=driveSAVE+'/'+Cenario+'/CDTxg35'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_CDTxg35_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('CDTxg35_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=CDTxg35_MON[:,:,:]
                TMmonw.units='days'
                TMmonw.long_name='Climatological Maximum number of consecutive days with Tmax>35 degC for each month'

                TMsazw = f.createVariable('CDTxg35_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=CDTxg35_SAZ[:,:,:]
                TMsazw.units='days'
                TMsazw.long_name='Climatological Maximum number of consecutive days with Tmax>35 degC for each season'

                TMcliw = f.createVariable('CDTxg35_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=CDTxg35_CLI[:,:]
                TMcliw.units='days'
                TMcliw.long_name='Climatological Maximum number of consecutive days with Tmax>35 degC'

                f.close()
                #+++++++++++++++++++




