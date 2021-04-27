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

#for cmod in range(0,nMod):
for cmod in [12]:

      	gsr,gsr_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'rsds',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
	gsr=gsr*(60.*60.*24)/1e6

	if gsr_flag==1:
		Lt,Ly,Lx=gsr.shape
	
		#++++++++++++++++++++++++++++
		# gsr simple avearge
		gsr_MON,gsr_SAZ,gsr_CLI=Simple_Average(gsr,ndias[cmod],Ano0i,Anoi,Anof)	

                # save GSR simple average
                dirsave=driveSAVE+'/'+Cenario+'/gsr'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_gsr_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('gsr_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=gsr_MON[:,:,:]
                TMmonw.units='MJ/m**2'
                TMmonw.long_name='Climatological averaged daily global solar radiation per month'

                TMsazw = f.createVariable('gsr_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=gsr_SAZ[:,:,:]
                TMsazw.units='MJ/m**2'
                TMsazw.long_name='Climatological averaged daily global solar radiation per season'

                TMcliw = f.createVariable('gsr_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=gsr_CLI[:,:]
                TMcliw.units='MJ/m**2'
                TMcliw.long_name='Climatological averaged daily global solar radiation per month'

                f.close()
		#+++++++++++++++++++
	

