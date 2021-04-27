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
from Compute_Extremes import Simple_Maximum
from Compute_Extremes import GT_numdays
from Compute_Variables import Compute_Vz

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


print(Cenario)

for cmod in range(0,nMod):


      	vh10,vh10_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'sfcWind',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)


	if vh10_flag==1:
		Lt,Ly,Lx=vh10.shape
	
		for z in [30.,60.]:
			Vz=Compute_Vz(vh10,z)

			#++++++++++++++++++++++++++++
			# Vz simple avearge
			Vz_MON,Vz_SAZ,Vz_CLI=Simple_Average(Vz,ndias[cmod],Ano0i,Anoi,Anof)	

	                # save RH simple average
	                dirsave=driveSAVE+'/'+Cenario+'/Vh'+str(z)
        	        if os.path.isdir(dirsave)==False:
                	        os.mkdir(dirsave)

	                urlsave=dirsave+'/PI_Vh'+str(z)+'_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
        	        f = Dataset(urlsave,'w',format='NETCDF4')
	                f.createDimension('mm',12)
        	        f.createDimension('ss',4)
                	f.createDimension('yy',Ly)
	                f.createDimension('xx',Lx)
	
        	        TMmonw = f.createVariable('Vh'+str(z)+'_mon', 'd', ('mm','yy','xx'))
                	TMmonw[:,:,:]=Vz_MON[:,:,:]
	                TMmonw.units='m/s'
        	        TMmonw.long_name='Climatological averaged daily average wind speed at '+str(z)+'m (power-law interpolation) by month'

                	TMsazw = f.createVariable('Vh'+str(z)+'_saz', 'd', ('ss','yy','xx'))
	                TMsazw[:,:,:]=Vz_SAZ[:,:,:]
        	        TMsazw.units='m/s'
                	TMsazw.long_name='Climatological averaged daily average wind speed at '+str(z)+'m (power-law interpolation) by season'

	                TMcliw = f.createVariable('Vh'+str(z)+'_cli', 'd', ('yy','xx'))
        	        TMcliw[:,:]=Vz_CLI[:,:]
                	TMcliw.units='m/s'
	                TMcliw.long_name='Climatological averaged daily average wind speed at '+str(z)+'m (power-law interpolation)'

        	        f.close()
			#+++++++++++++++++++
	
	
        	        #++++++++++++++++++++++++++++
	                # vhz simple maximum
                	MaxVz_MON,MaxVz_SAZ,MaxVz_CLI=Simple_Maximum(Vz,ndias[cmod],Ano0i,Anoi,Anof)


                        # save Vz max
                        dirsave=driveSAVE+'/'+Cenario+'/MaxVh'+str(z)
                        if os.path.isdir(dirsave)==False:
                                os.mkdir(dirsave)

                        urlsave=dirsave+'/PI_MaxVh'+str(z)+'_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                        f = Dataset(urlsave,'w',format='NETCDF4')
                        f.createDimension('mm',12)
                        f.createDimension('ss',4)
                        f.createDimension('yy',Ly)
                        f.createDimension('xx',Lx)

                        TMmonw = f.createVariable('MaxVh'+str(z)+'_mon', 'd', ('mm','yy','xx'))
                        TMmonw[:,:,:]=MaxVz_MON[:,:,:]
                        TMmonw.units='m/s'
                        TMmonw.long_name='Climatological maximum daily average wind speed at '+str(z)+'m (power-law interpolation) by month'

                        TMsazw = f.createVariable('MaxVh'+str(z)+'_saz', 'd', ('ss','yy','xx'))
                        TMsazw[:,:,:]=MaxVz_SAZ[:,:,:]
                        TMsazw.units='m/s'
                        TMsazw.long_name='Climatological maximum daily average wind speed at '+str(z)+'m (power-law interpolation) by season'

                        TMcliw = f.createVariable('MaxVh'+str(z)+'_cli', 'd', ('yy','xx'))
                        TMcliw[:,:]=MaxVz_CLI[:,:]
                        TMcliw.units='m/s'
                        TMcliw.long_name='Climatological maximum daily average wind speed at '+str(z)+'m (power-law interpolation)'

                        f.close()
                        #+++++++++++++++++++




	                #+++++++++++++++++++++++++++++++
        	        # NUMBER OF DAYS WITH Vz>th
                	for th in [33.3]:
                        	Vzg_MON,Vzg_SAZ,Vzg_CLI=GT_numdays(Vz,th,ndias[cmod],Ano0i,Anoi,Anof)


	                        # save Vzg max
	                        dirsave=driveSAVE+'/'+Cenario+'/Vh'+str(z)+'g'+str(th)
        	                if os.path.isdir(dirsave)==False:
                	                os.mkdir(dirsave)

	                        urlsave=dirsave+'/PI_Vh'+str(z)+'g'+str(th)+'_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
        	                f = Dataset(urlsave,'w',format='NETCDF4')
                	        f.createDimension('mm',12)
                        	f.createDimension('ss',4)
	                        f.createDimension('yy',Ly)
        	                f.createDimension('xx',Lx)

                	        TMmonw = f.createVariable('Vh'+str(z)+'g'+str(th)+'_mon', 'd', ('mm','yy','xx'))
                        	TMmonw[:,:,:]=Vzg_MON[:,:,:]
	                        TMmonw.units='number of days'
        	                TMmonw.long_name='Average number of days per year with Vh'+str(z)+'>'+str(th)+' m/s by season'

                	        TMsazw = f.createVariable('Vh'+str(z)+'g'+str(th)+'_saz', 'd', ('ss','yy','xx'))
                        	TMsazw[:,:,:]=Vzg_SAZ[:,:,:]
	                        TMsazw.units='number of days'
        	                TMsazw.long_name='Average number of days per year with Vh'+str(z)+'>'+str(th)+' m/s by season'

                	        TMcliw = f.createVariable('Vh'+str(z)+'g'+str(th)+'_cli', 'd', ('yy','xx'))
                        	TMcliw[:,:]=Vzg_CLI[:,:]
	                        TMcliw.units='number of days'
        	                TMcliw.long_name='Average number of days per year with Vh'+str(z)+'>'+str(th)+' m/s'

                	        f.close()
                        	#+++++++++++++++++++


			
