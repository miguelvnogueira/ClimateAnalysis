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

#for cmod in range(0,nMod):
for cmod in [17]:

      	vh10,vh10_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'sfcWind',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
	

	if vh10_flag==1:
		Lt,Ly,Lx=vh10.shape
	
		#++++++++++++++++++++++++++++
		# vh10 simple avearge
		vh10_MON,vh10_SAZ,vh10_CLI=Simple_Average(vh10,ndias[cmod],Ano0i,Anoi,Anof)	

                # save Vh10 simple average
                dirsave=driveSAVE+'/'+Cenario+'/Vh10'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_Vh10_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('Vh10_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=vh10_MON[:,:,:]
                TMmonw.units='m/s'
                TMmonw.long_name='Climatological monthly averaged daily average 10-m wind speed'

                TMsazw = f.createVariable('Vh10_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=vh10_SAZ[:,:,:]
                TMsazw.units='m/s'
                TMsazw.long_name='Climatological seasonal averaged daily average 10-m wind speed'

                TMcliw = f.createVariable('Vh10_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=vh10_CLI[:,:]
                TMcliw.units='m/s'
                TMcliw.long_name='Climatological yearly averaged daily average 10-m wind speed'

                f.close()
		#+++++++++++++++++++
	

                #++++++++++++++++++++++++++++
                # vh10 simple maximum
                Maxvh10_MON,Maxvh10_SAZ,Maxvh10_CLI=Simple_Maximum(vh10,ndias[cmod],Ano0i,Anoi,Anof)

                # save Vh10 simple Maximum
                dirsave=driveSAVE+'/'+Cenario+'/MaxVh10'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_MaxVh10_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('MaxVh10_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=Maxvh10_MON[:,:,:]
                TMmonw.units='m/s'
                TMmonw.long_name='Climatological maximum  daily average 10-m wind speed per month'

                TMsazw = f.createVariable('MaxVh10_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=Maxvh10_SAZ[:,:,:]
                TMsazw.units='m/s'
                TMsazw.long_name='Climatological maximum daily average 10-m wind speed per season'

                TMcliw = f.createVariable('MaxVh10_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=Maxvh10_CLI[:,:]
                TMcliw.units='m/s'
                TMcliw.long_name='Climatological maximum daily average 10-m wind speed'

                f.close()
                #+++++++++++++++++++

		#+++++++++++++++++++++++++++++++
		# NUMBER OF DAYS WITH Vh10>th
		for th in [5.5,10.8,20.8,33.3]:
			vh10g_MON,vh10g_SAZ,vh10g_CLI=GT_numdays(vh10,th,ndias[cmod],Ano0i,Anoi,Anof)

	                # save Num Dias Excedendo threshold
                	dirsave=driveSAVE+'/'+Cenario+'/Vh10g'+str(th)
                	if os.path.isdir(dirsave)==False:
                        	os.mkdir(dirsave)

	                urlsave=dirsave+'/PI_Vh10g'+str(th)+'_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
        	        f = Dataset(urlsave,'w',format='NETCDF4')
                	f.createDimension('mm',12)
	                f.createDimension('ss',4)
        	        f.createDimension('yy',Ly)
                	f.createDimension('xx',Lx)

	                TMmonw = f.createVariable('Vh10g'+str(th)+'_mon', 'd', ('mm','yy','xx'))
        	        TMmonw[:,:,:]=vh10g_MON[:,:,:]
                	TMmonw.units='number of days'
	                TMmonw.long_name='Average number of days per year with Vh10>'+str(th)+' m/s by month'

        	        TMsazw = f.createVariable('Vh10g'+str(th)+'_saz', 'd', ('ss','yy','xx'))
                	TMsazw[:,:,:]=vh10g_SAZ[:,:,:]
	                TMsazw.units='number of days'
        	        TMsazw.long_name='Average number of days per year with Vh10>'+str(th)+' m/s by season'

	                TMcliw = f.createVariable('Vh10g'+str(th)+'_cli', 'd', ('yy','xx'))
        	        TMcliw[:,:]=vh10g_CLI[:,:]
	                TMcliw.units='number of days'
        	        TMcliw.long_name='Average number of days per year with Vh10>'+str(th)+' m/s'

                	f.close()
		#+++++++++++++++++++++++++++++++++++++


                #+++++++++++++++++++++++++++++++
                # NUMBER OF DAYS WITH Vh10<th
                for th in [2]:
                        vh10l_MON,vh10l_SAZ,vh10l_CLI=LT_numdays(vh10,th,ndias[cmod],Ano0i,Anoi,Anof)

                        # save Num Dias abaixo threshold
                        dirsave=driveSAVE+'/'+Cenario+'/Vh10l'+str(th)
                        if os.path.isdir(dirsave)==False:
                                os.mkdir(dirsave)

                        urlsave=dirsave+'/PI_Vh10l'+str(th)+'_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                        f = Dataset(urlsave,'w',format='NETCDF4')
                        f.createDimension('mm',12)
                        f.createDimension('ss',4)
                        f.createDimension('yy',Ly)
                        f.createDimension('xx',Lx)

                        TMmonw = f.createVariable('Vh10l'+str(th)+'_mon', 'd', ('mm','yy','xx'))
                        TMmonw[:,:,:]=vh10l_MON[:,:,:]
                        TMmonw.units='number of days'
                        TMmonw.long_name='Average number of days per year with Vh10<'+str(th)+' m/s by month'

                        TMsazw = f.createVariable('Vh10l'+str(th)+'_saz', 'd', ('ss','yy','xx'))
                        TMsazw[:,:,:]=vh10l_SAZ[:,:,:]
                        TMsazw.units='number of days'
                        TMsazw.long_name='Average number of days per year with Vh10<'+str(th)+' m/s by season'

                        TMcliw = f.createVariable('Vh10l'+str(th)+'_cli', 'd', ('yy','xx'))
                        TMcliw[:,:]=vh10l_CLI[:,:]
                        TMcliw.units='number of days'
                        TMcliw.long_name='Average number of days per year with Vh10<'+str(th)+' m/s'

                        f.close()
                #+++++++++++++++++++++++++++++++++++++


		
