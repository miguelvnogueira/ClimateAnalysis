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
from Compute_UTCI import Compute_Tmrt
from Compute_UTCI import Compute_UTCI_approx
from Compute_Variables import Compute_RH

from Compute_Extremes import Simple_Average
from Compute_UTCI import Class_numdays

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



#for cmod in range(0,nMod):
for cmod in [1, 7, 11, 12, 17, 25, 26, 28, 33, 35, 37, 38, 39]:

	flags=np.zeros([9])
        tasmax,flags[0]=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmax',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        tasmin,flags[1]=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmin',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        ps,flags[2]=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'ps',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        huss,flags[3]=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'huss',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        rlus,flags[4]=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'rlus',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        rlds,flags[5]=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'rlds',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        rsus,flags[6]=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'rsus',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        rsds,flags[7]=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'rsds',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
 	vh10,flags[8]=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'sfcWind',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)


        if np.min(flags)==1:
                tasmea = (tasmax+tasmin)/2.
                del tasmax
                del tasmin
                RH,RH_flag=Compute_RH(tasmea,ps,huss)	
	        del ps
                del huss
		tasmea=tasmea-273.15
		Tmrt=Compute_Tmrt(rlus,rlds,rsus,rsds)
		Tmrt=Tmrt-273.15
		UTCI=Compute_UTCI_approx(tasmea,RH,Tmrt,vh10)
 		UTCI_flag=1
		del Tmrt
		del RH
		del vh10
		del tasmea
        else:
                UTCI_flag=0


	if UTCI_flag==1:
		Lt,Ly,Lx=UTCI.shape
	
		#++++++++++++++++++++++++++++
		# Tmax simple avearge
		UTCI_MON,UTCI_SAZ,UTCI_CLI=Simple_Average(UTCI,ndias[cmod],Ano0i,Anoi,Anof)	
	
		# save TMAX simple average
                dirsave=driveSAVE+'/'+Cenario+'/UTCI'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_UTCI_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('UTCI_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=UTCI_MON[:,:,:]
                TMmonw.units='deg C'
                TMmonw.long_name='Climatological monthly averaged UTCI'

                TMsazw = f.createVariable('UTCI_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=UTCI_SAZ[:,:,:]
                TMsazw.units='deg C'
                TMsazw.long_name='Climatological seasonal averaged UTCI'

                TMcliw = f.createVariable('UTCI_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=UTCI_CLI[:,:]
                TMcliw.units='deg C'
                TMcliw.long_name='Climatological yearly averaged UTCI'

                f.close()
		#+++++++++++++++++++
	

		#+++++++++++++++++++++++++++++++
		# NUMBER OF DAYS in each class of UTCI
		for th in ['ehs','vhs','shs','mhs','nhs','ncs','mcs','scs','vcs','ecs']:
			if th=='ehs':
				Threshu=9999.0
				Threshd=46.0
				utcileg='extreme heat stress'
                        elif th=='vhs':
                                Threshu=46.0
                                Threshd=38.0
                                utcileg='very strong heat stress'
                        elif th=='shs':
                                Threshu=38.0
                                Threshd=32.0
                                utcileg='strong heat stress'
                        elif th=='mhs':
                                Threshu=32.0
                                Threshd=26.0
                                utcileg='moderate heat stress'
                        elif th=='nhs':
                                Threshu=26.0
                                Threshd=9.0
                                utcileg='no heat stress'
                        elif th=='ncs':
                                Threshu=9.0
                                Threshd=0.0
                                utcileg='slight cold stres'
                        elif th=='mcs':
                                Threshu=0.0
                                Threshd=-13.0
                                utcileg='moderate cold stres'
                        elif th=='scs':
                                Threshu=-13.0
                                Threshd=-27.0
                                utcileg='strong cold stres'
                        elif th=='vcs':
                                Threshu=-27.0
                                Threshd=-40.0
                                utcileg='very strong cold stres'
                        elif th=='ecs':
                                Threshu=-40.0
                                Threshd=-9999.0
                                utcileg='extreme cold stres'





			Cl_MON,Cl_SAZ,Cl_CLI=Class_numdays(UTCI,Threshd,Threshu,ndias[cmod],Ano0i,Anoi,Anof)

	                # save Num Dias Excedendo threshold
                	dirsave=driveSAVE+'/'+Cenario+'/UTCI'+str(th)
                	if os.path.isdir(dirsave)==False:
                        	os.mkdir(dirsave)

	                urlsave=dirsave+'/PI_UTCI'+str(th)+'_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
        	        f = Dataset(urlsave,'w',format='NETCDF4')
                	f.createDimension('mm',12)
	                f.createDimension('ss',4)
        	        f.createDimension('yy',Ly)
                	f.createDimension('xx',Lx)

	                TMmonw = f.createVariable('UTCI'+str(th)+'_mon', 'd', ('mm','yy','xx'))
        	        TMmonw[:,:,:]=Cl_MON[:,:,:]
                	TMmonw.units='number of days'
	                TMmonw.long_name='Average number of days per year with '+utcileg+' for each month'

        	        TMsazw = f.createVariable('UTCI'+str(th)+'_saz', 'd', ('ss','yy','xx'))
                	TMsazw[:,:,:]=Cl_SAZ[:,:,:]
	                TMsazw.units='number of days'
        	        TMsazw.long_name='Average number of days per year with '+utcileg+' for each season'

	                TMcliw = f.createVariable('UTCI'+str(th)+'_cli', 'd', ('yy','xx'))
        	        TMcliw[:,:]=Cl_CLI[:,:]
	                TMcliw.units='number of days'
        	        TMcliw.long_name='Average number of days per year with '+utcileg+' degC for the entire period'

                	f.close()
		#+++++++++++++++++++++++++++++++++++++



