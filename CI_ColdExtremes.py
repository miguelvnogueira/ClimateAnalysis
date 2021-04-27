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
from Compute_Extremes import CW_numdays
from Compute_Extremes import GT_numdays
from Compute_Extremes import LT_numdays
from Compute_Extremes import CDLT_numdays

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
#for cmod in [12]:
	print(GCMs[cmod]+'/'+RCMs[cmod])


        url=driveSAVE+'/historical/tasmin_P10/PI_tasminP10_EUR-11_'+GCMs[cmod]+'_historical_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(1971)+str(2000)+'.nc' 
        ncfile=Dataset(url,'r') 
        tasmin_p10=np.array(ncfile.variables['tasmin_p10']) 
        ncfile.close() 
	
      	tasmin,tasmin_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmin',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
	tasmin=tasmin-273.15

	if tasmin_flag==1:
		Lt,Ly,Lx=tasmin.shape
	
		#++++++++++++++++++++++++++++
		# Tmin simple avearge
		Tn_MON,Tn_SAZ,Tn_CLI=Simple_Average(tasmin,ndias[cmod],Ano0i,Anoi,Anof)	

                # save TMAX simple average
                dirsave=driveSAVE+'/'+Cenario+'/Tn'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_Tn_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('Tn_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=Tn_MON[:,:,:]
                TMmonw.units='deg C'
                TMmonw.long_name='Climatological monthly averaged daily minimum 2m air temperature'

                TMsazw = f.createVariable('Tn_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=Tn_SAZ[:,:,:]
                TMsazw.units='deg C'
                TMsazw.long_name='Climatological seasonal averaged daily minimum 2m air temperature'

                TMcliw = f.createVariable('Tn_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=Tn_CLI[:,:]
                TMcliw.units='deg C'
                TMcliw.long_name='Climatological yearly averaged daily minimum 2m air temperature'

                f.close()
		#+++++++++++++++++++
	

		#+++++++++++++++++++++++++++++++
		# NUMBER OF DAYS WITH Tn>th
		for th in [20]:
			Tng_MON,Tng_SAZ,Tng_CLI=GT_numdays(tasmin,th,ndias[cmod],Ano0i,Anoi,Anof)

	                # save Num Dias Excedendo threshold
                	dirsave=driveSAVE+'/'+Cenario+'/Tng'+str(th)
                	if os.path.isdir(dirsave)==False:
                        	os.mkdir(dirsave)

	                urlsave=dirsave+'/PI_Tng'+str(th)+'_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
        	        f = Dataset(urlsave,'w',format='NETCDF4')
                	f.createDimension('mm',12)
	                f.createDimension('ss',4)
        	        f.createDimension('yy',Ly)
                	f.createDimension('xx',Lx)

	                TMmonw = f.createVariable('Tng'+str(th)+'_mon', 'd', ('mm','yy','xx'))
        	        TMmonw[:,:,:]=Tng_MON[:,:,:]
                	TMmonw.units='number of days'
	                TMmonw.long_name='Average number of days per year with Tmin>'+str(th)+' degC for each month'

        	        TMsazw = f.createVariable('Tng'+str(th)+'_saz', 'd', ('ss','yy','xx'))
                	TMsazw[:,:,:]=Tng_SAZ[:,:,:]
	                TMsazw.units='number of days'
        	        TMsazw.long_name='Average number of days per year with Tmin>'+str(th)+' degC for each season'

	                TMcliw = f.createVariable('Tng'+str(th)+'_cli', 'd', ('yy','xx'))
        	        TMcliw[:,:]=Tng_CLI[:,:]
	                TMcliw.units='number of days'
        	        TMcliw.long_name='Average number of days per year with Tmin>'+str(th)+' degC for the entire period'

                	f.close()
		#+++++++++++++++++++++++++++++++++++++

                #+++++++++++++++++++++++++++++++
                # NUMBER OF DAYS WITH Tn<th
                for th in [0]:
                        Tnl_MON,Tnl_SAZ,Tnl_CLI=LT_numdays(tasmin,th,ndias[cmod],Ano0i,Anoi,Anof)

                        # save Num Dias Excedendo threshold
                        dirsave=driveSAVE+'/'+Cenario+'/Tnl'+str(th)
                        if os.path.isdir(dirsave)==False:
                                os.mkdir(dirsave)

                        urlsave=dirsave+'/PI_Tnl'+str(th)+'_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                        f = Dataset(urlsave,'w',format='NETCDF4')
                        f.createDimension('mm',12)
                        f.createDimension('ss',4)
                        f.createDimension('yy',Ly)
                        f.createDimension('xx',Lx)

                        TMmonw = f.createVariable('Tnl'+str(th)+'_mon', 'd', ('mm','yy','xx'))
                        TMmonw[:,:,:]=Tnl_MON[:,:,:]
                        TMmonw.units='number of days'
                        TMmonw.long_name='Average number of days per year with Tmin<'+str(th)+' degC for each month'

                        TMsazw = f.createVariable('Tnl'+str(th)+'_saz', 'd', ('ss','yy','xx'))
                        TMsazw[:,:,:]=Tnl_SAZ[:,:,:]
                        TMsazw.units='number of days'
                        TMsazw.long_name='Average number of days per year with Tmin<'+str(th)+' degC for each season'

                        TMcliw = f.createVariable('Tnl'+str(th)+'_cli', 'd', ('yy','xx'))
                        TMcliw[:,:]=Tnl_CLI[:,:]
                        TMcliw.units='number of days'
                        TMcliw.long_name='Average number of days per year with Tmin<'+str(th)+' degC for the entire period'

                        f.close()
                #+++++++++++++++++++++++++++++++++++++


		#++++++++++++++++++++++++++++++++++++++++++++++++
		# MAXIMUM NUMBER OF CONSECUTIVE DAYS WITH Tn<7
                CDTnl7_MON,CDTnl7_SAZ,CDTnl7_CLI=CDLT_numdays(tasmin,7,ndias[cmod],Ano0i,Anoi,Anof)

                # save TMAX simple average
                dirsave=driveSAVE+'/'+Cenario+'/CDTnl7'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_CDTnl7_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('CDTnl7_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=CDTnl7_MON[:,:,:]
                TMmonw.units='days'
                TMmonw.long_name='Climatological Maximum number of consecutive days with Tmin<7 degC for each month'

                TMsazw = f.createVariable('CDTnl7_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=CDTnl7_SAZ[:,:,:]
                TMsazw.units='days'
                TMsazw.long_name='Climatological Maximum number of consecutive days with Tmin<7 degC for each season'

                TMcliw = f.createVariable('CDTnl7_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=CDTnl7_CLI[:,:]
                TMcliw.units='days'
                TMcliw.long_name='Climatological Maximum number of consecutive days with Tmin<7 degC'

                f.close()
                #+++++++++++++++++++




		# NUMBER OF DAYS in COLD WAVE
        	CWD_MON,CWD_SAZ,CWD_CLI=CW_numdays(tasmin,tasmin_p10,ndias[cmod],Ano0i,Anoi,Anof)

                # save 
                dirsave=driveSAVE+'/'+Cenario+'/CWD'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_CWD_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('CWD_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=CWD_MON[:,:,:]
                TMmonw.units='days'
                TMmonw.long_name='Average number of days in coldwave per year for each month'

                TMsazw = f.createVariable('CWD_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=CWD_SAZ[:,:,:]
                TMsazw.units='days'
                TMsazw.long_name='Average number of days in coldwave per year for each season'

                TMcliw = f.createVariable('CWD_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=CWD_CLI[:,:]
                TMcliw.units='days'
                TMcliw.long_name='Average number of days in coldwave per year'

                f.close()
                #+++++++++++++++++++
	
		
