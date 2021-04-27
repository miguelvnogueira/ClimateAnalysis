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
from Compute_Extremes import Simple_Sum
from Compute_Extremes import GT_numdays
from Compute_Extremes import CDGT_numdays
from Compute_Extremes import CDLT_numdays
from Compute_Extremes import MAXac_5days
from Compute_Extremes import GT_Sum
from Compute_Extremes import AVCDGT_numdays

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
#for cmod in [28]:
	print(GCMs[cmod]+'/'+RCMs[cmod])
	
      	pr,pr_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'pr',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
	pr=pr*60*60*24


	if pr_flag==1:
		Lt,Ly,Lx=pr.shape
		pr[pr>1000]=np.NaN
	
		#++++++++++++++++++++++++++++
		# Pr simple avearge
		pr_MON,pr_SAZ,pr_CLI=Simple_Average(pr,ndias[cmod],Ano0i,Anoi,Anof)	

                # save TMAX simple average
                dirsave=driveSAVE+'/'+Cenario+'/Prd'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_Prd_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('Prd_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=pr_MON[:,:,:]
                TMmonw.units='mm/day'
                TMmonw.long_name='Climatological monthly averaged daily mean precipitation'

                TMsazw = f.createVariable('Prd_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=pr_SAZ[:,:,:]
                TMsazw.units='mm/day'
                TMsazw.long_name='Climatological seasonal averaged daily mean precipitation'

                TMcliw = f.createVariable('Prd_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=pr_CLI[:,:]
                TMcliw.units='mm/day'
                TMcliw.long_name='Climatological yearly averaged daily mean precipitation'

                f.close()
		#+++++++++++++++++++
	

                #++++++++++++++++++++++++++++
                # Pr simple sum
                prac_MON,prac_SAZ,prac_CLI=Simple_Sum(pr,ndias[cmod],Ano0i,Anoi,Anof)    

                # save pr simple sum
                dirsave=driveSAVE+'/'+Cenario+'/Prac'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_Prac_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('Prac_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=prac_MON[:,:,:]
                TMmonw.units='mm'
                TMmonw.long_name='Climatological monthly accumulated precipitation per year'

                TMsazw = f.createVariable('Prac_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=prac_SAZ[:,:,:]
                TMsazw.units='mm'
                TMsazw.long_name='Climatological seasonal accumulated precipitation per year'

                TMcliw = f.createVariable('Prac_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=prac_CLI[:,:]
                TMcliw.units='mm'
                TMcliw.long_name='Climatological average yearly accumulated precipitation'

                f.close()
                #+++++++++++++++++++

                #+++++++++++++++++++++++++++
                # Max Pr accumulation in 5 days
                prac_MON,prac_SAZ,prac_CLI=MAXac_5days(pr,ndias[cmod],Ano0i,Anoi,Anof)

                # save pr simple sum
                dirsave=driveSAVE+'/'+Cenario+'/MaxPrac5d'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_MaxPrac5d_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('MaxPrac5d_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=prac_MON[:,:,:]
                TMmonw.units='mm'
                TMmonw.long_name='Maximum precipitation accumultation during 5 days by month'

                TMsazw = f.createVariable('MaxPrac5d_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=prac_SAZ[:,:,:]
                TMsazw.units='mm'
                TMsazw.long_name='Maximum precipitation accumultation during 5 days by season'

                TMcliw = f.createVariable('MaxPrac5d_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=prac_CLI[:,:]
                TMcliw.units='mm'
                TMcliw.long_name='Maximum precipitation accumultation during 5 days by month'

                f.close()
                #+++++++++++++++++++



                #+++++++++++++++++++++++++++++++
                # PERCENTAGE OF RAINFALL COMING FROM DAYS WITH pr>th
                for th in [10,50]:
                        PPrg_MON,PPrg_SAZ,PPrg_CLI=GT_Sum(pr,th,ndias[cmod],Ano0i,Anoi,Anof)
 
                        # save % precipitacao em  Dias Excedendo threshold
                        dirsave=driveSAVE+'/'+Cenario+'/PPrg'+str(th)
                        if os.path.isdir(dirsave)==False:
                                os.mkdir(dirsave)

                        urlsave=dirsave+'/PI_PPrg'+str(th)+'_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                        f = Dataset(urlsave,'w',format='NETCDF4')
                        f.createDimension('mm',12)
                        f.createDimension('ss',4)
                        f.createDimension('yy',Ly)
                        f.createDimension('xx',Lx)

                        TMmonw = f.createVariable('PPrg'+str(th)+'_mon', 'd', ('mm','yy','xx'))
                        TMmonw[:,:,:]=PPrg_MON[:,:,:]
                        TMmonw.units='%'
                        TMmonw.long_name='Percentage of total precipitation from days with Precipitation>'+str(th)+' mm/day by month'

                        TMsazw = f.createVariable('PPrg'+str(th)+'_saz', 'd', ('ss','yy','xx'))
                        TMsazw[:,:,:]=PPrg_SAZ[:,:,:]
                        TMsazw.units='%'
                        TMsazw.long_name='Percentage of total precipitation from days with Precipitation>'+str(th)+' mm/day by season'

                        TMcliw = f.createVariable('PPrg'+str(th)+'_cli', 'd', ('yy','xx'))
                        TMcliw[:,:]=PPrg_CLI[:,:]
                        TMcliw.units='%'
                        TMcliw.long_name='Percentage of total precipitation from days with Precipitation>'+str(th)+' mm/day'

                        f.close()
                #+++++++++++++++++++++++++++++++++++++



		#+++++++++++++++++++++++++++++++
		# NUMBER OF DAYS WITH pr>th
		for th in [1,10,20,50]:
			Prg_MON,Prg_SAZ,Prg_CLI=GT_numdays(pr,th,ndias[cmod],Ano0i,Anoi,Anof)

	                # save Num Dias Excedendo threshold
                	dirsave=driveSAVE+'/'+Cenario+'/Prg'+str(th)
                	if os.path.isdir(dirsave)==False:
                        	os.mkdir(dirsave)

	                urlsave=dirsave+'/PI_Prg'+str(th)+'_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
        	        f = Dataset(urlsave,'w',format='NETCDF4')
                	f.createDimension('mm',12)
	                f.createDimension('ss',4)
        	        f.createDimension('yy',Ly)
                	f.createDimension('xx',Lx)

	                TMmonw = f.createVariable('Prg'+str(th)+'_mon', 'd', ('mm','yy','xx'))
        	        TMmonw[:,:,:]=Prg_MON[:,:,:]
                	TMmonw.units='number of days'
	                TMmonw.long_name='Average number of days per year with Precipitation>'+str(th)+' mm/day by month'

        	        TMsazw = f.createVariable('Prg'+str(th)+'_saz', 'd', ('ss','yy','xx'))
                	TMsazw[:,:,:]=Prg_SAZ[:,:,:]
	                TMsazw.units='number of days'
        	        TMsazw.long_name='Average number of days per year with Precipitation>'+str(th)+' mm/day by season'

	                TMcliw = f.createVariable('Prg'+str(th)+'_cli', 'd', ('yy','xx'))
        	        TMcliw[:,:]=Prg_CLI[:,:]
	                TMcliw.units='number of days'
        	        TMcliw.long_name='Average number of days per year with Precipitation>'+str(th)+' mm/day'

                	f.close()
		#+++++++++++++++++++++++++++++++++++++


		#++++++++++++++++++++++++++++++++++++++++++++++++
		# MAXIMUM NUMBER OF CONSECUTIVE DAYS WITH pr>1
                CDPrg1_MON,CDPrg1_SAZ,CDPrg1_CLI=CDGT_numdays(pr,1,ndias[cmod],Ano0i,Anoi,Anof)

                # save CDPrg1
                dirsave=driveSAVE+'/'+Cenario+'/CDPrg1'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_CDPrg1_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('CDPrg1_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=CDPrg1_MON[:,:,:]
                TMmonw.units='days'
                TMmonw.long_name='Climatological Maximum number of consecutive days with precipitaiton>1 mm/day by month'

                TMsazw = f.createVariable('CDPrg1_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=CDPrg1_SAZ[:,:,:]
                TMsazw.units='days'
                TMsazw.long_name='Climatological Maximum number of consecutive days with precipitaiton>1 mm/day by season'

                TMcliw = f.createVariable('CDPrg1_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=CDPrg1_CLI[:,:]
                TMcliw.units='days'
                TMcliw.long_name='Climatological Maximum number of consecutive days with precipitaiton>1 mm/day'

                f.close()
                #+++++++++++++++++++


                #++++++++++++++++++++++++++++++++++++++++++++++++
                # MAXIMUM NUMBER OF CONSECUTIVE DAYS WITH pr<1
                CDPrl1_MON,CDPrl1_SAZ,CDPrl1_CLI=CDLT_numdays(pr,1,ndias[cmod],Ano0i,Anoi,Anof)

                # save CDPrl1
                dirsave=driveSAVE+'/'+Cenario+'/CDPrl1'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_CDPrl1_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('CDPrl1_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=CDPrl1_MON[:,:,:]
                TMmonw.units='days'
                TMmonw.long_name='Climatological Maximum number of consecutive days with precipitaiton<1 mm/day by month'

                TMsazw = f.createVariable('CDPrl1_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=CDPrl1_SAZ[:,:,:]
                TMsazw.units='days'
                TMsazw.long_name='Climatological Maximum number of consecutive days with precipitaiton<1 mm/day by season'

                TMcliw = f.createVariable('CDPrl1_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=CDPrl1_CLI[:,:]
                TMcliw.units='days'
                TMcliw.long_name='Climatological Maximum number of consecutive days with precipitaiton<1 mm/day'

                f.close()
                #+++++++++++++++++++



                #++++++++++++++++++++++++++++++++++++++++++++++++
                # AVERAGE NUMBER OF CONSECUTIVE DAYS WITH pr>1
                AVCDPrg1_MON,AVCDPrg1_SAZ,AVCDPrg1_CLI=AVCDGT_numdays(pr,1,ndias[cmod],Ano0i,Anoi,Anof)

                # save AVCDPrg1
                dirsave=driveSAVE+'/'+Cenario+'/AVCDPrg1'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_AVCDPrg1_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('AVCDPrg1_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=AVCDPrg1_MON[:,:,:]
                TMmonw.units='days'
                TMmonw.long_name='Climatological Average number of consecutive days with precipitaiton>1 mm/day by month'

                TMsazw = f.createVariable('AVCDPrg1_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=AVCDPrg1_SAZ[:,:,:]
                TMsazw.units='days'
                TMsazw.long_name='Climatological Average number of consecutive days with precipitaiton>1 mm/day by season'

                TMcliw = f.createVariable('AVCDPrg1_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=AVCDPrg1_CLI[:,:]
                TMcliw.units='days'
                TMcliw.long_name='Climatological Average number of consecutive days with precipitaiton>1 mm/day'

                f.close()
                #+++++++++++++++++++

		
