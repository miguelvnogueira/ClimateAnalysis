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
#for cmod in [28]:
	print(GCMs[cmod]+'/'+RCMs[cmod])
	
      	pr,pr_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'pr',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
	pr=pr*60*60*24


	if pr_flag==1:
		Lt,Ly,Lx=pr.shape
		pr[pr>1000]=np.NaN
	
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
	
