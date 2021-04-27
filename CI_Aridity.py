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

Cenario='historical'

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

      	pet,pet_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'pet',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        pr,pr_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'pr',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        pr=pr*60*60*24
	if pet_flag==1 and pr_flag==1:
	
		print(np.mean(pr))
		print(np.mean(pet))
		AI=pr/pet

		print(np.mean(AI))
		sys.exit()
		del pr
		del pet

		Lt,Ly,Lx=AI.shape
	
		#++++++++++++++++++++++++++++
		# AI simple avearge
		AI_MON,AI_SAZ,AI_CLI=Simple_Average(AI,ndias[cmod],Ano0i,Anoi,Anof)	

                # save TMAX simple average
                dirsave=driveSAVE+'/'+Cenario+'/AI'
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_AI_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12)
                f.createDimension('ss',4)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('AI_mon', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=AI_MON[:,:,:]
                TMmonw.units='AI'
                TMmonw.long_name='Climatological averaged Aridity Index by month'

                TMsazw = f.createVariable('AI_saz', 'd', ('ss','yy','xx'))
                TMsazw[:,:,:]=AI_SAZ[:,:,:]
                TMsazw.units='adim'
                TMsazw.long_name='Climatological averaged Aridity Index by season'

                TMcliw = f.createVariable('AI_cli', 'd', ('yy','xx'))
                TMcliw[:,:]=AI_CLI[:,:]
                TMcliw.units='adim'
                TMcliw.long_name='Climatological Aridity Index'

                f.close()
		#+++++++++++++++++++
	
	
