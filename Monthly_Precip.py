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


Cenario='rcp85bc'

driveCORDEX,driveSAVE,Ano0i,Ano0f,Anoi,Anof,Cenarioleg=Cenario_Props(Cenario)

GCMs=np.load(driveSAVE+'/Lista_GCMs.npy')
RCMs=np.load(driveSAVE+'/Lista_RCMs.npy')
EXPs=np.load(driveSAVE+'/Lista_EXPs.npy')
VERs=np.load(driveSAVE+'/Lista_VERs.npy')
ndias=np.load(driveSAVE+'/Lista_ndias.npy')
nMod=len(GCMs)


                
# CRIAR PASTA DO CENARIO SE NAO EXISTIR                
dirsave=driveSAVE+'/Monthly'
if os.path.isdir(dirsave)==False:                        
	os.mkdir(dirsave)


print(Cenario)

for cmod in range(0,nMod):
#for cmod in [35]:


      	pr,pr_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'pr',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
	pr=pr*60*60*24

	if pr_flag==1:
		Lt,Ly,Lx=pr.shape
		PrM=np.zeros([30*12,Ly,Lx])
                NM=np.zeros([30*12,Ly,Lx])

		# Pr simple avearge
	        Ano=Ano0i
        	Mes=1
        	Dia=1
		cMes=0
        	for t in range(0,Lt):
                	if Ano>=Anoi and Ano<=Anof:
                  	     	PrM[cMes,:,:]=PrM[cMes,:,:]+pr[t,:,:]
				NM[cMes,:,:]=NM[cMes,:,:]+1
	
				
              		Dia=Dia+1
                	if ndias[cmod]==360:
                        	if Dia>30:
                                	Dia=1
                                	Mes=Mes+1
					if Ano>=Anoi and Ano<=Anof:
						cMes=cMes+1
                	elif ndias[cmod]==365: # NO LEAP YEARS
                        	if Dia > calendar.monthrange(2001,Mes)[1]:
                                	Dia=1
                                	Mes=Mes+1
					if Ano>=Anoi and Ano<=Anof:
						cMes=cMes+1
                	elif ndias[cmod]==366: # WITH LEAP YEARS
                        	if Dia > calendar.monthrange(Ano,Mes)[1]:
                                	Dia=1
                                	Mes=Mes+1
					if Ano>=Anoi and Ano<=Anof:
						cMes=cMes+1
                	if Mes>12:
                        	Mes=1
                        	Ano=Ano+1
	

		PrM=PrM/NM
		PrM[NM==0]=np.NaN

                # save TMAX simple average
                dirsave=driveSAVE+'/Monthly/PrM/'+Cenario
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_PrM_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12*30)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('PrM', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=PrM[:,:,:]
                TMmonw.units='mm'
                TMmonw.long_name='Monthly averaged precipitation'
                f.close()
		#+++++++++++++++++++
	

		
