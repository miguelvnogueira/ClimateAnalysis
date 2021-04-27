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


Cenario='rcp85ec'

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

        tasmax,tasmax_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmax',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        tasmax=tasmax-273.15

        tasmin,tasmin_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmin',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
        tasmin=tasmin-273.15


	if tasmax_flag==1 and tasmin_flag==1:
	        tasmea = (tasmax+tasmin)/2.
		del tasmax
		del tasmin

		Lt,Ly,Lx=tasmea.shape

		TM=np.zeros([30*12,Ly,Lx])
                NM=np.zeros([30*12,Ly,Lx])

		# Temp simple avearge
	        Ano=Ano0i
        	Mes=1
        	Dia=1
		cMes=0
        	for t in range(0,Lt):
                	if Ano>=Anoi and Ano<=Anof:
                  	     	TM[cMes,:,:]=TM[cMes,:,:]+tasmea[t,:,:]
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
	

		TM=TM/NM
		TM[NM==0]=np.NaN

		# save TMEA simple average
                dirsave=driveSAVE+'/Monthly/TM/'+Cenario
                if os.path.isdir(dirsave)==False:
                        os.mkdir(dirsave)

                urlsave=dirsave+'/PI_TM_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
                f = Dataset(urlsave,'w',format='NETCDF4')
                f.createDimension('mm',12*30)
                f.createDimension('yy',Ly)
                f.createDimension('xx',Lx)

                TMmonw = f.createVariable('TM', 'd', ('mm','yy','xx'))
                TMmonw[:,:,:]=TM[:,:,:]
                TMmonw.units='deg C'
                TMmonw.long_name='Monthly averaged daily mean temperature'
                f.close()
		#+++++++++++++++++++
	

		
