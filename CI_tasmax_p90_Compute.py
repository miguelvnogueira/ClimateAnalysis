import numpy as np
from netCDF4 import Dataset
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import calendar
import os

Cenario='historical'

sys.path.append("../Funcoes")

from Load_Models import Cenario_Props
from Load_Models import Ler_Modelo

driveCORDEX,driveSAVE,Ano0i,Ano0f,Anoi,Anof,Cenarioleg=Cenario_Props(Cenario)


GCMs=np.load(driveSAVE+'/Lista_GCMs.npy')
RCMs=np.load(driveSAVE+'/Lista_RCMs.npy')
EXPs=np.load(driveSAVE+'/Lista_EXPs.npy')
VERs=np.load(driveSAVE+'/Lista_VERs.npy')
ndias=np.load(driveSAVE+'/Lista_ndias.npy')
nMod=len(GCMs)

for cmod in range(0,nMod):
    print(GCMs[cmod]+'-'+RCMs[cmod])

    tasmax,tasmax_flag=Ler_Modelo(driveCORDEX,cmod,Ano0i,Ano0f,'tasmax',GCMs[cmod],RCMs[cmod],EXPs[cmod],VERs[cmod],Cenarioleg)
    tasmax=tasmax-273.15
    if tasmax_flag==1:			
		
	Lt,Ly,Lx=tasmax.shape

	if ndias[cmod]>=365:	
		ndmax=365
	else:
		ndmax=360
		
	data=np.zeros([ndmax,30*31,Ly,Lx])*np.NaN
	
	Ano=Ano0i
	Mes=1
	Dia=1
	Diaj=1
	for t in range(0,Lt):
		if Ano>=Anoi and Ano<=Anof:
			if ndias[cmod]==366 and Mes==2 and Dia==29:
				'skip'
			else:
				ti=np.max([t-15,0])
				tf=np.min([t+16,Lt])
			
				Lta=len(tasmax[ti:tf,0,0])
				data[Diaj-1,(Ano-Anoi)*31:(Ano-Anoi)*31+Lta,:,:]=tasmax[ti:tf,:,:]
				
		if ndias[cmod]==366 and Mes==2 and Dia==29:
			'skip'
		else:
			Diaj=Diaj+1		

		Dia=Dia+1
		if ndias[cmod]==360:
			if Dia>30:
				Dia=1
				Mes=Mes+1
		elif ndias[cmod]==365: # NO LEAP YEARS
			if Dia > calendar.monthrange(2001,Mes)[1]:
				Dia=1
				Mes=Mes+1
                elif ndias[cmod]==366: # WITH LEAP YEARS
                        if Dia > calendar.monthrange(Ano,Mes)[1]:
				Dia=1
				Mes=Mes+1
		else:
			print(ndias[cmod])
			sys.exit()
		

		if Mes>12:
			Mes=1
			Ano=Ano+1
			Diaj=1


	
	P90=np.zeros([ndmax,Ly,Lx])

	for diaj in range(0,ndmax):
		P90[diaj,:,:]=np.nanpercentile(data[diaj,:,:,:],90,axis=0)
		




	
	# CRIAR PASTA DO CENARIO SE NAO EXISTIR
	dirsave=driveSAVE+'/'+Cenario
	if os.path.isdir(dirsave)==False:
	        os.mkdir(dirsave)



	# save 90th percentiles
	dirsave=driveSAVE+'/'+Cenario+'/tasmax_P90'
        if os.path.isdir(dirsave)==False:
                os.mkdir(dirsave)

	urlsave=dirsave+'/PI_tasmaxP90_EUR-11_'+GCMs[cmod]+'_'+Cenario+'_'+EXPs[cmod]+'_'+RCMs[cmod]+'_'+VERs[cmod]+'_day_'+str(Anoi)+str(Anof)+'.nc'
        f = Dataset(urlsave,'w',format='NETCDF4')
        f.createDimension('mm',ndmax)
        f.createDimension('yy',Ly)
        f.createDimension('xx',Lx)

       	TMmonw = f.createVariable('tasmax_p90', 'd', ('mm','yy','xx'))
        TMmonw[:,:,:]=P90[:,:,:]
        TMmonw.units='deg C'

	f.close()

	del P90
	del TMmonw
	del data
   
