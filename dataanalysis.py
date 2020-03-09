#Created on Wed Jan 22 14:15:11 2020
#@author: aarondomenzain

import numpy as np
import math
import matplotlib.pyplot as plt 


#Carga los archivos de datos
data=np.loadtxt("dataGAG.txt")
param=np.loadtxt("Parameters.txt")

if param.shape[0]==12:
	[dprot,Kprot,KRNA,Kprot,CgRNAin,CmRNA,GagIn,N,BdG,pfold01,pallost,dt]=param
if param.shape[0]==11:
	[dprot,Kprot,KRNA,Kprot,CgRNAin,CmRNA,GagIn,N,BdG,pfold01,pallost]=param

graf=data 


#Cálculo de promedios e incertidumbre

	#Número de puntos temporales
Nt=graf.shape[0]
	#Mitad entera de puntos demporales
N=int(Nt/2)


#Promedios tomados a partir de la mitad del tiempo total de simulación
avggagpop=np.mean(graf[int(Nt/2):,1])
avgcgrna=np.mean(graf[int(Nt/2):,2])
avgselect=np.mean(graf[int(0.6*Nt):,6])
avgcgags0=np.mean(graf[int(Nt/2):,7])
avgcgags1=avselect=np.mean(graf[int(Nt/2):,8])

#Desviación estándar 
stdgagpop=np.std(graf[int(Nt/2):,1])
stdcgrna=np.std(graf[int(Nt/2):,2])
stdselect=np.std(graf[int(Nt/2):,6])
stdcgags0=np.std(graf[int(Nt/2):,7])
stdcgags1=avselect=np.std(graf[int(Nt/2):,8])


#Construcción de arreglos para graficar
avgagpop=[avggagpop]*Nt
avcgrna=[avgcgrna]*Nt
avselect=[avgselect]*Nt
avcgags0=[avgcgags0]*Nt
avcgags1=[avgcgags1]*Nt


#state=[t*dt,gagpop,concgrna,vir,vlp,allostvir,selec,Cgags0,Cgags1]

#FIGURES
#plt.figure(0)
#plt.plot(graf[1:,1],graf[1:,4],label="Viriones")
#plt.legend()
#plt.title("Población de dímeros")

#plt.figure(0)
#plt.plot(graf[1:,1],graf[1:,5],label="VLP")
#plt.legend()


plt.figure(1)
plt.plot(graf[0:,0],graf[0:,1],label="Gag")
plt.plot(graf[0:,0],avgagpop,'k', dashes=[3,3])
plt.legend()
plt.title("Población de proteínas")
plt.savefig('gagpop.png')


plt.figure(2)
plt.plot(graf[0:,0],graf[0:,2],label="gRNA")
plt.plot(graf[0:,0],avcgrna, 'k', dashes=[3,3])
plt.legend()
plt.title("Concentración")

plt.figure(2)
plt.plot(graf[0:,0],graf[0:,6],label="Selectividad")
plt.plot(graf[0:,0],avselect, 'k', dashes=[3,3])
plt.legend()
plt.savefig('Concentracion.png')

plt.figure(4)
plt.plot(graf[0:,0],graf[0:,7],label="Gags en estado compacto")
plt.plot(graf[0:,0],avcgags0, 'k', dashes=[3,3])
plt.legend()

plt.figure(4)
plt.plot(graf[0:,0],graf[0:,8],label="Gags en estado extendido")
plt.plot(graf[0:,0],avcgags1, 'k', dashes=[3,3])
plt.legend()
plt.savefig('GagsComp-Ext.png')


plt.figure(3)
plt.plot(graf[0:,0],graf[0:,5],label="Viriones")
plt.legend()
plt.title("Población de viriones alostéricos")


#plt.title("Selectividad")

#plt.figure(3)
#plt.plot(graf[1:,1],graf[1:,7],label="VLP")
#plt.legend()
#plt.title("Población de virus alostéricos")



y='y'#input('Display parameters? Type yes/no or y/n:')
print()
if y=='y' or y=='yes':
	print()
	print('PARAMETERS')
	print()
	print('proteindiameter/length step:',dprot)
	print('KRNA:',KRNA/(dprot*80))
	print("Kprot",Kprot)
	print("Molecules of CmRNA:",CmRNA/(dprot*30))
	print("BdG:",BdG)
	print("Pfold01:",pfold01)
	print("Pallost:",pallost)
	print()

	print('AVERAGES,UNCERTAINTY')
	print('Average selectivity: ',round(np.mean(avselect),4),round(stdselect/(math.sqrt(N)),4))

	print('Average CgRNA: ',round(np.mean(avcgrna),4),round(stdcgrna/(math.sqrt(N)),4))
	print('Average Gag population: ', int(avgagpop[0]) ,round(stdgagpop/(math.sqrt(N)),4))
	print('Average concentration of compact Gag: ',round(avcgags0[0],4),round(stdcgags0/(math.sqrt(N)),4))
	print('Average concentration of extended Gag: ',round(avcgags1[0],4),round(stdcgags1/(math.sqrt(N)),4))
	print()

    
if param.shape[0]==12:
	print("Time of computing in s: ",dt)
	print()

x=input("Saved plots in png. Display? Type yes/no or y/n:")
if x=='yes' or x=='y':
	plt.show()














