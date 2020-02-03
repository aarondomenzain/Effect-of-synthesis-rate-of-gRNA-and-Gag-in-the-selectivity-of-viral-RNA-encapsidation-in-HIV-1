#Created on Wed Jan 22 14:15:11 2020
#@author: aarondomenzain

import numpy as np
import random
import math

#LISTA DE PARÁMETROS
Kprot=0.1 #tasa de producción de Gag
KRNA=0.01 #tasa de producción de gRNA
#tmax=100 #número de iteraciones en el tiempo
N=10000 #tamaño inicial de la recta numérica discreta
CgRNA=100 #cantidad de gRNA inicial
CmRNA=N-CgRNA #cantidad de mRNA inicial
GagIn=500 #cantidad de gag inicial
rprot=5 #radio promedio de proteína en unidades de pasos de random walk
vir=0 #cantidad de viriones
vlp=0 #cantidad de vlp
pfold01=0.1
BdG=11 #valor experimental BdG=11
Pfold=[[1-pfold01,pfold01],[pfold01*math.exp(-BdG),1-pfold01*math.exp(-BdG)]] 


#LOCALIZACION INICIAL DE GAG
#Argumentos: gagin=cantidad inicial de gag
#Salida: Matriz de estados Gag
def begin(gagin,cgrna,cmrna):
    #asigna posiciones y estados iniciales aleatorios
    n=cgrna+cmrna
    gag=np.zeros([gagin,2])
    for i in range (gagin):
        gag[i]=[random.randint(0,n),random.randint(0,1)] 
    #devuelve la matriz de posiciones y estados de gag
    return gag
    


#Random Walk
    #Hace que las posiciones de cada walker de un salto de tamaño d hacia la izq o der.
    #Argumentos: Gag, n=tamaño del dominio espacial, d=distancia del paso
    #Salida: Matriz de estados Gag
def randomwalk(gag,n,d):
    for elem in gag: 
        elem[0]=elem[0] + d*random.randrange(-1,2,2) # salta una distancia d hacia la derecha o a la izquierda
        #condiciones periódicas
        if elem[0]>n:
            elem[0]=elem[0]-n
        if elem[0]<0:
            elem[0]=elem[0]+n
    return gag




#Encapsidación
    #Compara distancias entre proteínas gag para determinar si sucede encapsidación
    #Si la encapsidación ocurre en mRNA, se formará un VLP
    #Si la encapsidación ocurre en gRNA, se formará un VIRIÓN
    #Argumentos: Gag ,d=radio efectivo de encapsidación , cgrna= cantidad de grna , vir=# de viriones , vlp=# de vlps
def encaps(gag,d,cmrna,cgrna,vir,vlp):
    par=[] #vector vacío que contendrá las etiquetas de las gags que formaron pareja
    n=cmrna+cgrna
    ngags=gag.shape[0]

    for i in range (ngags-1):
        for j in range (i+1,ngags):
            r=abs(gag[i][0]-gag[j][0]) #distancia entre gags
            m=min(gag[i][0],gag[j][0]) #posición de gag más cercana al gRNA                       
            
            #redefine la distancia
            if r >n/2:
                r=n-r

            #compara distancia con radio de encapsidación, además de estados extendidos
            if r<=d and gag[i][1]==1 and gag[j][1]==1: 
                if i not in par and j not in par:
                    if m<=cgrna:
                        vir=vir+1
                    else:
                        vlp=vlp+1
                    par+=[i,j]
    #sustrae de la matriz Gag las proteinas que formaron parejas. El tercer índice indica el axis.
    gag=np.delete(gag,par,0) 

    return ([gag,vir,vlp])




#Folding
    #Modifica el estado binario de las Gags de acuerdo con una probabilidad de folding por unidad de tiempo 
    # pfold=matriz de probabilidad de transición

def folding(gag,pfold):
	for elem in gag: 
		i=int(elem[1]) #foldstate
		r=random.random()

		if r<pfold[i][1-i]:
			elem[1]=1-i
	return gag

#Synthesis
    #Producción de grna a una tasa constante krna
    #Producción de gags a una tasa proporcional a la cantidad de gRNA

def synthesis(gag,cgrna,cmrna,kprot,krna,dt):
    cgrna=cgrna+dt*kprot
    dgag=dt*kprot*cgrna

    return ([cgrna,cmrna,dgag])

    #Asigna posiciones a las gag sintetizadas y las agrega al arreglo de gags
def localization(gag,dgag,cgrna,cmrna):
    n=cgrna+cmrna
    if abs(dgag-round(dgag)) <= 0.1:
        for i in range (round(dgag)):
            #Estado de nueva proteína
            newprot=[random.randint(0,n),random.randint(0,1)]
            #Agrega el nuevo estado al arreglo
            gag+=newprot 




def foldtest(gag):
    v=np.zeros(gag.shape[0])
    for i in range (gag.shape[0]):
        v[i]=gag[i][1] #vector ordenado de estados de gags 
    return v





############################ PROGRAMA PRINCIPAL ###########################
Gag=begin(GagIn,CgRNA,CmRNA) #GagIn: cantidad inicial de gags
nt=10000
#S=np.zeros([10,GagIn])
#v=np.zeros(Gag.shape[0])

print(GagIn)
#Iterador temporal
for t in range (nt):
    Gag=randomwalk(Gag,N,1)
    Gag=folding(Gag,Pfold)
    #S[t]=foldtest(Gag)
    [Gag,vir,vlp]=encaps(Gag,2*rprot,CmRNA,CgRNA,vir,vlp) #2do argumento = radio de encapsidación
    print(t+1,Gag.shape[0],vir,vlp,vir/(vir+vlp),vlp/(vir+vlp),CgRNA/(CgRNA+CmRNA),CmRNA/(CgRNA+CmRNA))
