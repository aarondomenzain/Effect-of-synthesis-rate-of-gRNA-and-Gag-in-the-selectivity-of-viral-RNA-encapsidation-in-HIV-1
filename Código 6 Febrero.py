#Created on Wed Jan 22 14:15:11 2020
#@author: aarondomenzain

import numpy as np
import random
import math

#LOCALIZACION INICIAL DE GAG
#Define el arreglo inicial de proteínas gag, en donde el i-ésimo renglón contiene la posición y el estado de folding de la i-ésima proteína.
def begin(gagin,cgrna,cmrna):
    #asigna posiciones y estados iniciales aleatorios
    n=cgrna+cmrna
    gag=np.zeros([gagin,2])
    for i in range (gagin):
        gag[i]=[random.randint(0,n),random.randint(0,1)] 
    #devuelve la matriz de posiciones y estados de gag
    return gag
    
#Random Walk
    #Cambia en +-d unidades la posición de los random walker de manera aleatoria
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
    #vector vacío que contendrá las etiquetas de las gags que formaron pareja
    par=[] 
    #Tamaño del dominio espacial
    n=cmrna+cgrna
    #Número de gags
    ngags=gag.shape[0]

    for i in range (ngags-1):
        for j in range (i+1,ngags):
            #Distancia entre gags
            r=abs(gag[i][0]-gag[j][0])
            #Posición de gag más cercana al gRNA
            m=min(gag[i][0],gag[j][0])                        
            
            #Redefine la distancia considerando las condiciones periódicas
            if r >n/2:
                r=n-r

            #Compara distancia con radio de encapsidación, además de estados extendidos
            if r<=d and gag[i][1]==1 and gag[j][1]==1: 
                if i not in par and j not in par:
                    if m<=cgrna:
                        vir=vir+1
                        #Captura de grna de 10 rprot de longitud
                        if cgrna<10*rprot:
                            cgrna=0
                        else:
                            cgrna=cgrna-10*rprot
                    else:
                        vlp=vlp+1

                    par+=[i,j]
    #Sustrae de la matriz Gag las proteinas que formaron parejas. 
    gag=np.delete(gag,par,axis=0) 
    return ([gag,cgrna,vir,vlp])




#Folding
    #Modifica el estado binario de las Gags de acuerdo con una probabilidad de folding por unidad de tiempo 
    # pfold=matriz de probabilidad de transición

def folding(gag,pfold):
    for elem in gag: 
        i=int(elem[1]) #foldstate
        r=random.random()
        if i==0 or i==1:
            if r<pfold[i][1-i]:
                elem[1]=1-i
        #test for debugging
        else:
            print("folding state error i = ",i)
            
    return gag


#Synthesis
    #Producción de grna a una tasa constante krna
    #Producción de gags a una tasa proporcional a la cantidad de gRNA existente
    #dt es el tiempo transcurrido (en segundos) en cada salto
def synthesis(gagin,cgrnain,vir,vlp,kprot,krna,n,dt): 
    
    if n==0:
        cgrna=cgrnain
        gag0=0
    else:
        #Cantidad de gag sintetizada en el tiempo anterior
        cgrna=cgrnain+int((n-1)*dt*krna)-10*vir*rprot
        gag0=(n-1)*dt*kprot*cgrna 
    #Incremento de grna    
    deltagrna=n*dt*krna 
    #Suma la parte entera de deltagrna
    cgrna=cgrnain+int(deltagrna)-10*rprot*vir
    #Cantidad de gag sintetizada en el tiempo actual
    gag1=n*dt*kprot*cgrna
    #Número de gags sintetizadas
    deltagag=gag1-gag0
    #print(deltagag,int(deltagag))
    #Población de gags al tiempo n
    gagpop=gagin+int(deltagag)-2*(vir+vlp)
    #print("synthesized = ",cgrna)
    return ([cgrna,gagpop]) #Devuelve las nuevas cantidades sintetizadas



#Localization
	#Asigna posiciones a las gag sintetizadas y las agrega al arreglo de gags
def localization(gag,gagpop,cgrna,cmrna):
    n=cgrna+cmrna
    dgag=gagpop-gag.shape[0]
    #if abs(deltagag-round(deltagag))<0.1:
    for i in range (dgag):
        #Nueva proteína con estado aleatorio
        newprot=[random.randint(0,n),random.randint(0,1)]
        #Agrega la nueva proteína estado al arreglo
        gag=np.append(gag,[newprot],0)
    return gag
        





#LISTA DE PARÁMETROS
    
#Tasa de producción de Gag
Kprot=0.1 
#Tasa de producción de gRNA
KRNA=1
#Tamaño inicial de la recta numérica discreta 
N=100000 
#Cantidad de gRNA inicial
CgRNAin=0 
#Cantidad de mRNA inicial
CmRNA=N-CgRNAin 
#Cantidad de gag inicial
GagIn=0
#radio promedio de proteína en unidades de pasos de random walk
rprot=5 
#cantidad de viriones
vir=0 
#cantidad de vlp
vlp=0 
#Energía libre entre estados extendido-compacto en unidades kbT (valor experimental BdG=11)
BdG=1 
#Matriz de probabilidad de transición
pfold01=0.1
Pfold=[[1-pfold01,pfold01],[pfold01*math.exp(-BdG),1-pfold01*math.exp(-BdG)]] 

############################ PROGRAMA PRINCIPAL #####################################

#Declaración de arreglo inicial
Gag=begin(GagIn,CgRNAin,CmRNA) 

#Iterador temporal
nt=5000
for t in range (nt):
    [CgRNA,Gagpop]=synthesis(GagIn,CgRNAin,vir,vlp,Kprot,KRNA,t,1)
    Gag=localization(Gag,Gagpop,CgRNA,CmRNA)
    Gag=randomwalk(Gag,CgRNA+CmRNA,1)
    Gag=folding(Gag,Pfold)
    [Gag,CgRNA,vir,vlp]=encaps(Gag,2*rprot,CmRNA,CgRNA,vir,vlp) 
    state=[t,Gag.shape[0],CgRNA,vir,vlp]
    #print(t+1,Gag.shape[0],CgRNA,vir,vlp) #,vir/(vir+vlp),vlp/(vir+vlp),CgRNA/(CgRNA+CmRNA),CmRNA/(CgRNA+CmRNA))
    print(state)