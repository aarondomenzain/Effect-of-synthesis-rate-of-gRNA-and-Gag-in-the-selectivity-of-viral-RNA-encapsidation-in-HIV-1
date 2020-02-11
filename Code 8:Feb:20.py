#Created on Wed Jan 22 14:15:11 2020
#@author: aarondomenzain

import numpy as np
import random
import math
import matplotlib.pyplot as plt 


#LOCALIZACION INICIAL DE GAG
#Asigna posiciones y estados iniciales aleatorios
#Prepara las condiciones iniciales de cantidad de gags, CgRNA, CmRNA, viriones, vlp, dintinguiendo entre éstos últimos los formados por extensión alostérica.
def begin(gagin,cgrnain,cmrna):
    #Cantidad inicial de CgRNA
    cgrna=cgrnain 
    #Tamaño inicial del dominio espacial
    n=cgrna+cmrna
    #Tamaño inicial del arreglo de gags
    gag=np.zeros([gagin,2])
    
    #Genera la matriz de posiciones y estados de gag
    for i in range (gagin):
        gag[i]=[random.randint(0,n),random.randint(0,1)] 
    
    #Población inicial de viriones, vlp y viriones formados por extensión alostérica
    vir=0
    vlp=0
    allostvir=0
    allostvlp=0
    gagprod=0
    #Step en segundos
    dt=0.0000488/(dprot**2)
    return (gag,vir,vlp,allostvir,allostvlp,cgrna,gagprod,dt)
    
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



#Captura de RNA de 10*dprot de longitud. 
    #Módulo anidado en Encapsidación.
def capture(rna):
    rna=rna-10*dprot
    return (rna)

#Extensión alostérica
    #Modifica un estado de par compacto-extendido a extendido-extendido según una probabilidad pfold01.  
def allosteric(ai,aj,pfold):
    allost=0
    if ai+aj==1:
        #Test
        #print("Allosteric extension is probable")
        p=random.random()
        
        if p<=pfold[0][1]:
            if ai==0:
                ai=ai+1 
            else:
                aj=aj+1
            allost=1
            #print("ALLOSTERIC EXTENSION HAPPENED!!")
        
        if p<=pfold[1][0]:
            if ai==1:
                ai=ai-1
            else:
                aj=aj-1
            #print("Allosteric compaction happened...")
            
    return (ai,aj,allost)
#Encapsidación
    #Compara distancias entre proteínas gag para determinar si sucede encapsidación
    #Si la encapsidación ocurre en mRNA, se formará un VLP
    #Si la encapsidación ocurre en gRNA, se formará un VIRIÓN
    #Argumentos: Gag ,d=radio efectivo de encapsidación , cgrna= cantidad de grna , vir=# de viriones , vlp=# de vlps
def encaps(gag,d,cmrna,cgrna,vir,vlp,allostvir,allostvlp):
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
            
            if r<=d:
                (gag[i][1],gag[j][1],allost)=allosteric(gag[i][1],gag[j][1],Pfold)
                
                if gag[i][1]+gag[j][1]==2: 
                    if i not in par and j not in par:
                        if m<=cgrna and cgrna>=10*dprot:
                            vir=vir+1
                            cgrna=capture(cgrna)
                        #Test de captura
                            #print("CgRNA captured")                                

                        #Conteo de viriones por extensión alostérica
                            allostvir=allostvir+allost

                        else:
                            vlp=vlp+1
                            #Conteo de vlp por extensión alostérica
                            allostvlp=allostvlp+allost
                        #Captura de mRNA
                        #cmrna=capture(cmrna)
                        #print("CmRNA captured")    
                        
                        #Vector de las etiquetas de las parejas que se formaron
                        par+=[i,j]  
                        
    #Sustrae de la matriz Gag las proteinas que formaron parejas.
    gag=np.delete(gag,par,axis=0)
    
    #Población actual de gags
    #gagprod=gag.shape[0]
    
    return ([gag,cgrna,cmrna,vir,vlp,allostvir,allostvlp])


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
            print("folding state error i=",i)
            
    return gag


#Synthesis
    #Producción de grna a una tasa constante krna
    #Producción de gags a una tasa proporcional a la cantidad de gRNA existente
def synthesis(gagin,gagprod,cgrna,cmrna,vir,vlp,kprot,krna,dt): 
    #Cantidad de CgRNA sintetizada
    cgrna=cgrna+dt*krna
    #Cantidad de CmRNA 
    cmrna=cmrna 
    #Incremento en la población de gags
    deltagag=dt*kprot*cgrna
    #Población de gags
    gagprod=gagprod+deltagag
    return ([cgrna,cmrna,gagprod,deltagag]) #Devuelve las nuevas cantidades sintetizadas



#Localization
	#Asigna posiciones a las gag sintetizadas y las agrega al arreglo de gags
def localization(gag,gagprod,cgrna,cmrna):
    #print("gagprod no localizada:",gagprod)
    if int(gagprod)>=1:
    
        #print("gagprod a localizar:",int(gagprod))
    
    #Tamaño del dominio espacial
        n=int(cgrna)+cmrna
    #Número de gags sintetizadas a agregar
        dgag=int(gagprod)
        
        for i in range (dgag):
        #Nueva proteína con estado aleatorio
            newprot=[random.randint(0,n),random.randint(0,1)]
        #Agrega la nueva proteína estado al arreglo
            gag=np.append(gag,[newprot],0)
        gagprod=gagprod-int(gagprod) #or gagprod=gagprod%1
        #print("gagprod restante despues de localizar:",gagprod)
    return (gag,gagprod)
        



#LISTA DE PARÁMETROS

#Diámetro promedio de proteína en unidades de pasos de random walk
dprot=5     
#Tasa de producción de Gag
Kprot=1
#Tasa de producción de gRNA (Número de unidades de gRNA de longitud igual a 10 diámetros de gag)
KRNA=1000*(10*dprot)
#Tamaño inicial de la recta numérica discreta 
N=10000
#Cantidad de gRNA inicial
CgRNAin=0
#Cantidad de mRNA inicial
CmRNA=N-CgRNAin 
#Cantidad de gag inicial
GagIn=0
#Energía libre entre estados extendido-compacto en unidades kbT (valor experimental BdG=11)
BdG=11 
#Matriz de probabilidad de transición
pfold01=0.1
Pfold=[[1-pfold01,pfold01],[pfold01*math.exp(-BdG),1-pfold01*math.exp(-BdG)]] 




############################ PROGRAMA PRINCIPAL #####################################

#Declaración de arreglo inicial
(Gag,vir,vlp,allostvir,allostvlp,CgRNA,Gagprod,dt)=begin(GagIn,CgRNAin,CmRNA) 

#Iterador temporal
nt=20000

#state=[]

virpop=[]
vlppop=[]
popgag=[]
graf=np.zeros([1,8])

for t in range (nt):
    [CgRNA,CmRNA,Gagprod,deltagag]=synthesis(GagIn,Gagprod,CgRNA,CmRNA,vir,vlp,Kprot,KRNA,dt)
    (Gag,Gagprod)=localization(Gag,Gagprod,CgRNA,CmRNA)
    
    #Si no hay gags, no hay nada que hacer
    if Gag.shape[0]<1:
        continue
    
    else:
        Gag=randomwalk(Gag,CgRNA+CmRNA,1)
        Gag=folding(Gag,Pfold)
        (Gag,CgRNA,CmRNA,vir,vlp,allostvir,allostvlp)=encaps(Gag,dprot,CmRNA,CgRNA,vir,vlp,allostvir,allostvlp) #2do argumento = radio de encapsidación
    
    if Gag.shape[0]>=1 and (t)%1000==0 or t==nt:
        state=[t,round(t*dt,5),Gag.shape[0],round(CgRNA/(CgRNA+CmRNA),5),vir,vlp,allostvir,allostvlp]
        graf=np.append(graf,[state],axis=0)
        print(state)

print("Elapsed time in s:",nt*dt)
#print(graf)

plt.figure(0)
plt.plot(graf[1:,1],graf[1:,4],label="Viriones")
plt.legend()
plt.title("Población de virus")

#plt.figure(0)
plt.plot(graf[1:,1],graf[1:,5],label="VLP")
plt.legend()
#plt.title("Población de VLP")
plt.show()
#print()
#print(graf[1][1])


