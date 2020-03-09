#Created on Wed Jan 22 14:15:11 2020
#@author: aarondomenzain

import numpy as np
import random
import math
from timeit import default_timer as timer

start = timer()


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
    gagprod=0
    concgrna=0
    Ngags0=0
    Ngags1=0
    #Step en segundos (calculado a partir del coeficiente de difusión en agua a 37 C para una partícula de 40nm)
    dt=0.0000488/(dprot**2)
    return (gag,vir,vlp,allostvir,cgrna,concgrna,gagprod,dt,Ngags0,Ngags1)
    
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
    rna=rna-80*dprot
    return (rna)

#Extensión alostérica
    #Modifica un estado de par compacto-extendido a extendido-extendido según una probabilidad pfold01.  
def allosteric(ai,aj,pfold):
    allost=0
    #if ai+aj==1:
    #    print("Allosteric not happening. States:",(ai,aj))
    if ai+aj==0:
        #Test
        print("Allosteric extension is probable")
        p=random.random()
        
        if p<=pallost:
            ai=1
            aj=1
            allost=1
            print("ALLOSTERIC EXTENSION HAPPENED!!")


            #Test for debugging
            if allost!=1:
                print("ERROR IN ALLOSTERIC MODULE: Allosteric extension happened but allost counter is equal to:", allost)
        
        #if ai+aj==0:
        #    print("Allosteric compaction happened...")

    #Test for debugging
    if abs(ai+aj)>2:
        print("ERROR IN ALLOSTERIC MODULE: non binary state")
        print("States=",(ai,aj))  

    return (ai,aj,allost)


#Counter
    #Cuenta el número de gags en estado compacto o extendido
def count(gag,ngags0,ngags1):
    ngags0=0
    ngags1=0

    for prot in gag:
    
        if prot[1]==0:
            ngags0=ngags0+1
        if prot[1]==1:
            ngags1=ngags1+1

    cgags0=ngags0/(ngags0+ngags1)
    cgags1=ngags1/(ngags0+ngags1)
    return(ngags0,ngags1,cgags0,cgags1)





#Encapsidación
    #Compara distancias entre proteínas gag para determinar si sucede encapsidación
    #Si la encapsidación ocurre en mRNA, se formará un VLP
    #Si la encapsidación ocurre en gRNA, se formará un VIRIÓN
    #Argumentos: Gag ,d=radio efectivo de encapsidación , cgrna= cantidad de grna , vir=# de viriones , vlp=# de vlps
def encaps(gag,d,cmrna,cgrna,vir,vlp,allostvir):
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
            
            #Redefine la distancia considerando las condiciones periódicas del dominio espacial
            if r >n/2:
                r=n-r
            
            #Compara distancia con radio de encapsidación, además de estados extendidos
            
            if r<=d:
                if m<=cgrna and cgrna>=80*dprot:
                    (gag[i][1],gag[j][1],allost)=allosteric(gag[i][1],gag[j][1],Pfold)


                    #Si no existe suficiente
                    #if allost==1 and cgrna<80*dprot:
                        #print("Not enough gRNA for encapsidation after allosteric extension")
                        #allost=0
                
                
                if gag[i][1]+gag[j][1]==2: 
                    if i not in par and j not in par:
                        #Verifica si hay suficiente gRNA para encapsidar
                        if m<=cgrna and cgrna>=80*dprot:

                            cgrna=capture(cgrna)
                            
                        #Test de captura
                            #print("CgRNA captured")                                

                        #Conteo de viriones por extensión alostérica
                            allostvir=allostvir+allost

                            vir=vir+allost+1

                        else:
                            vlp=vlp+1
                            #Conteo de vlp por extensión alostérica

                        #Captura de mRNA
                        #cmrna=capture(cmrna)
                        #print("CmRNA captured")    
                        
                        #Vector de las etiquetas de las parejas que se formaron
                        par+=[i,j]  
                        
    #Sustrae de la matriz Gag las proteinas que formaron parejas.
    gag=np.delete(gag,par,axis=0)
    
    #Población actual de gags
    #gagprod=gag.shape[0]
    
    return ([gag,cgrna,cmrna,vir,vlp,allostvir])


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
        
        for _ in range (dgag):
        #Nueva proteína con estado aleatorio
            newprot=[random.randint(0,n),random.randint(0,1)]
        #Agrega la nueva proteína estado al arreglo
            gag=np.append(gag,[newprot],0)
        gagprod=gagprod-int(gagprod) #or gagprod=gagprod%1
        #print("gagprod restante despues de localizar:",gagprod)
    return (gag,gagprod)
        



#LISTA DE PARÁMETROS
print()
print("Ingrese los parámetros")
#Diámetro promedio de proteína en unidades de pasos de random walk
dprot=3    
#Cantidad de moléculas de mRNA inicial, en donde cada molécula es aproximadamente 30 diámetros de gag extendida
CmRNA=(30*dprot)*int(input("Cantidad de moléculas de mRNA inicial: "))
#Cantidad de gRNA inicial
CgRNAin=0
#Cantidad de gag inicial
GagIn=0
#Tamaño inicial del dominio espacial
N=CgRNAin+CmRNA
#Tasa de producción de Gag
Kprot=float(input("Tasa de producción de Gag: "))
#Tasa de producción de gRNA (Número de unidades de gRNA de longitud igual a 80 diámetros de gag extendida)
KRNA=(80*dprot)*float(input("Tasa de producción de gRNA en moléculas/s: "))
#Energía libre entre estados extendido-compacto en unidades kbT (valor experimental BdG=11)
BdG=float(input("Energía libre entre estados extendido-compacto en unidades kbT (valor experimental BdG=11): ")) 
#Matriz de probabilidad de transición
pfold01=float(input("Probabilidad de transición compacto -> extendido: "))
Pfold=[[1-pfold01,pfold01],[pfold01*math.exp(-BdG),1-pfold01*math.exp(-BdG)]] 
#Probabilidad de extensión alostérica
pallost=float(input("Probabilidad de extensión alostérica cerca del gRNA: "))




############################ PROGRAMA PRINCIPAL #####################################

#Declaración de arreglo inicial
(Gag,vir,vlp,allostvir,CgRNA,concgrna,Gagprod,dt,Ngags0,Ngags1)=begin(GagIn,CgRNAin,CmRNA) 

#Número de iteraciones en el tiempo
tsim=float(input("Tiempo de simulación en s: "))
nt=int(tsim/dt)


data=np.zeros([1,9])

for t in range (nt):
    [CgRNA,CmRNA,Gagprod,deltagag]=synthesis(GagIn,Gagprod,CgRNA,CmRNA,vir,vlp,Kprot,KRNA,dt)
    (Gag,Gagprod)=localization(Gag,Gagprod,CgRNA,CmRNA)
    
    #Si no hay gags, no hay nada que hacer
    if Gag.shape[0]<1:
        continue
    
    else:
        Gag=randomwalk(Gag,CgRNA+CmRNA,1)
        Gag=folding(Gag,Pfold)
        (Gag,CgRNA,CmRNA,vir,vlp,allostvir)=encaps(Gag,dprot,CmRNA,CgRNA,vir,vlp,allostvir) #2do argumento = radio de encapsidación
        (Ngags0,Ngags1,Cgags0,Cgags1)=count(Gag,Ngags0,Ngags1)

    if Gag.shape[0]>=1 and (t)%1000==0 or t==nt:
        
        #La concentración se mide en partes por millón, por lo que la longitud se divide en el número de unidades dprot que le corresponda a cada molécula de RNA
        concgrna=round((CgRNA/80*dprot)/((CgRNA/80*dprot)+(CmRNA/30*dprot)),5)
        gagpop=Gag.shape[0]
        if vir+vlp==0:
            selec=0
        else:
            selec=vir/(vir+vlp)

        state=[t*dt,gagpop,concgrna,vir,vlp,allostvir,selec,Cgags0,Cgags1]
        data=np.append(data,[state],axis=0)
        print(round(100*t/nt,3),round(t*dt,5),gagpop,concgrna,vir,vlp,allostvir,round(selec,3),round(Cgags0,3),round(Cgags1,3))

#Elapsed time in s
dt = timer() - start

#Lista de parámetros
param=[dprot,Kprot,KRNA,Kprot,CgRNAin,CmRNA,GagIn,N,BdG,pfold01,pallost,dt]

#Guarda el archivo de datos
np.savetxt("dataGAG.txt", data)
np.savetxt("Parameters.txt",param)



print ("Simulated in %f s" % dt)













