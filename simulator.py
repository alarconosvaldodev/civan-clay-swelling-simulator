
"""Simulacion de yacimientos que presenta hinchamiento de arcillas

INICIO"""
#----------------------------------------------------------------------------

#PLUGINS
import numpy 
import matplotlib.pyplot as plt
import math

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#Nomenclatura
N=18

"Dimension de las matrices"
sigmapAST_1=numpy.zeros(N)
sigmap_1=numpy.zeros(N)
dsigmapASTdt_1=numpy.zeros(N)
dsigmapdt_1=numpy.zeros(N)

sigmapAST=numpy.zeros(N)
sigmap=numpy.zeros(N)
densidadpf=numpy.zeros(N)
R=numpy.ones(N)
K=numpy.zeros(N)
KtK=numpy.ones(N)
Kk=numpy.ones(N)
phisw=numpy.zeros(N+1)
phi=numpy.ones(N)
phip=numpy.zeros(N+1)
dphiswdt=numpy.ones(N)
dphipdt=numpy.ones(N)
#dphidt=numpy.ones(N)
P=numpy.ones(N)
u=numpy.zeros(N)
dpdx=numpy.zeros(N)
dsigmapASTdt=numpy.zeros(N)
dsigmapdt=numpy.zeros(N)
f=numpy.zeros(N)
F=numpy.zeros(N)
Pn=numpy.ones(N+1)
Pn2=numpy.ones(N+1)
Pn2s=numpy.ones(N)
K1=numpy.zeros(N)
K2=numpy.zeros(N)
U1=numpy.zeros(N)
U2=numpy.zeros(N)
"Condiciones iniciales"
t=0
n=0
sigmapAST[0]=23.09#3.7*10**(-1)  # i=1,2,...,N
sigmap[0]=  0   # i=1,2,...,N
densidadpf[0]=0 # i=1,2,...,N
"Boundary conditions"
Pinicial=4000

u[0]=0.000232939#7.1*10**-1
densidadpf[0]=0
"Rate data"
B=0.000749135#1.2*10**(-5)
k1=3.12475#6.4
k2=2.2324e-12#1*10**(-9)
k3=0.000212#9.5*10**(-2)
k4=7.93e-3#7.9*10**-3
k5=0.0001
k6=KtK[0]=1.7
dosAB=1e-4#1*10**(-4)
"Parameters and propierties"
R1=4.4*10**-2
phi0=0.19
K0=0.040#0.0004
D=100#0.0833#2.54
L=1000#0.16 #3000psi
dpdxcr=0
dpdxcrAST=0
viscosidad=0.8
Temp=186#f
densidadp=62.428#1
densidadl=156.07#2.5
Dx=L/N
Bc=1.127e-3
tfinal=360  #dias
Dt=100      #dias
T= tfinal/Dt
Bref=186
Vb=Dx*D
K
Pin=-2500
Pout=4100
for i in range(N):
    Pn[i]=Pn[i]*Pinicial
    Pn2[i]=Pn[i]
    Pn2s[i]=Pn2[i]
    u[i]=u[i]+u[0]
    sigmap[i]=sigmap[0]
    sigmapAST[i]=sigmapAST[0]   
"COMIENZO"
while t<=tfinal:

   t=t+Dt
   TEBM=1
   EBM=2
   toleP=1
   deltaP=2  
   while EBM >= TEBM:       
    t=t+Dt  
    iteracion=0          
    while deltaP>toleP: 
        for i in range(N):
            "Step 3"       #delta promedio de p 
       
            if i==0:
                dpdx[i]=(4*Pn2s[1]-Pn2s[2]-3*Pin)/(2*Dx)            
            elif i==N-1:                
                dpdx[N-1]=(3*Pout-4*Pn2s[N-2]+Pn2s[N-3])/(2*Dx)
            else:     
                dpdx[i]=(Pn2s[i+1]-Pn2s[i-1])/(2*Dx)       #i=2,3,...,(N-1)
         
            if abs(dpdx[i])<=dpdxcr:
                U1[i]=0
            else:
                U1[i]=1
            if abs(dpdx[i])<=dpdxcrAST:
                U2[i]=0
            else:
                U2[i]=1
        
        for i in range(N):   
            #generated and rate of particle release
            sigmapAST_1[i]=sigmapAST[i]-Dt*k3*sigmapAST[i]*(1-math.exp(-k4*t**0.5))*math.exp(-k5*sigmap[i])*U2[i]*(-dpdx[i]-(-dpdxcrAST))  #i=1,2,...,N
            dsigmapASTdt_1[i]=-k3*sigmapAST[i]*(1-math.exp(-k4*t**0.5))*math.exp(-k5*sigmap[i])*U2[i]*(-dpdx[i]-(-dpdxcrAST))  #i=1,2,...,N  

            "Step 4"
            #Deposition and rate of deposition
            sigmap_1[i]=sigmap[i]+Dt*(k1*u[i]*densidadpf[i]*phi[i]-k2*sigmap[i]*U1[i]*(-dpdx[i]-dpdxcr))
            dsigmapdt_1[i]=k1*u[i]*densidadpf[i]*phi[i]-k2*sigmap[i]*U1[i]*(-dpdx[i]-(-dpdxcr))
           
        "Step 5" 
            #phi:porosidad del medio
            #phio: Porosidad inicial
            #phisw:Porosidad
            #phip:Porosidad 
            #esfp: masa de partÃ­culas depositadas por unidad de volumen
            #densp: Densidad de partícula
            #2AB: constante fenomenolÃ³gica por hinchamiento
            #B:constante fenomenológica por absorción líquida
            #Effective porosities and permeabilities
                
        for i in range(N):
            
            Kk[i]=(k6)+(1-(k6))*math.exp(-dosAB*t**0.5)
            phip[i]=sigmap_1[i]/densidadp
            phisw[i]= -phi0*(1-Kk[i]**(1/3))
            phi[i]=phi0-phip[i]-phisw[i]  #i=1,2,...,N
            K[i]=K0*(phi[i]/phi0)**3
        
        Sp=0.0001
        #Solv the pressure equation
        for i in range(1,N):
            dphiswdt[i]=(phisw[i+1]-phisw[i-1])/(Dt)
            dphipdt[i]=(phip[i+1]-phip[i-1])/(Dt)
        
        for i in range(1,N):

                sumdphiswdt=numpy.sum(dphiswdt[i])/(N-2)

                sumdphipdt=numpy.sum(dphipdt[i])/(N-2)
        mdphiswdt=sumdphiswdt/(N-2)
        mdphipdt=sumdphiswdt/(N-2)
        dphidt=-mdphipdt-mdphiswdt
        
    
        for i in range(N-1):
            
            K1[i]=2*Bc*K[i]*K[i+1]/(K[i+1]+K[i])
            
        for i in range(1,N):
            K2[i]=2*Bc*K[i]*K[i-1]/(K[i-1]+K[i])
        
        
        for i in range(N):
            
            f[i]=viscosidad*(dphidt+(1/densidadp)*(dsigmapdt_1[i]+dsigmapASTdt_1[i])+Sp/densidadl)
           
            F[i]=(-Dx**2)*f[i]
        
            if i==0:
                F[i]=(-Dx**2)*f[i]-K1[i]*Pin
            if i==N-1:
                F[i]=(-Dx**2)*f[i]-K2[i]*Pout
              			
        W=numpy.eye(N,N,k=-1)
        break 
        E=numpy.eye(N,N,k=1)
        C=numpy.identity(N)
        A=numpy.zeros((N,N))
        for i in range(N):
            for k in range(1,N):
                             
                    W[i,k-1]=K2[i]*W[i,k-1]
                    
        for i in range(N):
            C[i,i]=-(K1[i]+K2[i]+f[i])
       
        for i in range(N):
            for k in range(N):
                    E[i,k]=K1[i]*E[i,k]
                            					
        for i in range(N):
            for k in range(N):
                    A[i,k]=C[i,k]+W[i,k]+E[i,k]		
        B=numpy.zeros((N,1))		        
        for i in range(N):
            B[i,0]=F[i]
        A_=numpy.linalg.inv(A)
        PSOL=numpy.ones((N,1))
        PSOL=A_@B  
        Pn3=numpy.zeros((N,1))
        for i in range(N):
             Pn3[i,0]=Pn2s[i]                            
        deltaP1=numpy.zeros(N)                           
        for i in range(N):
            deltaP1[i]=math.fabs(PSOL[i]-Pn3[i])#PARA RESTAR COMPONENTE A COMPONENTE 
            deltaP=numpy.zeros(N)
        for i in range(N):
            deltaP=deltaP1[i].max()           
        for i in range(N):
            Pn2s[i]=PSOL[i,0] #SE CONVIERTE EN MATRIZ PARA CAL         
      
        "Step 9"
        u_1=numpy.zeros(N)
        "Volumen de flujo"
        for i in range(N):
            if i==0:
                 u_1[i]=-K1[i]*(4*Pn2s[1]-Pn2s[2]-3*Pin)/(2*Dx*viscosidad) 
            elif i==N-1:
                 u_1[i]=-K2[i]*(3*Pout-4*Pn2s[N-2]+Pn2s[N-3])/(2*Dx*viscosidad)
            else:   
                 u_1[i] = -K1[i]*(Pn2s[i+1]-Pn2s[i-1])/(2*Dx*viscosidad)   #i=2,3,...,(N-1)        
              
        "Step 11"
        "Concentración másica"        
        b=numpy.zeros(N)
        c=numpy.zeros(N)
        d=numpy.zeros(N)
        Gamma=Dt/(2*Dx)
        
        for i in range(N):
                 
        
                 if i==0:
                      b[i] = phi[i]+Dt*dphidt+Gamma*(u_1[i]-u_1[i+1])
                 if i==N-1:
                     b[i] = phi[i]+Dt*dphidt+Gamma*(u_1[i-1]-u_1[i])  #i=2,3,...,(N-1)
                 else:
                     b[i] = phi[i]+Dt*dphidt+Gamma*(u_1[i+1]-u_1[i-1])   #i=2,3,...,(N-1)
        
        for i in range (N):              
            c[i]= Gamma*u_1[i]
            
        for i in range (N):    
        
            if i==0:
                d[i]= phi[i]*densidadpf[i]+Gamma*u_1[i]*densidadpf[i]+Dt*(dsigmapdt_1[i]+dsigmapASTdt_1[i])
            if i==N-1:
                d[i]= phi[i]*densidadpf[i]+Dt*(dsigmapdt_1[i]+dsigmapASTdt_1[i])   
               
            d[i]= phi[i]*densidadpf[i]+Dt*(dsigmapdt_1[i]+dsigmapASTdt_1[i])
           
        W=numpy.eye(N,N,k=-1)
        for i in range(N):
            for k in range(1,N):
                    W[i,k-1]=-u_1[i]*W[i,k-1]
        C=numpy.identity(N)
        for i in range(N):
            C[i,i]=-(b[i]*C[i,i]+d[i])
        E=numpy.eye(N,N,k=1)
        for i in range(N):
            for k in range(N):
                    E[i,k]=c[i]*E[i,k]                    
        X=numpy.zeros((N,N))					
        for i in range(N):
            for k in range(N):
                    X[i,k]=C[i,k]+W[i,k]+E[i,k]				        
        DD=numpy.zeros((N,1))        
        for i in range(N):
            DD[i,0]=d[i]                         
        Rho_inv=numpy.linalg.inv(X)
        Zz=numpy.ones((N,1))
        Zz=Rho_inv@DD 
        
        for zz in range (N):
            densidadpf[zz]=math.fabs(Zz[zz,0])
            u_1[i]=u[i]
        
        break 
    dVdt=numpy.zeros(N)
    suma=numpy.zeros(N)
    summa=numpy.zeros(N)
    q=-(3*Pin+Pn2s[3]-4*Pn2s[2])/(viscosidad*Dx*2/K0)  
    for i in range(N): 
                  
            dVdt[i]=Vb/t
            phin12=phi[i]*((deltaP))
            dphid=((phin12*0.0006)/Bref) 
            suma[i]= (dVdt[i]*dphid) + suma[i]                                           
    for i in range(N):
            summaACU=suma[i].max()                          
    EBM=math.fabs((summaACU/q))
    if EBM>TEBM:
        Dt=Dt+Dt
        print("no cumple criterio EBM")
    else:
        print("si cumple criterio EBM")
        break
   
"Gráficas"

Distanciax=numpy.zeros(N)
for i in range(N):
    if N==0:
        Distanciax[i]=Dx
    else:
        Distanciax[i]=Distanciax[i-1]+Dx

  
plt.plot(Distanciax,K/K[0],'o-')
plt.title("Distancia contra Permeabilidad (K/K0)") 
plt.xlabel("Distancia [ft]")               
plt.ylabel("K/K0") 
plt.show()

plt.plot(Distanciax, Pn2s,'o-')
plt.title("Distancia contra Presión final")   
plt.xlabel("Distancia [ft]")           
plt.ylabel("Presión [Psi]") 
print("PRESIÓN EN EL ESPACIO")
print(Pn2s)
plt.show()

plt.plot(Distanciax, densidadpf,'o-')
plt.title("Distancia contra Densidad de la partícula final")   
plt.xlabel("Distancia [ft]")               
plt.ylabel("Densidad [lb/ft3]") 
print("DENSIDAD DE LA PARTÍCULA")
print(densidadpf)
plt.show()

print("PERMEABILIDAD EFECTIVA")
print(K)

print("POROSIDAD EFCTIVA")
print(phi)

