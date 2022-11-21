#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3


import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
plt.style.use('bmh')

def runSim():
    L = 40     #length of duct (m)
    W = .07      #width of duct (m)
    thick = .01 #thickness of sheet metal
    Cpg = 1200  #specific heat of gas (J/kg
    Cpw = 800   #specific heat of steel J/kg
    gamma = 1.4 #ratio of gas specific heats
    MWgas = 28.7 #gas molecular weight
    rhow = 8000 #density of steel kg/m**3
    mdot_scfm = 220 #SCFM
    mdot = mdot_scfm*np.power(.3048,3)*1.17/60 #mass flow in kg/s
    Rgas = 8314/MWgas
    viscgas = 1.8e-5; #%Gas Viscosity
    kgas = .024; #%W/mK
    Re = mdot/(viscgas*W)
    Nu = .027*np.power(Re,0.8)*np.power(0.7,0.35)

    hconv=100 #%W/m^2K
    hconv = kgas*Nu/W #
    hamb = 5 #%W/m^2K

    T_init = 300 #
    T_amb = 300 #
    T_hot = 1300
    T_hot_period = 600 #
    T_hot_freq = 2*3.141/T_hot_period


    timesteps = 1800; #%number of time seconds;
    delta_t = 1;
    #time = 0 : delta_t : delta_t*(timesteps-1); %time step vector
    time  = np.linspace(0, delta_t*timesteps, timesteps+1)


    N = 10 #%Number of subdivisition
    dx = L/N

    alpha = W*hconv*dx/(mdot*Cpg)
    C1 = (2.0-alpha)/(2.0+alpha)
    C2 = (2.0*alpha)/(2.0+alpha)

    betac = hconv/(rhow*Cpw*thick)
    betaamb = hamb/(rhow*Cpw*thick)

    A=np.zeros([N,N])
    B=np.zeros([N,2])

    for i in range(0,N):
      B[i,0]=betac*0.5*(np.power(C1,i)+np.power(C1,(i-1)));
      B[i,1]=betaamb
      for j in range(0,i-1):
        A[i,j]=0.5*betac*C2*(np.power(C1,(i-1-j))+np.power(C1,(i-j)));
      A[i,i]=C2*betac*0.5-betac-betaamb;




    C = np.zeros([1,N])

    for i in range(0,N):
        C[0,i]=C2*np.power(C1,(N-i-1))

    #%Add all Xs to Youtput
    #C=[C;eye(N)];
    C = np.append(C,np.eye(N), axis=0)



    #%C[N-1]=1;

    D = np.zeros([N,2])
    D[0,0] = np.power(C1,N)

    U = np.zeros([timesteps+1,2])
    for i in range(0,timesteps):
        U[i,0] = T_hot
        U[i,1] = T_amb;
        #print(N," ",U[i,:])


    X0=np.ones([N,1])*T_init;

    #X0=linspace(T_amb,T_hot,N)

    Cgas = np.zeros([N,N]);
    for i in range(0,N):
        for j in range(0,i):
            Cgas[i,j] = C2*np.power(C1,(i-j))

    Cgas = np.append(np.zeros([1,N]),Cgas, axis=0) #%spacer for initial gas temperature

    C=np.append(C,Cgas,axis=0)

    Dgas = np.zeros([N+2,2]);
    for i in range(0,N):
        Dgas[i,0]=np.power(C1,i)
        
    
    Dgas[0,0]= 1

    D=np.append(D,Dgas,axis=0)
    
    print(A.shape,B.shape,C.shape,D.shape)

    sys=signal.StateSpace(A,B,C,D);


    tout, yout, xout =signal.lsim(sys, U,time,X0)
    
    return yout
    #plot(t,X,t,Y,t,U(:,1))

    #Y(length(t))

    #sisosys=ss(A,B(:,1),C,D(1))
    #bode(sisosys)

    #Ywall=Y(:,2:N+1);
    #Ygas=Y(:,N+2:2*N+2);




if __name__ == '__main__':
    runSim()
