#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import colebrook_pipeflow as cPipe

### CONSTANTS AND 
gravity     = 9.81          #m/s^2      @300 C 
R_g         = 8315         #J/mol-K
P_ref       = 101325        #Atmospheric Pressure



### MATERIAL PROPERTIES for NaK ###
rho_NaK     = 806          #kg/m^3     @300 C interpolate from Wikipedia
drho_dT_NaK = -.0002356 #kg/m^3-K    @300 C slope of Wikipedia Data
cp_NaK      = 912         #J/kg-K     @300 C 
mu_NaK      = 2.35e-4     #Pa-s       @300 C
dmu_dT_NaK  = -4.18e-7    #Pa-s/K     @300 C linear fit between 300-550 C  
k_NaK  =    22.4            #W/m-K        @300 C
T_freeze_NaK = -12.8        #C from Wikipedia 
Pr_NaK = cp_NaK*mu_NaK/k_NaK


### MATERIAL PROPERTIES for N2 ###
MW_N2   =   28              #g/mol
rho_N2     = P_ref*MW_N2/(R_g*(300+273.15))    #kg/m^3   @300 C ifrom Ideal Gas Law
rho_N2star = P_ref*MW_N2/(R_g*(377+273.15))    #kg/m^3   @550 C ifrom Ideal Gas Law
drho_dT_N2 = (rho_N2star- rho_N2)/(377-300)    #kg/m^3-K slope between 300-377 C  from Ideal Gas Law
cp_N2      = 1040         #J/kg-K     @300 C    From Engineering Toolbox
mu_N2      = 2.8e-5       #Pa-s       @300 C    From Engineering Toolbox
mu_N2star  = 3.3e-5       #Pa-s       @377 C    From Engineering Toolbox
dmu_dT_N2  = (mu_N2star-mu_N2)/(377-300)    #Pa-s/K     @300 C linear fit between 300-377 C  
k_N2       = .045         #W/m-K        @300    From Engineering Toolbox
Pr_N2 = cp_N2*mu_N2/k_N2



deltaL  = 10


#pipeDiam = .2 #insideDiameter of the pipe
#pipeA = 0.25*np.pi*np.power(pipeDiam,2)  #pipe Cross Sectional Area

QdotShutdown = 180e3 #Heat Rejection
QdotNominal = 30e3 #Heat Rejection

ThotNominal = 371 #sodium temperatue in deg C under nominal conditions
ThotShutdown = 377 #sodium temperatue in deg C under shutdown conditions
Tambient = 20 #ambient temperatue in deg C

nPipePts = 20
nHpts = 6
pipeDiams = cPipe.myLogSpace(.03125,1,nPipePts)
deltaHs   = cPipe.myLogSpace(1,100,nHpts)

TreturnShutdownGrid = np.zeros([nPipePts,nHpts])
TreturnNominalGrid = np.zeros([nPipePts,nHpts])
TreturnShutdownGridN2 = np.zeros([nPipePts,nHpts])
TreturnNominalGridN2 = np.zeros([nPipePts,nHpts])
dGrid = np.zeros([nPipePts,nHpts])
dHGrid = np.zeros([nPipePts,nHpts])


for j, deltaH in enumerate(deltaHs):    
    for i, pipeDiam in enumerate(pipeDiams):
        totalPipeLength = 2*(deltaH+deltaL)
        mdotShutdown, Ushutdown, TreturnShutdown, rhoShutdown, muShutdown, f0           = cPipe.solveBuoyantFlowLoop(QdotShutdown, ThotShutdown, rho_NaK, drho_dT_NaK, cp_NaK, mu_NaK, dmu_dT_NaK, 1.0, pipeDiam, totalPipeLength, deltaH,.03)
        mdotNominal, UNominal, TreturnNominal, rhoNominal, muNominal, f0                = cPipe.solveBuoyantFlowLoop(QdotNominal, ThotNominal, rho_NaK, drho_dT_NaK,  cp_NaK, mu_NaK, dmu_dT_NaK, 1.0, pipeDiam, totalPipeLength, deltaH,.03)
        mdotShutdownN2, UshutdownN2, TreturnShutdownN2, rhoShutdownN2, muShutdownN2, f0 = cPipe.solveBuoyantFlowLoop(QdotShutdown, ThotShutdown, rho_N2, drho_dT_N2,  cp_N2, mu_N2, dmu_dT_N2, 7, pipeDiam, totalPipeLength, deltaH,.03)
        mdotNominalN2, UNominalN2, TreturnNominalN2, rhoNominalN2, muNominalN2, f0      = cPipe.solveBuoyantFlowLoop(QdotNominal, ThotNominal, rho_N2, drho_dT_N2,    cp_N2, mu_N2, dmu_dT_N2, 7, pipeDiam, totalPipeLength, deltaH,.03)
        print(f'Diameter: {pipeDiam:.3f}\t MassFlow: {mdotShutdown:.3f}\t{mdotShutdownN2:.3f}\t Return Temperature: {TreturnShutdown:.3f}\t{TreturnShutdownN2:.3f}')
        dGrid[i,j]=pipeDiam 
        dHGrid[i,j]=deltaH
        TreturnShutdownGrid[i,j]=max(TreturnShutdown,T_freeze_NaK)
        TreturnNominalGrid[i,j]=max(TreturnNominal,T_freeze_NaK)
        TreturnShutdownGridN2[i,j]=max(TreturnShutdownN2,0)
        TreturnNominalGridN2[i,j]=max(TreturnNominalN2,0)
        
fig = plt.figure()
ax3 = fig.add_subplot(421)
for j in range(0,len(deltaHs),round(len(deltaHs)/4)):
    ax3.semilogx(pipeDiams, TreturnShutdownGrid[:,j], '.-', label=r' $\Delta $_$H}$: '+f'{deltaHs[j]:.2f}m')
plt.xlabel('Pipe Diameter (m)')
plt.ylabel('T$_{return}$  ($^\circ$C)')
plt.title('Shutdown T$_{return}$ ($^\circ$C)')
plt.legend()

ax3 = fig.add_subplot(422)
for j in range(0,len(deltaHs),round(len(deltaHs)/4)):
    ax3.semilogx(pipeDiams, TreturnNominalGrid[:,j], '.-', label=r' $\Delta $_$H}$: '+f'{deltaHs[j]:.2f}m')
plt.xlabel('Pipe Diameter (m)')
plt.ylabel('T$_{return}$ ($^\circ$C)')
plt.title('Nominal T$_{return}$  ($^\circ$C)')
plt.legend()

ax3 = fig.add_subplot(423)
for j in range(0,len(deltaHs),round(len(deltaHs)/4)):
    ax3.semilogx(pipeDiams, TreturnShutdownGridN2[:,j], '.-', label=r' $\Delta $_$H}$: '+f'{deltaHs[j]:.2f}m')
plt.xlabel('Pipe Diameter (m)')
plt.ylabel('T$_{return}$  ($^\circ$C)')
plt.title('N2 Shutdown T$_{return}$ ($^\circ$C)')
plt.legend()

ax3 = fig.add_subplot(424)
for j in range(0,len(deltaHs),round(len(deltaHs)/4)):
    ax3.semilogx(pipeDiams, TreturnNominalGridN2[:,j], '.-', label=r' $\Delta $_$H}$: '+f'{deltaHs[j]:.2f}m')
plt.xlabel('Pipe Diameter (m)')
plt.ylabel('T$_{return}$ ($^\circ$C)')
plt.title('N2 Nominal T$_{return}$  ($^\circ$C)')
plt.legend()


ax1 = fig.add_subplot(425)
plt.contourf(dGrid,dHGrid,TreturnShutdownGrid)
ax1.set_xlabel('Diameter (m)')
ax1.set_ylabel('DeltaH (m)')
plt.colorbar()

ax2 = fig.add_subplot(426)
plt.contourf(dGrid,dHGrid,TreturnNominalGrid)
ax2.set_xlabel('Diameter (m)')
ax2.set_ylabel('DeltaH (m)')
plt.colorbar()

ax1 = fig.add_subplot(427)
plt.contourf(dGrid,dHGrid,TreturnShutdownGridN2)
ax1.set_xlabel('Diameter (m)')
ax1.set_ylabel('DeltaH (m)')
plt.colorbar()

ax2 = fig.add_subplot(428)
plt.contourf(dGrid,dHGrid,TreturnNominalGridN2)
ax2.set_xlabel('Diameter (m)')
ax2.set_ylabel('DeltaH (m)')
plt.colorbar()

fig.tight_layout()
plt.show()



print('done')