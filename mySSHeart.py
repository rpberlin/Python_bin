import numpy as np
import matplotlib.pyplot as plt
plt.style.use('bmh')

class mySSHeart:
    def __init__(self,HRmin, HRmax, Pmax):
        self.delayTau = 4.582
        self.decayTau = 75
        self.HRmin = HRmin
        self.HRmax = HRmax
        self.Pmax =  Pmax
        self.HRvsP_b = .12
        self.HRvsP_m = .6
        self.loadCoeff = 3e-5
        self.loadExp = .8

    def __str__(self):
        ret = f"Pmax: {self.Pmax:.2f} HR: {self.HRmin:.2f}-{self.HRmax:.2f}\n"
        ret = ret+f"Tau Delay: {self.delayTau:.3f}: Decay: {self.decayTau:.3f}\n"
        ret = ret+f"Slope: {self.HRvsP_m:.3f} Intercept {self.HRvsP_b:.3f}\n"
        ret = ret+f"LoadCoeff: {self.loadCoeff:.3e} LoadExp: {self.loadExp:.3f}\n"
        ret = ret+f"300W Steady State HR: {self.getHRfor300W():.2f} (bpm)"
        return ret


    def HRfitter(self, HRactual,power):
        HR0 = HRactual[0]
        nSearches = 1000;
        nLineSearchSteps = 35;
        convergenceTolerance = 1e-7
        delta = 1e-6;
        lineSearchStepSize0=1;
        Xorig = self.getXfromHRclassParams()
        X0 = Xorig
        Xnorm = [.1, .1, .1, .1, 1, 1e-2, 1e-2, 1e-7, .001]
        nParams = len(X0)


        Hhat0 = self.HRsim(power, HR0)
        err0 = self.HRerr(HRactual, Hhat0)
        bestErr = err0;
        bestHhat = Hhat0;
        bestX = X0;
        
        errHist = []
        bestErrHist = []
        bestkHist = []
        convHist = []
        errHist.append(err0)
        bestErrHist.append(err0)
        
        for j in range(0,nSearches):
            errVec=np.zeros(nParams)
            
            #CalculateGradient
            for i in range(0,nParams):
                Xg=X0;
                Xg[i]=Xg[i]+delta*Xnorm[i]
                self.setHRclassParamsFromX(Xg)
                Hhat = self.HRsim(power,HR0)
                err = self.HRerr(HRactual,Hhat)
                errVec[i] = err - bestErr 
            
            #Do Line Search
            deltaVec0 = -1.0*errVec*Xnorm/(bestErr*delta)
            lineErr=np.zeros(nLineSearchSteps)
            bestk = 0;
            lineSearchStepSizes = [lineSearchStepSize0*np.power(1.5,k) for k in range(0,nLineSearchSteps)]
            k=0
            while k < nLineSearchSteps:
                Xls=X0+lineSearchStepSizes[k]*deltaVec0;
                self.setHRclassParamsFromX(Xls)
                Hhat = self.HRsim(power,HR0)
                lineErr[k] = self.HRerr(HRactual,Hhat)
                errHist.append(lineErr[k])
                if lineErr[k] < bestErr:
                    bestErr = lineErr[k]
                    bestHhat = Hhat;
                    bestX = Xls;
                    bestk=k
        
                if lineErr[k] > bestErr:
                    k=nLineSearchSteps;

                k=k+1;

            if bestk > 0 and bestk < nLineSearchSteps -1 :

                lineSearchAlpha, bestFitErr = self.parabolaMinFromThreePoints(lineSearchStepSizes[bestk-1:bestk+2],lineErr[bestk-1:bestk+2])
                Xalpha=X0+lineSearchAlpha*deltaVec0
                self.setHRclassParamsFromX(Xalpha)
                Hhat = self.HRsim(power,HR0)
                alphaErr = self.HRerr(HRactual,Hhat)
                errHist.append(alphaErr)
                if alphaErr < bestErr:
                    bestErr = alphaErr;
                    bestHhat =Hhat;
                    bestX = Xalpha;
                    

        
            
            errHist.append(bestErr)
            bestErrHist.append(bestErr)
            bestkHist.append(bestk)
            if all(item == 0 for item in bestkHist[-10:]):
                break
            
            convergence = (max(bestErrHist[-10:])-min(bestErrHist[-10:]))/bestErrHist[0]
            convHist.append(convergence)
            if convergence < convergenceTolerance:
                break
            
            print(j,bestk,bestErr, convergence)
            X0=bestX;
            self.setHRclassParamsFromX(X0)

        
        return convHist


    def HRerr(self,HRactual,HRhat):
        return np.power(sum(np.power(np.power(HRactual,3)-np.power(HRhat,3),2))/len(HRactual),1/6)
    
    def HRerrPlot(self,HRactual,HRhat):
        nPts = len(HRactual)
        HRerr = self.HRerr(HRactual,HRhat)
        time = np.linspace(0,nPts-1,nPts)
        plt.figure()
        plt.plot(time,HRactual,label='Actual')
        plt.plot(time,HRhat,label='Fitted')
        plt.xlabel('Time (-)')
        plt.ylabel('HR (bpm)')
        plt.title(f'HR @ 300W: {self.getHRfor300W():.2f} (bpm)')
        plt.legend()
        plt.show()
        
    def HRerrPlotWithBounds(self,HRactual,HRhat,power):
        nPts = len(HRactual)
        HRerr = self.HRerr(HRactual,HRhat)
        time = np.linspace(0,nPts-1,nPts)
        ones = np.ones(nPts)
        loadNorm = self.loadNormCalc(power)
        loadNormHR = self.loadNormCalc(power)*(self.HRmax - self.HRmin)
        plt.figure()
        plt.plot(time,HRactual,label='Actual')
        plt.plot(time,HRhat,label='Fitted')
        plt.plot(time,self.HRmax*ones,label = 'HR Maximum')
        plt.plot(time,loadNormHR+self.HRmin,label="HR Minimum")
        plt.xlabel('Time (-)')
        plt.ylabel('HR (bpm)')
        plt.title(f'HR @ 300W: {self.getHRfor300W():.2f} (bpm)')
        plt.legend()
        plt.show()
        
    def HRbeforeAfterPlot(self,HRactual,HRhat1,HRhat2,errHist):
        nPts = len(HRactual)
        HRerr1 = self.HRerr(HRactual,HRhat1)
        HRerr2 = self.HRerr(HRactual,HRhat2)
        time = np.linspace(0,nPts-1,nPts)
        plt.figure()
        plt.subplot(121)
        plt.plot(time,HRactual,label='Actual')
        plt.plot(time,HRhat1,label='Before')
        plt.plot(time,HRhat2,label='After')
        plt.xlabel('Time (-)')
        plt.ylabel('HR (bpm)')
        plt.title(f'HR @ 300W: {self.getHRfor300W():.2f} (bpm)')
        plt.title
        plt.legend()
        plt.subplot(122)
        plt.semilogy(errHist)
        plt.xlabel('Line Sweep Iterations')
        plt.ylabel('HR Error (bpm)')
        plt.title(f'HR Error Before: {HRerr1:.2f} After: {HRerr2:.2f}')
        plt.show()
    
    def loadNormCalc(self,power):
        loadNorm = np.array([self.loadCoeff*np.power(power_i/self.Pmax,self.loadExp) for power_i in power])
        loadNorm = loadNorm.cumsum()
        return loadNorm


    def dHRdt(self, power_i,HR,loadNorm_i):
        steadyHR = self.HRvsP_b+self.HRvsP_m*power_i/self.Pmax+loadNorm_i;
        HRnorm = (HR-self.HRmin)/(self.HRmax-self.HRmin);
        HRdot = (self.HRmax-self.HRmin)*(steadyHR-HRnorm)/self.decayTau;
        return HRdot

    def HRsim(self, power, HR0):
        nPts = len(power)
        Hhat = np.zeros(nPts)
        loadNorm = self.loadNormCalc(power)
        Hhat[0]=HR0
        for i in range(1,nPts):
            delta_HR = self.dHRdt(power[i],Hhat[i-1],loadNorm[i-1])
            Hhat[i]=Hhat[i-1]+delta_HR;
        return Hhat
    
    def setHRclassParamsFromX(self,Xtemp):
        #Xorig = [self.HRmin, self.HRmax, self.Pmax, self.delayTau, self.decayTau, self.HRvsP_b, self.HRvsP_m, self.loadCoeff, self.loadExp]
        self.HRmin = Xtemp[0]
        self.HRmax = Xtemp[1]
        self.Pmax  = Xtemp[2]
        self.delayTau  = Xtemp[3]
        self.decayTau  = Xtemp[4]
        self.HRvsP_b   = Xtemp[5]
        self.HRvsP_m   = Xtemp[6]
        self.loadCoeff = Xtemp[7]
        self.loadExp   = Xtemp[8]
        return
        
    def getXfromHRclassParams(self):
        #Xorig = [self.HRmin, self.HRmax, self.Pmax, self.delayTau, self.decayTau, self.HRvsP_b, self.HRvsP_m, self.loadCoeff, self.loadExp]
        Xtemp = np.zeros(9)
        Xtemp[0] = self.HRmin
        Xtemp[1] = self.HRmax 
        Xtemp[2] = self.Pmax
        Xtemp[3] = self.delayTau
        Xtemp[4] = self.decayTau
        Xtemp[5] = self.HRvsP_b
        Xtemp[6] = self.HRvsP_m
        Xtemp[7] = self.loadCoeff 
        Xtemp[8] = self.loadExp
        return Xtemp
    
    def parabolaMinFromThreePoints(self,X,Y):
        X1 = X[0]
        X2 = X[1]
        X3 = X[2]
        Y1 = Y[0]
        Y2 = Y[1]
        Y3 = Y[2]
        denom = (X1 - X2) * (X1 - X3) * (X2 - X3)
        A     = (X3 * (Y2 - Y1) + X2 * (Y1 - Y3) + X1 * (Y3 - Y2)) / denom
        B     = (X3*X3 * (Y1 - Y2) + X2*X2 * (Y3 - Y1) + X1*X1 * (Y2 - Y3)) / denom
        C     = (X1 * X3 * (X1 - X3) * Y1 + X3 * X1 * (X3 - X1) * Y2 + X1 * X1 * (X1 - X1) * Y3) / denom
        bestX = -B / (2*A)
        bestY = C - B*B / (4*A)
        return bestX, bestY
    
    def getHRfor300W(self):
        return self.HRmin+(self.HRmax-self.HRmin)*(self.HRvsP_b+300*self.HRvsP_m/self.Pmax)
        
        
















