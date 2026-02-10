import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
plt.style.use('bmh')




class mySSACD:
    def __init__(self,ACD=None,sPathRefData=None,zElevRefData=None,Powerrefdata=None,mass_rider = 85, mass_bike = 8):
        self.ACD = 0.4
        self.mass_rider = mass_rider
        self.mass_bike = mass_bike
        self.sPathRefData = sPathRefData
        self.zElevRefData = zElevRefData
        self.Powerrefdata = Powerrefdata
        self.speedRefData = self.get_speed_from_sPath()


    def get_speed_from_sPath(self):
        speed = np.zeros_like(self.sPathRefData)
        nPts = len(self.sPathRefData)
        speed[0] = self.sPathRefData[1]-self.sPathRefData[0]
        for i in range(1,nPts):
            speed[i] = self.sPathRefData[i]-self.sPathRefData[i-1]
        return speed

    def __str__(self):
        ret = f"Best Fit ACD: {self.ACD:.3f}\n"
        return ret
    

    def ACDfitter_scipy(self):
        x0 = self.ACD
        res = minimize(self.ACDerr_standalone, x0, method='nelder-mead',
               options={'xatol': 1e-4, 'disp': True})
        return res

    def speed_err(self,speed_ref,speed_hat):
        return np.power(sum(np.power(speed_ref - speed_hat,2))/len(speed_hat),1/2)

    def plot_fitted_ACD(self):
        self.speed_fit_plot(self.speedRefData,self.ACD[0])
        return 

    def speed_fit_plot(self,speed_ref, ACD):
        speed_hat = self.ACDsim(ACD)
        nPts = len(speed_ref)
        fit_error = self.speed_err(speed_ref,speed_hat)
        time = np.linspace(0,nPts-1,nPts)
        plt.figure()
        plt.plot(time,3.6*speed_ref,label='Reference')
        plt.plot(time,3.6*speed_hat,label='Fitted')
        plt.xlabel('Time (-)')
        plt.ylabel('Speed (km/h)')
        plt.title(f'ACD = {ACD:.3f}  Error: {fit_error:.3f} (kmh_rms)')
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
        plt.plot(time/60,HRactual,label='Actual')
        plt.plot(time/60,HRhat,label='Fitted')
        plt.plot(time/60,self.HRmax*ones,label = 'HR Maximum')
        plt.plot(time/60,loadNormHR+self.HRmin,label="HR Minimum")
        plt.xlabel('Time (min)')
        plt.ylabel('HR (bpm)')
        plt.title(f'HR @ 300W: {self.getHRfor300W():.2f} (bpm) :: Power at 185bpm: {self.getSSPowFor185BPM():.2f}')
        plt.legend()
        plt.show()

    def HRerrPlotWith300WBounds(self,HRactual,HRhat,power):
        nPts = len(HRactual)
        HRerr = self.HRerr(HRactual,HRhat)
        time = np.linspace(0,nPts-1,nPts)
        ones = np.ones(nPts)
        loadNorm = self.loadNormCalc(power)
        loadNormHR = self.loadNormCalc(power)*(self.HRmax - self.HRmin)
        plt.figure()
        plt.plot(time/60,HRactual,label='Actual')
        plt.plot(time/60,HRhat,label='Fitted')
        plt.plot(time/60,self.getHRfor300W()*ones,label = 'HR 300W')
        plt.plot(time/60,loadNormHR+self.HRmin,label="HR Minimum")
        plt.xlabel('Time (min)')
        plt.ylabel('HR (bpm)')
        plt.title(f'HR @ 300W: {self.getHRfor300W():.2f} (bpm) :: Power at 185bpm: {self.getSSPowFor185BPM():.2f}')
        plt.legend()
        plt.show()

    def speed_beforeAfterPlot(self,speed_hat):
        nPts = len(HRactual)
        HRerr1 = self.HRerr(HRactual,HRhat1)
        HRerr2 = self.HRerr(HRactual,HRhat2)
        time = np.linspace(0,nPts-1,nPts)
        plt.figure()
        plt.subplot(121)
        plt.plot(time/60,HRactual,label='Actual')
        plt.plot(time/60.,HRhat1,label='Before')
        plt.plot(time/60.,HRhat2,label='After')
        plt.xlabel('Time (min)')
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


    def ACDsim(self, ACD):
        speed_hat = np.zeros_like(self.speedRefData)
        speed_hat[0]=self.speedRefData[0]
        mass_total = self.mass_bike+self.mass_rider
        self.ACD = ACD
        C1 = 0.25 * 1.18 * self.ACD 
        C2 = 0.5 * mass_total

        for i in range(1,len(self.speedRefData)):
            delta_t = 1
            delta_x = self.sPathRefData[i]-self.sPathRefData[i-1]
            power_term = delta_t * 0.5*(self.Powerrefdata[i] + self.Powerrefdata[i-1])
            grav_term = mass_total * 9.81 * (self.zElevRefData[i] - self.zElevRefData[i-1])
            V1_sq = speed_hat[i-1]*speed_hat[i-1]
            drag_term  = C1*V1_sq*delta_x
            KE_term    = C2*V1_sq 
            numerator = max(.01,power_term - grav_term - drag_term + KE_term)
            denominator = C1 *delta_x + C2 
            speed_hat[i] = np.sqrt(numerator/denominator)

            delta_speed = speed_hat[i] - speed_hat[i-1]
            a=1
        return speed_hat

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

    def getSSPowFor185BPM(self):
        return (self.Pmax/self.HRvsP_m)*(((185-self.HRmin)/(self.HRmax-self.HRmin))-self.HRvsP_b)
    
    def ACDerr_standalone(self, X0):
        self.ACD = X0
        speed_hat = self.ACDsim(self.ACD)
        speed_err = self.speed_err(self.speedRefData,speed_hat)
        print('ACD:',self.ACD[0],'Err: ',speed_err)
        return speed_err 
