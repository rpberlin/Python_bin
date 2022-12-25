#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('bmh')

def myLogSpace(fmin,fmax,nFreqs):
    alpha = np.exp(np.log(fmax/fmin)/(nFreqs-1))
    freqs = np.zeros(nFreqs)
    freqs[0]=fmin
    for i in range(1,nFreqs):
        freqs[i]=alpha*freqs[i-1]
    return freqs



def  runChirp():

    maxTime = 10
    nChirpPts = 5000
    deltaTime = maxTime/nChirpPts
    #print(deltaTime,'fMax:', 0.5/deltaTime)
    timeChirp = np.linspace(0,maxTime,nChirpPts)
    yChirp = np.zeros(nChirpPts)
    for i in range(0,nChirpPts):
        f = 20+6.5*timeChirp[i]-0.35*timeChirp[i]*timeChirp[i]
        #f = 10+2*np.sin(2*np.pi*timeChirp[i])
        yChirp[i]=np.sin(2*np.pi*f*timeChirp[i])

    return timeChirp, yChirp

def getMorletKernel(f,dt):
    k = f*f
    negK = -1.0*k
    minGauss = .0001
    maxTime = np.sqrt(-1*np.log(minGauss)/k)
    #print(maxTime)
    nMax = int(maxTime/dt)+1
    nMin = -1*nMax
    nKernelPts = 2*nMax+1
    kernel = np.zeros([nKernelPts,2])
    kernelTime = np.linspace(nMin*dt,nMax*dt,nKernelPts)
    #print(kernelTime.shape,kernelTime[0])
    gaussSum = 0
    for i in range(0,nKernelPts):
        gauss = np.exp(negK*kernelTime[i]*kernelTime[i])
        real = gauss*np.cos(2*np.pi*f*kernelTime[i])
        imag = gauss*np.sin(2*np.pi*f*kernelTime[i])
        kernel[i,0] = real
        kernel[i,1] = imag
        gaussSum += gauss

    kernel = 2.0*kernel/gaussSum


    return kernel, nKernelPts




def runWavelet(time, ysignal):

    nPts = len(time)
    dt = time[1]-time[0]
    fmax = 0.9*0.5/dt
    fmin = 10*fmax/nPts
    nFreqs = 300
    freqs = myLogSpace(fmin,fmax,nFreqs)
    waveReal = np.zeros([nPts,nFreqs]) #i: time j: frequency => waveMag[iTime, jFreq]
    waveImag = np.zeros([nPts,nFreqs]) #i: time j: frequency => waveMag[iTime, jFreq]
    waveMag = np.zeros([nPts,nFreqs]) #i: time j: frequency => waveMag[iTime, jFreq]
    waveFreqMesh = np.zeros([nPts,nFreqs])
    waveTimeMesh = np.zeros([nPts,nFreqs])
    maxFFT = np.zeros(nFreqs)

    for j, freq in enumerate(freqs):
        for i in range(0,nPts):
            waveFreqMesh[i,j]=freq
            waveTimeMesh[i,j]=time[i]

    for j,freq in enumerate(freqs):
        kernel, nKernelPts = getMorletKernel(freq,dt)
        iPlus = int((nKernelPts-1)/2)
        #print(j, freq, nKernelPts)
        if nPts > nKernelPts:
            for i in range(iPlus,nPts-iPlus):
                ysub = ysignal[i-iPlus:i+iPlus+1]
                waveReal[i,j] = np.dot(kernel[:,0],ysub)
                waveImag[i,j] = np.dot(kernel[:,1],ysub)
                waveMag[i,j] = np.sqrt(waveReal[i,j]*waveReal[i,j]+waveImag[i,j]*waveImag[i,j])
                maxFFT[j] = max(maxFFT[j],waveMag[i,j])
    return waveFreqMesh,waveTimeMesh,waveReal, waveImag, waveMag, freqs, maxFFT


def plotWavelet(waveFreqMesh,waveTimeMesh,waveReal, waveImag, waveMag, freqs, maxFFT):

    fig = plt.figure(figsize=(14, 5))
    ax2 = fig.add_subplot(121)
    orig_map= plt.cm.get_cmap('viridis')
    plt.contourf(waveTimeMesh,waveFreqMesh,waveMag,cmap=orig_map)
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Frequency (Hz)')
    cbar = plt.colorbar()

    ax2 = fig.add_subplot(122)
    ax2.semilogx(freqs,maxFFT)
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylabel('Max Wavelet Intensity (-)')
    plt.grid(True, which ='both')

    plt.show()
    return












def plotTandX(T,X,species_to_track):
    nPhiPts, nSpecies = X.shape
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    line1 = plt.plot(phiPts,T,label='Temperature')
    ax1.set_xlabel('Equivalence Ratio (-)')
    ax1.set_ylabel('Temperature (K)')
    plt.title('Temperature')
    plt.legend()

    ax2 = fig.add_subplot(122)
    for i, specie in enumerate(species_to_track):
        line2 = plt.semilogy(phiPts,X[:,i],label=specie)
    ax2.set_ylim([1e-6,1])

    plt.legend()
    plt.title('Equilibrium Mole Fractions')
    plt.show()





if __name__ == '__main__':
    timeChirp, yChirp = runChirp()
    waveFreqMesh,waveTimeMesh,waveReal, waveImag, waveMag, freqs, maxFFT = runWavelet(timeChirp,yChirp)
    plotWavelet(waveFreqMesh,waveTimeMesh,waveReal, waveImag, waveMag, freqs, maxFFT)
