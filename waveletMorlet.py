#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('bmh')

def  runChirp():

    maxTime = 10
    nChirpPts = 1000
    deltaTime = maxTime/nChirpPts
    #print(deltaTime,'fMax:', 0.5/deltaTime)
    timeChirp = np.linspace(0,maxTime,nChirpPts)
    yChirp = np.zeros(nChirpPts)
    for i in range(0,nChirpPts):
        f = 20+1.5*timeChirp[i]-0.125*timeChirp[i]*timeChirp[i]
        #f = 10+2*np.sin(2*np.pi*timeChirp[i])
        yChirp[i]=np.sin(2*np.pi*f*timeChirp[i])

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    line1 = plt.plot(timeChirp,yChirp,label='chirp')

    ax1.set_xlabel('Time (s)')
    ax1.set_xlabel('Time(s)')
    ax1.set_ylabel('Signal')
    plt.title('Inputs')
    plt.legend()
    plt.show()
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

    kernel = kernel/gaussSum
    #fig = plt.figure()
    #ax1 = fig.add_subplot(121)
    #line1 = plt.plot(kernelTime,kernel[:,0],label='real')
    #line2 = plt.plot(kernelTime,kernel[:,1],label='imag')
    #plt.legend()
    #plt.show()

    return kernel, nKernelPts




def runWavelet(time, ysignal):

    nPts = len(time)
    dt = time[1]-time[0]
    fmax = 0.5/dt
    fmin = fmax/nPts
    waveReal = np.zeros([nPts,nPts]) #i: time j: frequency => waveMag[iTime, jFreq]
    waveImag = np.zeros([nPts,nPts]) #i: time j: frequency => waveMag[iTime, jFreq]
    waveMag = np.zeros([nPts,nPts]) #i: time j: frequency => waveMag[iTime, jFreq]
    waveFreqMesh = np.zeros([nPts,nPts])
    waveTimeMesh = np.zeros([nPts,nPts])

    for j in range(0,nPts):
        for i in range(0,nPts):
            waveFreqMesh[i,j]=fmin*(j+1)
            waveTimeMesh[i,j]=time[i]

    for j in range(0,nPts):
        freq = fmin*(j+1)
        kernel, nKernelPts = getMorletKernel(freq,dt)
        iPlus = int((nKernelPts-1)/2)
        #print(j, freq, nKernelPts)
        if nPts > nKernelPts:
            for i in range(iPlus,nPts-iPlus):
                ysub = ysignal[i-iPlus:i+iPlus+1]
                waveReal[i,j] = np.dot(kernel[:,0],ysub)
                waveImag[i,j] = np.dot(kernel[:,1],ysub)
                #print(i,j)

    waveMag = np.sqrt(np.multiply(waveReal,waveReal)+np.multiply(waveImag,waveImag))
    fig = plt.figure()
    plt.contourf(waveTimeMesh,waveFreqMesh,waveMag)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    cbar = plt.colorbar()
    #plt.imshow(waveMag)
    plt.show()
    return waveReal, waveImag, waveMag











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
    waveReal, waveImag, waveMag = runWavelet(timeChirp,yChirp)
