#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import numpy as np
import sys
import matplotlib.pyplot as plt

def TLVcalcFromSpectrumCSV(inputfilename):
    # Using readlines()
    file1 = open(inputfilename, 'r')
    Lines = file1.readlines()

    count = 0
    wavelengths  =[]
    intensities  =[]
    eyeWeights   =[]
    skinWeights  =[]
    sumIxTLVeye     = 0
    sumTLVeye       = 0
    sumIxTLVskin    = 0
    sumTLVskin      = 0
    sumIxlamba      = 0
    sumI            = 0
    maxI            = -1
    for line in Lines:
        count +=1
        if count > 1:
            pair = line.split(",")
            lamba = float(pair[0].strip())
            intensity = float(pair[1].strip())
            if intensity > maxI:
                maxI = intensity
            eyeTLV, skinTLV = TLVInterpolateTable(lamba)
            sumI += intensity
            sumTLVeye += eyeTLV
            sumTLVskin += skinTLV
            sumIxlamba += intensity*lamba
            sumIxTLVeye += intensity*eyeTLV
            sumIxTLVskin += intensity*skinTLV
            wavelengths.append(lamba)
            intensities.append(intensity)
            eyeWeights.append(eyeTLV)
            skinWeights.append(skinTLV)
            #print('Wavelength: ',lamba,' intensity: ',intensity,' eyeTLV: ',eyeTLV,' skinTLV: ',skinTLV)
    #print(wavelengths)
    nPts = len(wavelengths)
    deltaLamb = (wavelengths[-1]-wavelengths[0])/(nPts-1)
    lambaBar = sumIxlamba/sumI
    halfWidth = deltaLamb*sumI*2/maxI
    eyeTLVBar = sumIxTLVeye/sumI
    skinTLVBar = sumIxTLVskin/sumI
    eyeTLVatLambaBar, skinTLVatLambaBar = TLVInterpolateTable(lambaBar)

    squareWaveX = [lambaBar-0.5*halfWidth, lambaBar-0.5*halfWidth,lambaBar+0.5*halfWidth, lambaBar+0.5*halfWidth]
    squareWaveY = [0, 0.5*maxI, 0.5*maxI, 0]

    intensityRef =0
    for i in range(0,nPts-1):
        if  lambaBar >= wavelengths[i] and lambaBar < wavelengths[i+1]:
            x1 = wavelengths[i]
            x2 = wavelengths[i+1]
            y1 = intensities[i]
            y2 = intensities[i+1]
            intensityRef = y1+(y2-y1)*(lambaBar-x1)/(x2-x1)

    print(f'Mean Wavelength: {lambaBar:.2f} HalfWidth: {halfWidth:.2f} Eye TLV: {eyeTLVBar:.2f} Skin TLV: {skinTLVBar:.2f}')
    fig2 = plt.figure(2)
    #ax1 = fig2.add_subplot(122)
    ax1 = fig2.add_subplot(111)
    line4 = ax1.semilogy(wavelengths,eyeWeights,'k')
    line5 = ax1.semilogy(wavelengths,skinWeights,'r')
    line6 = ax1.scatter(lambaBar,eyeTLVBar)
    line7 = ax1.scatter(lambaBar,skinTLVBar)

    #plt.legend()

    #ax1 = fig2.add_subplot(121)
    ax1a = ax1.twinx()
    line1 = ax1a.plot(wavelengths,intensities)
    line2 = ax1a.plot(squareWaveX,squareWaveY,'-.',linewidth=0.5)
    line3 = ax1a.plot([lambaBar, lambaBar],[0,maxI])

    ax1a.set_xlim([lambaBar-1.5*halfWidth, lambaBar+1.5*halfWidth])
    #ax1a.set_ylim([0, maxI])
    ax1a.set_xlabel('Wavelength (nm)')
    ax1a.set_ylabel('Intensity (A.U.)')
    ax1.set_ylabel('Dose (mJ/cm2)')
    #ax1a.set_xlim([lambaBar-1.5*halfWidth, lambaBar+1.5*halfWidth])
    lines = [line1,line2,line3,line4,line5,line6,line7]
    labels = ['Spectral Intensity',f'HalfWidth: {halfWidth:.2f} (nm)',f'Mean Wavelength: {lambaBar:.2f} (nm)','eyeDose','skinDose',f'eyeTLV: {eyeTLVBar:.2f} (mJ/cm2)',f'skinTLV: {skinTLVBar:.2f} (mJ/cm2)']
    #abels = [l.get_label() for l in lines]
    plt.legend(labels,loc=2)
    plt.show()

def TLVInterpolateTable(wavelengthNM):
    nmVals = [180,  200,    207,    210,    230,    235,    240,    250,    260,    270,    280,    300,    310,    320,    340,    360,    380,    400]
    eyeVals =[1626, 1626,   1626,   1023,   46.8,   21.6,   10,     7,      4.6,    3,      3.4,    10,     200,    2900,   11000,  23000,  47000, 100000]
    skinVals=[10000, 10000, 3802,   2512,   158,    79.4,   39.8,   10,     10,     10,     10,     10,     200,    2900,   11000,  23000,  47000, 100000]
    nPts = len(nmVals)
    if wavelengthNM <= nmVals[0]:
        eyeTLV = eyeVals[0]
        skinTLV = skinVals[0]
        return eyeTLV, skinTLV

    if wavelengthNM >= nmVals[-1]:
        eyeTLV = eyeVals[-1]
        skinTLV = skinVals[-1]
        return eyeTLV, skinTLV

    for i in range(0,nPts-1):
        if wavelengthNM >= nmVals[i] and wavelengthNM < nmVals[i+1]:
            x1 = nmVals[i]
            x2 = nmVals[i+1]
            y1 = eyeVals[i]
            y2 = eyeVals[i+1]
            z1 = skinVals[i]
            z2 = skinVals[i+1]

            eyeTLV = y1*np.exp(np.log(y2/y1)*(wavelengthNM-x1)/(x2-x1))
            skinTLV = z1*np.exp(np.log(z2/z1)*(wavelengthNM-x1)/(x2-x1))

            return eyeTLV, skinTLV











if __name__ == '__main__':
    n = len(sys.argv)
    print("Total arguments passed:", n)
    for i in range(1, n):
        print(i," ",sys.argv[i])



    if n < 2:
        print("Correct Usage TLVcalc.py filename.csv")
    inputfilename = sys.argv[1];
    TLVcalcFromSpectrumCSV(inputfilename)
