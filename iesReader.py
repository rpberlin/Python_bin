#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import sys
import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
mWperLm = 1000/683.0


def iesReader(inputfilename):
    file1 = open(inputfilename, 'r')
    Lines = file1.readlines()
    nLines = len(Lines)
    for i, line in enumerate(Lines):
        if ("TILT" in line) or ("tilt" in line):
            idatastart = i+1

    vectorizedData = []
    for iLine in range(idatastart,nLines):
        for item in Lines[iLine].split(' '):
            item = item.strip()
            if item:
                vectorizedData.append(item)
    nVectorDataPts = len(vectorizedData)
    #print('vectorSize: ',nVectorDataPts)



    nLamps = int(vectorizedData[0])
    #print('nLamps: ',nLamps)

    lumensPerLamp = float(vectorizedData[1])
    #print('lumensPerLamp: ',lumensPerLamp)

    lumenMultiplier = float(vectorizedData[2])
    #print('lumenMultiplier: ',lumenMultiplier)
    if int(lumenMultiplier) < 0:
        lumenMultiplier = 1.0

    nominalMilliWattsOutput = mWperLm*nLamps*lumensPerLamp*lumenMultiplier
    #print('Total Output = ',totalmilliWattsOutput,' mW')

    nThetaPts = int(vectorizedData[3])
    nPhiPts = int(vectorizedData[4])
    #print('nThetaPts: ',nThetaPts,' nPhiPts: ',nPhiPts)

    photometricVal = int(vectorizedData[5])
    if photometricVal == 1:
        photometricType = 'Type C'
    elif photometricVal ==  2:
        photometricType = 'Other'
    else:
        photometricType = 'NOT FOUND'
    #print('photometricType: ',photometricType)

    unitsVal = int(vectorizedData[6])
    if unitsVal == 1:
        unitsType = 'feet'
    elif unitsVal ==2:
        unitsType = 'meters'
    else:
        unitsType = 'NOT_FOUND'
    #print('unitsType: ',unitsType)

    lampWidth = float(vectorizedData[7])
    lampLength = float(vectorizedData[8])
    lampHeight = float(vectorizedData[9])
    #print('Lamp::Length: ',lampLength,' Width: ',lampWidth,' Height: ',lampHeight)

    dummy1 = float(vectorizedData[10])
    dummy2 = float(vectorizedData[11])
    dummy3 = float(vectorizedData[12])
    #print('dummy1: ',dummy1,' dummy2: ',dummy2,' dummy3: ',dummy3)

    iStart = 13
    thetaVals =[]
    for i in range(0,nThetaPts):
        thetaVals.append(float(vectorizedData[i+iStart]))
    #print('thetaVals: ',thetaVals)

    phiVals =[]
    for i in range(0,nPhiPts):
        phiVals.append(float(vectorizedData[i+iStart+nThetaPts]))
    #print('phiVals: ',phiVals)

    meshgridTheta = np.zeros([nThetaPts,nPhiPts])
    meshgridPhi = np.zeros([nThetaPts,nPhiPts])
    photoNet = np.zeros([nThetaPts,nPhiPts])
    meshgridX = np.zeros([nThetaPts,nPhiPts])
    meshgridY= np.zeros([nThetaPts,nPhiPts])
    meshgridZ = np.zeros([nThetaPts,nPhiPts])
    idx = iStart+nThetaPts+nPhiPts
    #print('vectorizedData: ',vectorizedData[idx-5:idx+5])
    #print(nVectorDataPts, ' but needed are: ',idx+nThetaPts*nPhiPts)
    if(phiVals[0] < 0):
        #print("Phi Values Must Be Offset")
        tmpOffset = phiVals[0]
        for i, tmp in enumerate(phiVals):
            phiVals[i]=tmp-tmpOffset
    #print('phiVals:',phiVals)

    if(thetaVals[0] < 0):
        #print("Theta Values Must Be Offset")
        tmpOffset = thetaVals[0]
        for i, tmp in enumerate(thetaVals):
            thetaVals[i]=tmp-tmpOffset
    #print('thetaVals:',thetaVals)

    for j in range(0,nPhiPts):
        for i in range(0,nThetaPts):
            theta = float(thetaVals[i])
            phi = float(phiVals[j])
            intensity = lumenMultiplier*float(vectorizedData[idx])
            meshgridTheta[i,j] = theta
            meshgridPhi[i,j] = phi
            photoNet[i,j]= intensity
            #print(theta,phi, intensity)
            meshgridX[i,j] = intensity*np.sin(theta*np.pi/180.0)*np.sin(phi*np.pi/180.0)
            meshgridY[i,j] = intensity*np.sin(theta*np.pi/180.0)*np.cos(phi*np.pi/180.0)
            meshgridZ[i,j] = -1.0*intensity*np.cos(theta*np.pi/180.0)
            idx+=1
    #print('phiVals',phiVals,'lastphival: ',phiVals[-1])

    phiSpan = float(phiVals[-1]) - float(phiVals[0])
    if int(phiSpan)==90:
        #print('photoNet Must Be Mirrored')
        meshgridTheta,meshgridPhi,meshgridX,meshgridY,meshgridZ, photoNet,thetaVals, phiVals = mirrorPhotoNet(meshgridTheta,meshgridPhi,meshgridX,meshgridY,meshgridZ, photoNet,thetaVals,phiVals)

    phiSpan = float(phiVals[-1]) - float(phiVals[0])
    if int(phiSpan)==180:
        #print('photoNet Must Be Mirrored')
        meshgridTheta,meshgridPhi,meshgridX,meshgridY,meshgridZ, photoNet,thetaVals, phiVals = mirrorPhotoNet(meshgridTheta,meshgridPhi,meshgridX,meshgridY,meshgridZ, photoNet,thetaVals,phiVals)
    #print('photoNet: ',photoNet)
    return meshgridTheta,meshgridPhi,meshgridX,meshgridY,meshgridZ, photoNet,thetaVals, phiVals,nominalMilliWattsOutput, unitsType, photometricType, nLamps, nPhiPts,nThetaPts,lampWidth,lampHeight,lampLength, dummy1,dummy2,dummy3

def plotIESdata(meshgridTheta,meshgridPhi,meshgridX,meshgridY,meshgridZ, photoNet,thetaVals, phiVals,meshgridEyeX,meshgridEyeY,meshgridEyeIrr,meshGridFluence, maxEyeIrr, avgFluence):
    xspan = np.max(meshgridX)-np.min(meshgridX)
    yspan = np.max(meshgridY)-np.min(meshgridY)
    zspan = np.max(meshgridZ)-np.min(meshgridZ)
    pspan = np.max(photoNet)-np.min(photoNet)
    allspan = 0.5*max(xspan,yspan,zspan)
    orig_map= cm.get_cmap('viridis')

    C = np.linspace(photoNet.min(), photoNet.max(), photoNet.size).reshape(photoNet.shape)
    scamap = plt.cm.ScalarMappable(cmap='viridis')
    fcolors = scamap.to_rgba(photoNet)

    reversed_map = orig_map.reversed()
    #print(f'xspan {xspan} yspan {yspan} zspan {zspan}')

    figtmp = plt.figure()
    axtmp = figtmp.add_subplot(111,projection='3d')
    axtmp.plot_surface(meshgridX, meshgridY, meshgridZ,facecolors=fcolors)
    axtmp.view_init(elev=-20., azim=220)
    plt.title(f'{inputfilename} Output = {photoIntegral:.3f} (mW)')
    plt.savefig('tmp.png')
    plt.close(figtmp)

    fig = plt.figure(figsize=(15, 10))
    ax2 = fig.add_subplot(221)
    plt.contourf(meshgridPhi,meshgridTheta,photoNet,cmap=orig_map)
    ax2.set_xlabel('Azimuth (deg)')
    ax2.set_ylabel('Elevation From Negative Z-Axis (deg)')
    plt.title(f'Intensity Map Total Output = {photoIntegral:.3f} (mW)')
    cbar = plt.colorbar()
    cbar.set_label('(mW/sr)')



    ax = fig.add_subplot(222,projection='3d')
    ax.plot_surface(meshgridX, meshgridY, meshgridZ,facecolors=fcolors)
    ax.set_xlim(-allspan,allspan)
    ax.set_ylim(-allspan,allspan)
    ax.set_zlim(-2*allspan,0)
    plt.title(inputfilename)

    ax3 = fig.add_subplot(223)
    plt.contourf(meshgridEyeX,meshgridEyeY,meshgridEyeIrr,cmap=orig_map)
    ax3.set_xlabel('X-Direction (m)')
    ax3.set_ylabel('Y-Direction (m)')
    plt.title(f'Eye Level Intensity Map: Max = {maxEyeIrr:.3} (\u03bcW/cm$^2$)')
    cbar = plt.colorbar()
    cbar.set_label('(\u03bcW/cm$^2$)')

    ax4 = fig.add_subplot(224)
    plt.contourf(meshgridEyeX,meshgridEyeY,meshGridFluence,cmap=orig_map)
    ax3.set_xlabel('X-Direction (m)')
    ax3.set_ylabel('Y-Direction (m)')
    maxFluence = np.amax(meshGridFluence)
    plt.title(f'Irradiance Map \n Max = {maxFluence:.3f} (\u03bcW/cm$^2$)   Average = {avgFluence:.3f} (\u03bcW/cm$^2$)')
    cbar = plt.colorbar()
    cbar.set_label('(\u03bcW/cm$^2$)')


    plt.show()
    figFilename = inputfilename+'_iconized.png'
    #plt.savefig(figFilename)

    command1 = 'sips -i tmp.png'
    command2 = "DeRez -only icns tmp.png > tmpicns.rsrc"
    command3 = "Rez -append tmpicns.rsrc -o "+inputfilename
    command4 = "rm tmp.png tmpicns.rsrc"
    command5 = "SetFile -a C "+inputfilename
    #print(command1)
    os.system(command1)
    #print(command2)
    os.system(command2)
    #print(command3)
    os.system(command3)
    #print(command4)
    os.system(command4)
    #print(command5)
    os.system(command5)
    #np.savetxt("fooph.csv", photoNet, delimiter=",")
    #print('thetaVals: ',thetaVals)
    #print('phiVals: ', phiVals)
    #np.savetxt("footheta.csv", thetaVals, delimiter=",")

    #np.savetxt("foophi.csv", phiVals, delimiter=",")

def mirrorPhotoNet(oldmeshgridTheta,oldmeshgridPhi,oldmeshgridX,oldmeshgridY,oldmeshgridZ, oldphotoNet, oldthetaVals,oldphiVals):
    jPhiMirror = len(oldphiVals)-1
    phiFold = float(oldphiVals[jPhiMirror])
    nPhiPtsFull = 1+2*(len(oldphiVals)-1)
    nThetaPtsFull = len(oldthetaVals)
    meshgridThetaNew = np.zeros([nThetaPtsFull,nPhiPtsFull])
    meshgridPhiNew = np.zeros([nThetaPtsFull,nPhiPtsFull])
    photoNetNew = np.zeros([nThetaPtsFull,nPhiPtsFull])
    meshgridXNew = np.zeros([nThetaPtsFull,nPhiPtsFull])
    meshgridYNew= np.zeros([nThetaPtsFull,nPhiPtsFull])
    meshgridZNew = np.zeros([nThetaPtsFull,nPhiPtsFull])
    phiValsNew = np.zeros(nPhiPtsFull)
    thetaValsNew = np.zeros(nThetaPtsFull)

    for j in range(0,nPhiPtsFull):
        for i in range(0,nThetaPtsFull):

            if j>jPhiMirror:
                iOld = i
                jOld = 2*jPhiMirror-j
                flip= -1
                phiSubj = float(2*phiFold - float(oldphiVals[jOld]))
            else:
                iOld = i
                jOld = j
                flip = 1
                phiSubj = float(oldphiVals[jOld])

            thetaSubi = float(oldthetaVals[iOld])

            thetaValsNew[i] = thetaSubi
            phiValsNew[j] = phiSubj

            intensityij = float(oldphotoNet[iOld,jOld])
            meshgridThetaNew[i,j] = float(thetaSubi)
            meshgridPhiNew[i,j] = float(phiSubj)
            photoNetNew[i,j] = float(intensityij)
            meshgridXNew[i,j] = intensityij*np.sin(thetaSubi*np.pi/180.0)*np.sin(phiSubj*np.pi/180.0)
            meshgridYNew[i,j] = intensityij*np.sin(thetaSubi*np.pi/180.0)*np.cos(phiSubj*np.pi/180.0)
            meshgridZNew[i,j] = -1.0*intensityij*np.cos(thetaSubi*np.pi/180.0)

    return meshgridThetaNew,meshgridPhiNew,meshgridXNew,meshgridYNew,meshgridZNew, photoNetNew, thetaValsNew, phiValsNew

def fullSphereIntegral(photoNet,azimuthVals, elevationVals):
    nAzPts = len(azimuthVals)
    nElPts = len(elevationVals)
    photoIntegral = 0
    surfaceIntegral = 0

    for i in range(0,nAzPts):
        for j in range(0,nElPts):
            iminus = max(0,i-1)
            iplus = min(nAzPts-1,i+1)
            jminus = max(0,j-1)
            jplus = min(nElPts-1,j+1)
            az      = float(azimuthVals[i])
            azminus = float(azimuthVals[iminus])
            azplus  = float(azimuthVals[iplus])
            el      = float(elevationVals[j])
            elminus = float(elevationVals[jminus])
            elplus  = float(elevationVals[jplus])
            delta_az = 0.5*(az+azplus)-0.5*(az+azminus)
            delta_elev = 0.5*(el+elplus)-0.5*(el+elminus)
            dw1 = np.sin(elminus*np.pi/180)*delta_az*np.pi/180
            dw2 = np.sin(elplus*np.pi/180)*delta_az*np.pi/180
            dh = delta_elev*np.pi/180
            dA = 0.5*(dw1+dw2)*dh
            photoIntegral += photoNet[j,i]*dA*mWperLm
            surfaceIntegral += dA
    #print('Method One: photoIntegral: ',photoIntegral,' surfaceIntegral:',surfaceIntegral/(4*np.pi))
    surfaceIntegral = 0
    photoIntegral = 0
    for i in range(0,nAzPts-1):
        for j in range(0,nElPts-1):
            PointA = photoNet[j,i]
            PointB = photoNet[j+1,i]
            PointC = photoNet[j,i+1]
            PointD = photoNet[j+1,i+1]
            azminus = float(azimuthVals[i])
            azplus  = float(azimuthVals[i+1])
            elminus = float(elevationVals[j])
            elplus  = float(elevationVals[j+1])
            delta_az = azplus-azminus
            delta_elev = elplus-elminus
            dw1 = np.sin(elminus*np.pi/180)*delta_az*np.pi/180
            dw2 = np.sin(elplus*np.pi/180)*delta_az*np.pi/180
            dh = delta_elev*np.pi/180
            dA = 0.5*(dw1+dw2)*dh
            photoIntegral += mWperLm*dh*(dw1*(2*PointA + PointB + 2*PointC + PointD)+dw2*(PointA + 2*PointB + PointC + 2*PointD))/12
            surfaceIntegral += dh*(dw1*(2*1 + 1 + 2 + 1)+dw2*(1 + 2 + 1 + 2))/12
            #print(f'i: {i} j: {j} el: {el} az: {az} dA:{dA}')
    #print('Method Two: photoIntegral: ',photoIntegral,' surfaceIntegral:',surfaceIntegral/(4*np.pi))
    return photoIntegral, surfaceIntegral/(4*np.pi)


def calcEyeIntensityAndFluenceMaps(photoNet,thetaVals, phiVals,lampLength,lampWidth):
    verticalDistanceToEyePlane = 3*.3048 #vertica distance in feet converted to meters between lamp and eye-plane
    verticalDistanceToFluencePlane = 5*.3048 #vertica distance in feet converted to meters between lamp and eye-plane
    roomSqft = 400 #room Size in squareFeet
    eyeDX   = 0.03 #distance in meters between eye points
    lampDX = .05
    lampDY = .05
    roomLW = 0.5*np.sqrt(roomSqft)*.3048
    eyeHalfConeAngle = 40*np.pi/180 #HalfconeAngle of viewer
    nEyePts = int(np.ceil(2*roomLW/eyeDX))
    #print('nEyePts', nEyePts)
    if nEyePts%2 == 0:
        nEyePts+=1
    eyePts = np.linspace(-1.0*roomLW,roomLW,num=nEyePts)

    nLampPtsX = 1+int(lampLength/lampDX)
    nLampPtsY = 1+int(lampWidth/lampDY)
    nLampPts = nLampPtsX*nLampPtsY
    lampPtsX = np.linspace(-0.5*lampLength,0.5*lampLength,num=nLampPtsX)
    lampPtsY = np.linspace(-0.5*lampWidth,0.5*lampWidth,num=nLampPtsY)
    lampPts = np.zeros([nLampPts,2])
    idx = 0
    for i in range(0,nLampPtsX):
        for j in range(0,nLampPtsY):
            lampX = lampPtsX[i]
            lampY = lampPtsY[j]
            lampPts[idx,0] = lampX
            lampPts[idx,1] = lampY
            idx +=1
    #print('nLampPts: ',nLampPts,'lampPts: ',lampPts)



    meshgridEyeX = np.zeros([nEyePts,nEyePts])
    meshgridEyeY = np.zeros([nEyePts,nEyePts])
    meshgridEyeIrr = np.zeros([nEyePts,nEyePts])
    meshGridFluence = np.zeros([nEyePts,nEyePts])
    maxEyeIrr = 0
    avgFluence = 0
    for i in range(0,nEyePts):
        for j in range(0,nEyePts):
            eyeIrr =0
            localFluence= 0
            for k in range(0,nLampPts):
                deltaX = eyePts[i] - lampPts[k,0]
                deltaY = eyePts[j] - lampPts[k,1]
                deltaZeye = verticalDistanceToEyePlane
                deltaZfluence = verticalDistanceToFluencePlane
                deltaR = np.sqrt(deltaX*deltaX + deltaY*deltaY)
                elev  = np.arctan2(deltaR,deltaZeye)
                elevFluence  = np.arctan2(deltaR,deltaZfluence)
                azim = np.pi+np.arctan2(deltaY,deltaX)
                eyeAngle = np.pi*0.5-elev
                #print(f'i {i} j {j} deltaX {deltaX} deltaY {deltaY} deltaZ {deltaZ} deltaR {deltaR} elev {elev} azim {azim}')
                intensityFluence = (1/nLampPts)*interpolate2DPhoton(elevFluence*180/np.pi,azim*180/np.pi,thetaVals,phiVals,photoNet)
                distance2 = deltaR*deltaR+deltaZeye*deltaZeye
                distanceFluence2 = deltaR*deltaR+deltaZfluence*deltaZfluence
                localFluence += 0.1*intensityFluence/distanceFluence2
                if eyeAngle < eyeHalfConeAngle:
                    intensityEye = (1/nLampPts)*interpolate2DPhoton(elev*180/np.pi,azim*180/np.pi,thetaVals,phiVals,photoNet)
                    eyeIrr += 0.1*intensityEye/distance2
                else:
                    eyeIrr += 0

            avgFluence += localFluence

            meshgridEyeX[i,j] = deltaX
            meshgridEyeY[i,j] = deltaY
            meshgridEyeIrr[i,j]=eyeIrr
            meshGridFluence[i,j]=localFluence
            maxEyeIrr = max(maxEyeIrr,eyeIrr)
            #print(f'{i} in {nEyePts} {j} in {nEyePts} {k} in {nLampPts}')
    avgFluence = avgFluence/(nEyePts*nEyePts)
    #return meshGrids in meters and fluence and irradiance in uW/cm2
    return meshgridEyeX,meshgridEyeY,meshgridEyeIrr,meshGridFluence, maxEyeIrr, avgFluence

def calcFluenceMap(photoNet,thetaVals, phiVals):
    verticalDistanceToEyePlane = 3*.3048 #vertica distance in feet converted to meters between lamp and eye-plane
    roomSqft = 400 #room Size in squareFeet
    eyeDX   = 0.01 #distance in meters between eye points
    roomLW = 0.5*np.sqrt(roomSqft)*.3048
    eyeHalfConeAngle = 40*np.pi/180 #HalfconeAngle of viewer
    nEyePts = int(np.ceil(2*roomLW/eyeDX))
    #print('nEyePts', nEyePts)
    if nEyePts%2 == 0:
        nEyePts+=1
    eyePts = np.linspace(-1.0*roomLW,roomLW,num=nEyePts)
    meshgridEyeX = np.zeros([nEyePts,nEyePts])
    meshgridEyeY = np.zeros([nEyePts,nEyePts])
    meshgridEyeIrr = np.zeros([nEyePts,nEyePts])
    avgFluence = 0
    for i in range(0,nEyePts):
        for j in range(0,nEyePts):
            deltaX = eyePts[i]
            deltaY = eyePts[j]
            deltaZ = verticalDistanceToEyePlane
            deltaR = np.sqrt(deltaX*deltaX + deltaY*deltaY)
            elev  = np.arctan2(deltaR,deltaZ)
            azim = np.pi+np.arctan2(deltaY,deltaX)
            eyeAngle = np.pi*0.5-elev
            #print(f'i {i} j {j} deltaX {deltaX} deltaY {deltaY} deltaZ {deltaZ} deltaR {deltaR} elev {elev} azim {azim}')
            if eyeAngle < eyeHalfConeAngle:
                intensity = interpolate2DPhoton(elev*180/np.pi,azim*180/np.pi,thetaVals,phiVals,photoNet)
                #intensity =1
            else:
                intensity = 0
            distance2 = deltaR*deltaR+deltaZ*deltaZ
            eyeIrr = intensity/distance2

            meshgridEyeX[i,j] = deltaX
            meshgridEyeY[i,j] = deltaY
            meshgridEyeIrr[i,j]=eyeIrr

    return meshgridEyeX,meshgridEyeY,meshgridEyeIrr





def interpolate2DPhoton(theta0, phi0, thetaVals,phiVals, photoNet):
    nThetaPts = len(thetaVals)
    nPhiPts = len(phiVals)
    iLow = 0
    jLow = 0
    iHigh = 0
    jHigh = 0
    alpha = 1
    beta = 1
    if theta0 >= thetaVals[-1]:
        iLow = nThetaPts-1
        iHigh = nThetaPts-1
        alpha = 1
    elif theta0 < thetaVals[-1] and theta0 > thetaVals[0]:
        for i in range(1,nThetaPts):
            if theta0 >= thetaVals[i-1] and theta0 < thetaVals[i]:
                iLow = i-1
                iHigh = i
                alpha = (theta0-thetaVals[i-1])/(1e-6+thetaVals[i]-thetaVals[i-1])
    if phi0 >= phiVals[-1]:
        jLow = nPhiPts-1
        jHigh = nPhiPts-1
        beta = 1
    elif phi0 < phiVals[-1] and phi0 > phiVals[0]:
        for j in range(1,nPhiPts):
            if phi0 >= phiVals[j-1] and phi0 < phiVals[j]:
                jLow = j-1
                jHigh = j
                beta = (phi0-phiVals[j-1])/(1e-6+phiVals[j]-phiVals[j-1])
    return photoNet[iLow][jLow]*(1-alpha)*(1-beta) + photoNet[iHigh][jLow]*alpha*(1-beta) + photoNet[iLow][jHigh]*(1-alpha)*beta + photoNet[iHigh][jHigh]*alpha*beta










if __name__ == '__main__':
    n = len(sys.argv)
    print("Total arguments passed:", n)
    #for i in range(1, n):
    #    print(i," ",sys.argv[i])

    if n < 2:
        print("Correct Usage iesReader.py filename.ies")
    if n ==2:
        inputfilename = sys.argv[1];
        base, ext = os.path.splitext(inputfilename)
        if ext == '.ies':
            meshgridTheta,meshgridPhi,meshgridX,meshgridY,meshgridZ, photoNet,thetaVals, phiVals, nominalMilliWattsOutput, unitsType, photometricType, nLamps, nPhiPts,nThetaPts,lampWidth,lampHeight,lampLength, dummy1,dummy2,dummy3 = iesReader(inputfilename)
            photoIntegral, surfaceIntegral = fullSphereIntegral(photoNet,phiVals, thetaVals)
            print(f'Integrated Intensity: {photoIntegral:.2f} (mW) of Integrated Area: {surfaceIntegral:.2f}')
            meshgridEyeX,meshgridEyeY,meshgridEyeIrr,meshGridFluence, maxEyeIrr, avgFluence = calcEyeIntensityAndFluenceMaps(photoNet,thetaVals, phiVals,lampLength,lampWidth)
            plotIESdata(meshgridTheta,meshgridPhi,meshgridX,meshgridY,meshgridZ, photoNet,thetaVals, phiVals,meshgridEyeX,meshgridEyeY,meshgridEyeIrr,meshGridFluence, maxEyeIrr, avgFluence)
            print('TotalOutput: ',nominalMilliWattsOutput,' (claimed)',photoIntegral,'(integrated) (mW) AverageFluence: ', avgFluence,' (uW/cm2)  MaxEyeIrr: ', maxEyeIrr ,' ',maxEyeIrr ,' unitsType:',unitsType ,' photometricType:',photometricType,' nLamps: ',nLamps ,' WxHxL',lampWidth,'x',lampHeight,'x',lampLength,' Dummy: ',dummy1,dummy2,dummy3,' filename: ',inputfilename)

    if n > 2:
        for i in range (1,n):
            inputfilename = sys.argv[i]
            base, ext = os.path.splitext(inputfilename)
            if ext == '.ies':
                meshgridTheta,meshgridPhi,meshgridX,meshgridY,meshgridZ, photoNet,thetaVals, phiVals,nominalMilliWattsOutput, unitsType, photometricType, nLamps, nPhiPts,nThetaPts,lampWidth,lampHeight,lampLength, dummy1,dummy2,dummy3 = iesReader(inputfilename)
                photoIntegral, surfaceIntegral = fullSphereIntegral(photoNet,phiVals, thetaVals)
                meshgridEyeX,meshgridEyeY,meshgridEyeIrr,meshGridFluence, maxEyeIrr, avgFluence = calcEyeIntensityAndFluenceMaps(photoNet,thetaVals, phiVals,lampLength,lampWidth)
                print(f'File: {i}  TotalOutput: {nominalMilliWattsOutput:.2f} (mWnom) {photoIntegral:.2f} (mWint)\tAvgFluence: {avgFluence:.2f} (uW/cm2)\tMaxEyeIrr: {maxEyeIrr:.2f} (uW/cm2)  unitsType: {unitsType}  photometricType: {photometricType}  nLamps: {nLamps} WxHxL {lampWidth:.2f} x {lampHeight:.2f} x {lampLength:.2f} Dummy: {dummy1:.2f} {dummy2:.2f} {dummy3:.2f}  filename: {inputfilename}')




