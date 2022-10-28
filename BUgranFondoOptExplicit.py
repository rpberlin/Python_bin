#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import sys
import math
from parseGPX import parseGPX
from plot_PMdata import calcPMdata
from compare2GPXroutes import getQtyAtS
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import scipy

def granFondoOptExplicit(gpxFilename, basePower, myMass, bikeMass, CdA):

    sPath, zElev, power, cad, hr = importGPX(gpxFilename)
    mass = myMass+bikeMass;
    eta = 0.5*1.18*CdA;
    dt =1;
    slopeFilterDistance = 30;
    nPoints=length(zElev);
    HRmin = 73;
    HRmax = 187;
    Pmax = 300;
    NFrontMin = 36;
    NFrontMax = 52;
    NRearMin = 11;
    NRearMax = 30;
    MinCadence = 60; %Minimum Cadence at Minimum Gear Ratio
    MaxCadence = 85;
    MinSpeedKPH = 3.6*(1/60)*2.14*MinCadence*NFrontMin/NRearMax;
    MaxSpeedKPH = 3.6*(1/60)*2.14*MaxCadence*NFrontMax/NRearMin;

    %Filter Profile
    zFiltered = zElev;
    for i in range(0,nPts):
        sMin = sPath[i]-0.5*slopeFilterDistance
        sMax = sPath[i]+0.5*slopeFilterDistance
        zFiltered[i]=0.5*(getQtyAtS(sMin,sPath,zElev)+getQtyAtS(sMax,sPath,zElev))
    fig = plt.figure()
    ax1 = fig.add_subplot(611)
    line0 = plt.plot(sPath/1000),zElev,label='Raw'
    line1 = plt.plot(sPath/1000,zFiltered,label='Filtered')
    plt.show()



    %Simulate the Race
    maxDist = sPath(end);
    dist = 0;
    iStep = 1;
    Uo=0;
    speeds(1)=Uo;
    HRo = HRmin;
    HRX2 = 0;
    while dist < maxDist
        sPos(iStep) = dist;
        zPos(iStep) = getElevAtS(dist,sPath,zFiltered);
        z0=zPos(iStep);
        slopesPct(iStep)=100*(getElevAtS(dist+slopeFilterDistance,sPath,zFiltered)-z0)/slopeFilterDistance;
        [powerFactor, aeroFactor]= powerSlopeFactor(slopesPct(iStep));
        powers(iStep)=basePower*powerFactor;
        speeds(iStep)=Uo;
        Uo=speeds(iStep);
        HR(iStep) = HRo;
        HRo=HR(iStep);
        [HRdot, HRX2dot] = dHdt(HRo,HRmin, HRmax, powers(iStep), Pmax,HRX2);
        HR1=HRo+dt*HRdot;
        HRX2 = HRX2+HRX2dot;
        dzhat = getElevAtS(dist+Uo*dt,sPath,zFiltered)-z0;
        U1hat= abs(sqrt(Uo*Uo+(2/mass)*(powers(iStep)*dt-mass*9.81*dzhat-eta*aeroFactor*Uo*Uo*Uo*dt)));
        dzhat2 = getElevAtS(dist+U1hat*dt,sPath,zFiltered)-z0;
        U1= abs(sqrt(Uo*Uo+(2/mass)*(powers(iStep)*dt-mass*9.81*dzhat2-eta*aeroFactor*U1hat*U1hat*Uo*dt)));
        dist = dist+ U1*dt;
        Uo=U1;
        HRo=HR1;
        iStep = iStep+1;
        hSlope = waitbar(dist/maxDist);
    end
    iStep=iStep-1;
    close(hSlope)

    dS=diff(sPos);
    dWeightKPH=sum(dS.*dS)/sum(dt*dS)*3.6;
    %Find Climbs and Calculate Statistics
    climbSlopeThreshold = 2;
    minClimbVert = 15;
    climbVerts = zeros(1,1);
    climbSlopes = zeros(1,1);
    climbDistances = zeros(1,1);
    climbDurations = zeros(1,1);
    climbPowers  = zeros(1,1);
    climbKMStarts  = zeros(1,1);
    climbKMEnds  = zeros(1,1);
    climbStartsMin  = zeros(1,1);
    climbEnds  = zeros(1,1);
    climbVAMs = zeros(1,1);
    climbKMHstart = zeros(1,1);
    climbKMHend = zeros(1,1);
    nClimbs = 0;
    for i=2:iStep
    %start Nth Climb
        if slopesPct(i) > climbSlopeThreshold && (slopesPct(i-1) < climbSlopeThreshold | i == 2)
            climbStartZ = zPos(i);
            climbStartS = sPos(i);
            climbStarti = i;
        end
    %close Nth Climb
        if slopesPct(i) < climbSlopeThreshold && slopesPct(i-1) > climbSlopeThreshold
            climbEndZ = zPos(i);
            climbEndS = sPos(i);
            climbEndi = i;
            deltaZ = climbEndZ - climbStartZ;
            if deltaZ > minClimbVert
                nClimbs = nClimbs+1;
                climbVerts(nClimbs) = deltaZ;
                climbKMStarts(nClimbs) = climbStartS/1000;
                climbKMEnds(nClimbs) = climbEndS/1000;
                climbStartsMin(nClimbs) = climbStarti/60;
                climbEnds(nClimbs) = climbEndi/60;
                climbDistances(nClimbs) = climbEndS-climbStartS;
                climbDurations(nClimbs) = (climbEndi-climbStarti)/60;
                climbPowers(nClimbs) = sum(powers(climbStarti:climbEndi))/(climbEndi-climbStarti);
                climbSlopes(nClimbs) = 100*deltaZ/climbDistances(nClimbs);
                climbVAMs(nClimbs) = deltaZ*3600/(climbEndi-climbStarti);
                climbKMHstart(nClimbs) = climbStartS*3.600/climbStarti;
                climbKMHend(nClimbs) = climbEndS*3.600/climbEndi;
            end

        end
    end
    climbTable = array2table([climbKMStarts' climbKMEnds' climbKMHstart' climbKMHend' climbStartsMin' climbEnds' climbVerts' climbSlopes' climbPowers' climbDurations' climbVAMs'],'VariableNames',{'startKM' 'endKM' 'startKPH' 'endKPH' 'startMin' 'endMin' 'VertM' 'SlopesPCT' 'PowerW' 'TimeMIN' 'VAM'});


    totalVert = sum(max(diff(zFiltered),0))
    distance = sPath(end);
    kilometers = distance/1000;
    elapsedTime = iStep*dt;
    elapsedHours = elapsedTime/3600
    AverageSpeed = distance/elapsedTime
    AverageSpeedKPH = AverageSpeed*3.6
    AverageSpeedMPH = AverageSpeedKPH/1.609
    kCal = sum(powers)*dt/(0.2*4184)

    AvgPower = sum(powers)/iStep

    maxSlope = max(slopesPct)

    VAMfilter1 = 100;
    VAMfilter2 = 200;
    for i=1:iStep
        iMin1=max(1,i-VAMfilter1);
        iMax1=min(iStep,i+VAMfilter1);
        iMin2=max(1,i-VAMfilter2);
        iMax2=min(iStep,i+VAMfilter2);
        VAM1(i)=3600*max((zPos(iMax1)-zPos(iMin1))/(dt*(iMax1-iMin1)),0);
        VAM2(i)=3600*max((zPos(iMax2)-zPos(iMin2))/(dt*(iMax2-iMin2)),0);
    end

    maxVAM1 =max(VAM1)
    maxVAM2 =max(VAM2)


    figure(1)
    subplot(6,1,1)
    plot(sPath/1000,zElev,sPath/1000,zFiltered);
    ylabel('Elevation');
    subplot(6,1,2);
    plot(sPos/1000,speeds*3.6,'.')
    ylabel('Speed kmh');
    subplot(6,1,3);
    plot(sPos/1000,slopesPct,'.',sPos/1000,0*slopesPct)
    ylabel('Slope %');
    subplot(6,1,4)
    plot(sPos/1000,powers)
    ylabel('Power (W)');
    subplot(6,1,5)
    plot(sPos/1000,VAM1,'.',sPos/1000,VAM2,'-.')
    ylabel('VAM (m/hr)');
    subplot(6,1,6)
    plot(sPos/1000,HR)
    ylabel('Heart Rate (bpm)');
    xlabel('Distance (km)');
    ylim([HRmin HRmax]);

    figure(2)
    % ppslopes=[-.06:.01:.1];
    ppslopes=[floor(min(slopesPct))/100:0.01:ceil(max(slopesPct))/100];
    for i=1:size(ppslopes,2)
        [PPpowerFactor,PPaeroFactor] = powerSlopeFactor(100*ppslopes(i));
        PPPower(i) = basePower*PPpowerFactor;
        PPACd(i) = CdA*PPaeroFactor;
        PPU(i)=bikePowerSpeedStead(PPPower(i),mass,1.18*PPACd(i),ppslopes(i));
    end
    subplot (2,2,1)
    plot(ppslopes,PPPower)
    ylabel('SteadyState Power (W)');
    subplot(2,2,2)
    plot(ppslopes,PPACd)
    subplot(2,2,3)
    plot(ppslopes,PPU*3.6,ppslopes,ppslopes*0+MinSpeedKPH,ppslopes,ppslopes*0+MaxSpeedKPH)
    ylabel('SteadyState Speed (kmh)');
    subplot(2,2,4)
    plot(ppslopes,max(0,1000*PPU*3.6.*ppslopes));
    ylabel('SteadyState VAM')



    figure(3)
    powerCurve = powerCurveCalc(powers,dt);
    % plot(sPos,'.')
    subplot(3,1,1)
    semilogx(powerCurve(:,1)/60,powerCurve(:,2))
    ylabel('Power W');
    xlabel('Minutes');

    subplot(3,1,2)
    binEdges = -6.5:1:10;
    hSlope=histogram(slopesPct,binEdges);
    histogram('BinCounts',hSlope.Values/60,'BinEdges',binEdges)
    xlabel('Slope (%)');
    ylabel('time (min)');

    subplot(3,1,3)
    HRbinEdges = HRmax*(0.5:0.05:1.0);
    hhr=histogram(HR,HRbinEdges);
    histogram('BinCounts',hhr.Values/60,'BinEdges',HRbinEdges);
    xlabel('HR (BPM)');
    ylabel('time (min)');


    toc
    end



    function [sPath,zElev] = s_and_z_from_latLong(latitude,longitude,elev)
    Rearth = 6371000;
    nPts =length(latitude);
    sPath = zeros(nPts,1);
    zElev = zeros(nPts,1);
        for j = 2:nPts
            dx = Rearth*cos(deg2rad(latitude(j)))*deg2rad(longitude(j)-longitude(j-1));
            dy = Rearth*deg2rad(latitude(j)-latitude(j-1));
            sPath(j)=sPath(j-1)+sqrt(dx*dx+dy*dy);
        end
    zElev = elev;
    end


    function [U1, deltaTime, totalPower] = bikePowerSolver(basePower,Uo,deltaS, deltaZ,slope, mass, eta)
    %eta = 0.5*rho*Cd*A
    U1 =1;
    deficit = 1000;
    deltaU1 = 0;
    powerFactor = powerSlopeFactor(slope);
    totalPower = basePower*powerFactor;
    while abs(deficit) > .1
        U1 = max(U1 + 0.001*deltaU1,.000001);
    %     Ubar = sqrt(0.5*(U1*U1+Uo*Uo));
        Ubar = 0.5*(U1+Uo);
        dt = deltaS*0.5*(1/U1+1/Uo);
        EnergyIn = totalPower*dt+0.5*mass*Uo*Uo;
        EnergyOut = mass*9.81*deltaZ+eta*Ubar*Ubar*deltaS+0.5*mass*U1*U1;
        deficit = (EnergyIn-EnergyOut)/(abs(EnergyIn)+abs(EnergyOut));
    %     dEoutdU1 = 2*Ubar*eta*deltaS;
    %     deltaU1 = deficit/dEoutdU1;
        deltaU1 = U1*deficit;
        deltaZ;
        deltaS;
        Uo;
    end
    U1;
    deltaTime = deltaS/Ubar;

    end



    function [powerFactor,aeroFactor] = powerSlopeFactor(slope)

    powerTable = [ -4  0;-2 .3; 0 1.0;1  1.05; 2  1.1;3  1.15;4   1.2;6  1.2;8  1.2;10 1.3];
    aeroTable  = [ -4  1;-2 0.9; 0 0.9;1  0.9; 2  1.0;3  1.00;4   1.0;6  1.0;8  1.0;10 1.0];

    [nn , blah]= size(powerTable);
    powerFactor = powerTable(1,2);
    aeroFactor = aeroTable(1,2);
    for i=1:nn-1
        if slope >= powerTable(i,1)
        powerFactor = powerTable(i,2);
        aeroFactor = aeroTable(i,2);
            if slope < powerTable(i+1,1)
            x1 = powerTable(i,1);
            x2 = powerTable(i+1,1);
            y1 = powerTable(i,2);
            y2 = powerTable(i+1,2);
            z1 = aeroTable(i,2);
            z2 = aeroTable(i+1,2);
            powerFactor = y1+(slope-x1)*(y2-y1)/(x2-x1);
            aeroFactor  = z1+(slope-x1)*(z2-z1)/(x2-x1);
            end
        end
        if slope >=powerTable(end,1)
            powerFactor = powerTable(end,2);
            aeroFactor = aeroTable(end,2);
        end

    end
    %hardcode aeroFactor to Fulgaz value
    % aeroFactor=1;
    end


    function [powerCurve] = powerCurveCalc(powers,dt)
    filterSizes = [10,20,30,60,120,180,240,300,400,500,600,800,1000,1200,60*25,60*30,60*45,60*60,60*60*1.5,60*60*2];
    nIntervals = length(filterSizes);
    idx = 1;
    for  idx =1:nIntervals
        h = waitbar(i/nIntervals);
        filterWidth = min(length(powers),round(filterSizes(idx)/dt));
        powerCurve(idx,1)=filterWidth;
        powerCurve(idx,2)=1;

        for j=1:length(powers)-filterWidth-1
             powerbar=mean(powers(j:j+filterWidth-1));
             powerCurve(idx,2) = max(powerbar,powerCurve(idx,2));
        end
        idx=idx+1;
    end
    close(h)

    end





    % function [Hdot, X2dot] = dHdt(H, Hmin, Hmax,P, Pmax,X2)
    % a1 = .0113;
    % a2 = .0072;
    % a3 = .0041;
    % a4 = .0049;
    % a5 = 19.8002;
    %
    % X1 = (H-Hmin)/Hmax;
    % u = P/Pmax;
    % phi = a4*X1/(1+exp(-1*(X1-a5)));
    %
    % X1dot = -a1*X1+a2*X2+a2*u;
    % X2dot = -a3*X2+phi;
    % Hdot = X1dot*Hmax;
    %
    %
    %
    % end





if __name__ == '__main__':
    n = len(sys.argv)
    print("Total arguments passed:", n)
    for i in range(1, n):
        print(i," ",sys.argv[i])

    if n < 2:
        print("Correct Usage granFondoOptExplicit.py filename.gpx")

    basePower = 200
    myMass = 84
    bikeMass = 7.5
    CdA = 0.4
    inputfilename = sys.argv[1]

    if n ==2:


    granFondoOptExplicit(inputfilename, basePower, myMass, bikeMass, CdA)
