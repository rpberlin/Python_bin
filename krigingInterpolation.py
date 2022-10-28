import numpy as np
from scipy import linalg

def krigingInterpolation(Xobslocs,Yobs,Xinterps):
    Xobslocs = np.asmatrix(Xobslocs)
    Yobs = np.asmatrix(Yobs)
    Xinterps = np.asmatrix(Xinterps)

    Xobslocs = Xobslocs.T
    Yobs = Yobs.T
    Xinterps = Xinterps

    nObs, nXdims = Xobslocs.shape;
    nObsY, nYdims = Yobs.shape
    nInterps, nXdims2 = Xinterps.shape
    Yfull = np.vstack((Yobs,np.ones((1,nYdims))))

    Cij = np.zeros((nObs,nObs));
    for i in range(0,nObs):
        xi = Xobslocs[i,:]
        for j in range(i,nObs):
            xj = Xobslocs[j,:]
            Cij[i,j] = rbf_kernel(xi,xj)
            Cij[j,i] = Cij[i,j]

    Cij = np.hstack((Cij, np.ones((nObs,1))))
    #print(Cij.shape)
    bottomRow = np.ones((1,nObs))
    #print('brshape ',bottomRow.shape)
    bottomRow = np.hstack((bottomRow,np.asmatrix([0])))
    #print(bottomRow.shape)
    Cij = np.vstack((Cij, bottomRow))
    #print(Cij.shape,Cij)

    CijInv = linalg.pinv(Cij)

    lamba = CijInv*Yfull

    weights = lamba[0:nObs,:]

    Cio=np.ones((nInterps,nObs))
    for i in range(0,nInterps):
        xo = Xinterps[i,:]
        for j in range(0,nObs):
            xj = Xobslocs[j,:]
            Cio[i,j] = rbf_kernel(xo,xj)



        #if(mod(i,round(nInterps/100))==1)
        #pct = 100*i/nInterps;
        #krigT = cputime-startTime;
        #startTime = cputime;
    #print('Cio',Cio.shape)
    #print('Cijinv',CijInv.shape)
    #print('weights',weights.shape,weights , 'sumWeights = ',sum(weights))
    #print('lamba',lamba.shape)




    Yinterps = Cio*weights+lamba[-1]
    #YinterpVar = weights*weights.T
    variogram = np.empty((0, 2), float)
    for i in range(nObs):
        xi = Xobslocs[i]
        yi = Yobs[i]
        for j in range(i+1,nObs):
            xj = Xobslocs[j]
            yj = Yobs[j]
            delta = xi-xj
            d = sum(abs(delta))
            var = abs(yi-yj)
            #print('xi ',xi[0,0],' xj ',xj,' yi ',yi,' yj ', yj ,' delta ',delta,' var ',var,' d ',d)
            variogram= np.append(variogram,[[d[0,0], var[0,0]]],axis=0)



    #Yinterps = Cio*weights+lamba[len(lamba)-1];
    #size(Cio);
    #size(weights);

    #ond(Cij);
    #size(Yinterps);
    print('Yinterps',Yinterps.shape,Yinterps)
    #print('YinterpVar',YinterpVar.shape,YinterpVar)

    return Yinterps, variogram




def rbf_kernel(xi,xj):
    #print('xi',xi.shape,xi)
    #print('xj',xj.shape,xj)
    delta = xi-xj
    weight = sum(delta**2)
    return weight
