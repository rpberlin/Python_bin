import numpy as np
from scipy import linalg

def myRegress(y,X):

    X = np.asmatrix(X)
    X = X.T

    y = np.asmatrix(y)
    y = y.T


    Nsamps,Nvar = X.shape;
    print('Nsamsp Nvar', Nsamps, Nvar)

    X0 = np.ones((Nsamps,1))
    #print('X0',X0)
    X = np.hstack((X,X0))
    #print('X = ', X.shape, X)
    #print('y = ', y.shape, y)

    FtF = linalg.inv(X.T*X)
    coeffs = FtF*X.T*y

    H=X*FtF*X.T;
    yhat = H*y;
    resid = y-X*coeffs;

    #print('coeffs = ', coeffs.shape, coeffs)
    #print('yhat = ',yhat.shape, yhat)

    #MSE = (1/Nsamps)*(resid.T*resid);
    #C = MSE*FtF;

    #Se = np.sqrt(np.diag(C));

    #SST = y.T*(np.eye(Nsamps)-(1/Nsamps)*np.ones(Nsamps,Nsamps))*y;
    #SSR = y.T*(H-(1/Nsamps)*np.ones(Nsamps,Nsamps))*y;
    #SSE = y.T*(np.eye(Nsamps)-H)*y;
    #MSE2 = SSE/(Nsamps-Nvar+1);
    #Tstats = abs(coeffs)/Se;
    #return coeffs, Tstats, yhat, resid
    return coeffs, yhat, resid
