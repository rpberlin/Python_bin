#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

def dNdt(Lambda,rho, Beta, N, sub_lambdas,c):
    return (rho-Beta)*N/Lambda+np.dot(sub_lambdas,c)

def dcdt(Lambda, sub_betas, N, sub_lambdas, c):
    term1 = sub_betas*N/Lambda
    term2 = sub_lambdas*c

    return term1 - term2


def neutronIntegrator(Lambda, Beta, sub_lambdas, sub_betas, N0, c0):
    print(sub_betas)
    print(sum(sub_betas))
    rho =.001
    Ndot = dNdt(Lambda, rho,Beta,N0,sub_lambdas, c0)
    cdot = dcdt(Lambda, sub_betas,N0, sub_lambdas, c0)
    a=1

    A = np.zeros([7,7])
    B = np.zeros([7,1])
    C = np.eye(7)
    D = 0

    A[0,0] = -Beta/Lambda
    B[0,0] = rho/Lambda
    for i in range(0,6):
        A[0,i+1]   = sub_lambdas[i]
        A[i+1,0]   = sub_betas[i]/Lambda 
        A[i+1,i+1] = -sub_lambdas[i]

    print('A: ',A,'B: ',B)
    




























    pass


if __name__ == '__main__':
    Lambda = 4.48274e-7
    Beta = 0.00432039
    sub_lambdas = np.array([.0127023, 0.0301099, 0.112331, 0.327449, 1.22596, 8.14883])
    sub_betas = np.array([8.78147e-5, 8.16105e-5, 6.51854e-4, 1.7707e-3, 7.88203e-4, 2.05715e-4])

    neutronIntegrator(Lambda, Beta, sub_lambdas, sub_betas, 1.0 , 0.0*sub_lambdas)

    

