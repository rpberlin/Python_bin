
from sympy import Matrix, zeros
import matplotlib.pyplot as plt
plt.style.use('bmh')

#name L M T K
vars =[]
vars.append(['g',1,0,-2,0])
vars.append(['L',1,0,0,0])
vars.append(['m',0,1,0,0])
vars.append(['t',0,0,1,0])
vars.append(['mu',1,1,-1,0])
vars.append(['rho',-3,1,0,0])
vars.append(['A',2,0,0,0])
vars.append(['U',1,0,-1,0])
vars.append(['F',1,1,-2,0])
nVars = len(vars)
names = []
M = zeros(4,nVars)
for i in range(0,nVars):
    names.append(vars[i][0])
    for j in range(0,4):
        M[j,i]=vars[i][j+1]

kRank = M.rank()
nPi = nVars-kRank

RRef, v = M.rref()

for i in range(0,4):
    piStr = ""
    for j in range(0,nVars):
        exp=RRef[i,j]
        if exp != 0:
            piStr += names[j]+'^'+str(exp)+' '
    print (piStr)

print('names: ',names)
print('dimMatx: \n',M)
print('RRef: \n',RRef)
print('v: \n',v)
print('nPi: ',nPi)
