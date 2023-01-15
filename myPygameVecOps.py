def Pblend(P1,P2,fac):
    delta = (int(fac*(P2[0]-P1[0])),int(fac*(P2[1]-P1[1])))
    return (P1[0]+delta[0],P1[1]+delta[1])

def floatVectorAdd(vec1,vec2,factor=1.0):
    return (factor*vec2[0]+factor*vec1[0],factor*vec2[1]+factor*vec1[1])

def floatVectorScale(vec1,factor=1.0):
    return (factor*vec1[0],factor*vec1[1])

def floatVectorDiff(vec1,vec2,factor=1.0):
    return (factor*(vec2[0]-vec1[0]),factor*(vec2[1]-vec1[1]))

def intVectorAdd(vec1,vec2,factor=1):
    return (factor*vec2[0]+factor*vec1[0],factor*vec2[1]+factor*vec1[1])

def intVector(floatVector):
    return (int(floatVector[0]),int(floatVector[1]))
