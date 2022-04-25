import numpy as np

def cosspace(startP,endP,N=100):
    
    x = np.zeros(N); x[N-1] = endP
    x[0] = startP
    
    if endP <= startP:
        raise TypeError('End point must be greater than the start point')
    else:
        midP = (endP-startP)/2
        
        
    angleInc = np.pi/(N-1)
    
    curAngle = angleInc
    for idx in range(1,N-1):
        x[idx]=startP + midP*(1-np.cos(curAngle))
        curAngle += angleInc
 
    # for idx in range(1,int(np.ceil(N/2))):
    #     x[idx] = startP + midP*(1 - np.cos(curAngle))
    #     x[-(idx)] = x[-(idx-1)] - (x[idx+1]-x[idx])
    #     curAngle += angleInc
    
        
    return x
    