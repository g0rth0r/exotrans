import numpy as np

#Return mean molecular weight
def getMeanW(molecules):
    return np.prod(molecules,axis=1).sum()   
    
import const as c



def pressure(P0,H,z): #Pressure for atmopshere in hydrostatic equilibrium z in cm
    return P0 * np.exp(-z/H)

def pressure_reci(P0,H,P): #reciprical of the pressure eq. Will return z in fct of pressure
    return -H * np.log(P/P0)
    
def density(P,T): #return the atmosphere density based on the current pressure, based on IGL
    return (P*1e6)/(c.kb * T) #here pressure in Bar is converted in Barye

def num_density(n0,C): #return th density based on the mixing ratio of a gaz and the global atmosphere density
    return C * n0
    
def horzX(R,z): #get the slant path x coordinates based on the current height z and its upper atm layers
    x = [0] #empty array
    
    rad = z + R
    rad0 = rad[0]

    for r in rad[1:]:
        x.append(np.sqrt(np.square(r) - np.square(rad0)))
    
    return np.asarray(x)

 
#sum all the molecules by the densities
def mix(xsec,comp,mols):
    return np.sum(comp[mol]*xsec[mol] for mol in mols)

#rearrange the xsec array to have each col correspond to the correct temp    
def arrange_xsec(T,xsec,mols):
    
    res = {}
    
    for m in mols:
        
        dxsec = np.empty((len(xsec[T[0]]),len(T)))
    
        for i,t in enumerate(T):
            dxsec[:,i] = xsec[t][m]
    
        res[m] = dxsec
    
    return res
    

#validate temperature profile
def validateTemp(T,N):
    
    #check if we have a isothermal case
    if not isinstance(T,list):
        return True,[T]
    #check if temps profile matces atmosphere layers
    elif len(T) != N :
        print('Warning: Temperature profile does not match atmosphere layers. Will only use first temperature instead')
        return True,[T[0]]
    #return the unique temperatures for xsec retrieval
    else:
        return False,np.unique(T)
        
def slant(y,path):
    return 2 * np.trapz(y,x = path,axis=1)    

def ring(R,r,trans):
    return (2*np.pi * np.trapz(r*trans,x=r,axis=1))/R
