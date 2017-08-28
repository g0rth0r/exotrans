# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 16:59:53 2017
Main
@author: Seb.BN
"""
import time
import const as c
import numpy as np
import exomol as em
import broadening
import exofunc as exf


def exotrans(**kwargs):
    #KWARGS INPT PARAMETERS********************************************************
    start_time = time.clock()
    #DEFINE WAVELENGHT PARAMEERS AND RESOLUTION
    v_min = kwargs['v_min']  #wavenumber in cm^-1
    v_max = kwargs['v_max'] #wavenumber in cm^-1
    dv = kwargs['dv']
    
    #DEFINE ATMOSPHERIC COMPOSITION AND PRESSURE PROFILE
    mol_in = kwargs['mol_in']
    atm_in = kwargs['atm_in']
    strats = kwargs['strats'] #number of stratification of the atmosphere profile
    P_min = kwargs['P_min'] #Pressure in Bar, outer atmosphere
    P_max = kwargs['P_max'] #max pressure; fully opaque afterward
    
    #DEFINE PLANET PROPRETIES
    M = kwargs['M'] #in earth's mass units; correspond to Jupiter's
    R = kwargs['R'] #in earth radius unit
    
    #non iso thermal case. Temps must be given from inner-most to outermost atm 
    #layer
    T = kwargs['T']
    p_broad = kwargs['p_broad']
    lorz_width = kwargs['lorz_width']
    res = kwargs['res']
    #modifier that defines from which atmosphere layer we want to integrate.
    mod = kwargs['mod']
    #******************************************************************************
    
    #pack molecular data
    molname_list = [mol[0] for mol in mol_in] #list of molecule names
    broadname_list = [mol[0] for mol in atm_in]
    
    #rest of molecular data: mixing ratio, molecular mass (g), including H, He
    comp_data = [mol[1:] for mol in atm_in + mol_in] 
    
    atm_layers = np.arange(strats)
    
    #ONE-TIME CONSTANT CALCULATIONS
    g = c.G*M/np.square(R)
    u_atm = exf.getMeanW(comp_data)
    H = (c.kb * np.mean(T))/(u_atm * g)
    
    #define the pressure/height profile
    z_max = exf.pressure_reci(P_max,H,P_min)
    z = np.linspace(0,z_max, num=strats)
    P = exf.pressure(P_max,H,z)
    
    #calculation of the atmosphere's composition profile
    n0 = exf.density(P,np.array(T))
    
    #generate coposition dictionary 
    comp = {}
    for mol in mol_in:
        comp[mol[0]] = mol[1] * n0
    
    broad_mix = {}
    for mol in atm_in:
        broad_mix[mol[0]] = mol[1]
    
    #validate temperature profile; i.e. isothermal or not
    iso,uT = exf.validateTemp(T,strats)
    
    #retreaval of multiple cross secion data. Set val to False to ovveride exomol 
    #parameters validation
    xsec_list = em.getMultExoData(dv,v_min,v_max,uT,*molname_list,val=False)
    
    #rearrange the xsec in a dictionary of 2d arrays of xsec/atm layers. Molecules
    # are dict keys.
    ar_xsec = exf.arrange_xsec(T,xsec_list,molname_list)
    
    xsec = xsec_list[T[0]]['wavenum']
 
    #if wanted, apply pressure broadening to the xsec dictionary
    if p_broad:
        #gerate all the P-T points. WIll be used for broadening profiles.
        PT = list(zip(P,T))
    
        #get the list of available broadeners
        avail = broadening.get_all_data(molname_list,broadname_list)
        
        #broaden the xsec data.Params in order: 
        #ar_xsec: ditionary of 2D array of xsec data
        #avail: reference list of available boradeners
        #broad_mix: mixing ratio data of atmosphere broadeners
        #PT: PT profile of each atm layers
        #lorz_width: how many times should he Lorentz profle be wide in * HWHM
        #res: resolution of the Lorz. profle.
        ar_xsec = broadening.broad(ar_xsec,avail,broad_mix,PT,lorz_width,res)
        
    #calculate the slant paths for each layer of the asmosphere
    paths = []
    for i in atm_layers:
        paths.append(exf.horzX(R,z[i:]))
    
    #get the mixing+xsec data
    mix = exf.mix(ar_xsec,comp,molname_list)
    
    #define the ring area that is incident to the light flux
    ring = (np.pi*(((z[-1:]+R) ** 2) - ((z[0]+R) ** 2)))
        
    #trs = np.empty(xlen) #empty array to store calculated transmision coefs.
    xtrans = np.empty_like(mix)
    
    #itterate over each atmosphere layers and integrate over the slant paths
    #the the full mix array is given and the result is appeneded as a column
    #each itteration truncate the mix matrix by one layer; 
    for i in range(strats-mod):
        xtrans[:,i+mod] = exf.slant(mix[:,i+mod:],paths[i+mod])
    
    #get the transmission coef. based on the optical depth calculate from the
    #Beer-Lambert Law.        
    flux = np.exp(-1 * xtrans)
    
    #integrate on the ring surface that correspond to the exoplanet atmosphere
    #linb and get the % transmitted.
    trs = exf.ring(ring,z[mod:]+R,flux[:,mod:])
    
    #to return only if wanted
    #end_time = time.clock() - start_time
    
    return xsec,trs
    

if __name__ == "__main__":
    
    #example of species definition: (mame, mixing ratio, molecular mass (g))
    s1 = ('H2',0.85,2*1.6737236e-24)
    s2 = ('He',0.15,6.6464764e-24)
    s3 = ('12C-1H4',2e-4,2.66391311e-23)
    s4 = ('12C-16O',2e-4,4.65118646e-23)
    s5 = ('1H2-16O',5e-4,2.99150758e-23)
    
    #example of dictionary that needs to be passed to exotrans
    test1 = {'v_min' : 4255,
        'v_max' : 4366,
        'dv' : 0.1,
        'mol_in' : [s3,s4,s5],
        'atm_in' : [s1,s2],
        'strats' : 20,
        'P_min' : 1e-5,
        'P_max' : 5,
        'M' : 317.83 * c.mer,
        'R' : 11.209 * c.rer,
        'T' : [1250] * 20,
        'mod' : 0,
        'p_broad' : True,
        'lorz_width' : 5,
        'res' : 13}
    
    #cp1 = cProfile.Profile()
    #cp1.enable()
        
    x,y = exotrans(**test1)

