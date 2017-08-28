# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 16:46:25 2017

@author: Seb.BN
"""

import numpy as np
import os
import requests as rq
import matplotlib.pyplot as plt


#define dictionary of broadeners. Need to be updated if we want to work with 
#more broadeners. Format will be: 
# [Molecule][Default Lorz half wdth,Default Temp exponent,[Broadener]:[.broad file name, URL]]
#The default values are given per the .def files of each molecule on Exomol.
#for CO, the mean values of the parameters given by Faure were used.

broad_param = {'12C-1H4': [0.0488,0.400,{'He':['12C-1H4__He.broad','http://exomol.com/db/CH4/12C-1H4/12C-1H4__He.broad'],
                           'H2':['12C-1H4__H2.broad','http://exomol.com/db/CH4/12C-1H4/12C-1H4__H2.broad']}],
               '1H2-16O': [0.0700,0.500,{'He':['1H2-16O__He.broad','http://exomol.com/db/H2O/1H2-16O/1H2-16O__He.broad'],
                           'H2':['1H2-16O__H2.broad','http://exomol.com/db/H2O/1H2-16O/1H2-16O__H2.broad']}],
               '12C-16O': [2.794,0.610,{'H2':['fit-co.tab','http://exomol.com/db/CO/12C-16O/Faure/fit-co.tab']}]
              }

#working directory for broadening related files
broad_dir = './broad/'

#check if the working directory exists
def checkdir():
    if not os.path.isdir(broad_dir):
        print('Directory not found. Creating {}'.format(broad_dir))
        os.makedirs(broad_dir)

def check_file(filename):
    if os.path.isfile('{}{}'.format(broad_dir,filename)):
        return True
    else:
        return False        

#if needed, download the .broad file from Exomol
def get_broad_file(url,filename):
    print('Downloading {}...'.format(filename))
    r = rq.get(url)
    if r.status_code == 200:
        with open('{}{}'.format(broad_dir,filename),"wb") as f:
            f.write(r.content)
    else:
        print('Error retrieving Broadening data')
        r.raise_for_status()

def get_all_data(mols,broads):
    checkdir()
    #will store the available data in an array as tuple (molecule, broadener)
    available = {}
    for mol in mols:
        available[mol] = []
        #check if we have data for the broadener
        if mol in broad_param:
            for b in broads:
                if b in broad_param[mol][2]:
                    if not check_file(broad_param[mol][2][b][0]):
                        get_broad_file(broad_param[mol][2][b][1],broad_param[mol][2][b][0])
                    available[mol].append((mol,b))
                else:
                    print('No data for broadener {} with {}. Skipping.'.format(b,mol))
        else:
            print('No broadening data for {} available. Skipping.'.format(mol))

    return available
        
#return the Lorentzian pressure broadening half-width at half-maximum for a 
#P-T point. Parameters are: 
#available: list of available boradener for each molecules
#PT: the current P-T point
#mixing: reference dictionary of mixing ratios
#usedefault (bool): use default gamme and n values
def lorz_param(PT,available,mixing,useDefault):
    #reference values from ExoMol
    T_ref = 296 #kelvin
    P_ref = 1 #bar
    n = []
    gamma_b = []
    
    for sp in available:
        n.append(broad_param[sp[0]][1])
        gamma_b.append(broad_param[sp[0]][0] * (mixing[sp[1]] * PT[0]))
        
    return ((T_ref/PT[1])**np.mean(n)) * (PT[0]/P_ref) * np.sum(gamma_b)

def L(x,x0, gamma): #Return Lorentzian line shape at x with HWHM gamma, centered on x0
    return gamma / (np.pi * ((x-x0)**2 + gamma**2))

#broadens the xsec data givne in a dictionary form where key = molecule,
#using the list of avaiable broadeners (avail), broadeners mixing ratio 
#broad_data) and the PT profile of the atmosphere
def broad(d_xsec,avail,broad_data,PT,lorz_width,res):
    #empty dict to be passed return at the end
    b_xsec = {}
    #iterate over all the molecules types
    for mol in d_xsec:
        
        gamma_l0 = lorz_param(PT[0],avail[mol],broad_data,True)
        x0 = np.linspace(0-lorz_width*gamma_l0,lorz_width*gamma_l0,res)
        
        lay = 0
        temp = np.empty_like(d_xsec[mol].T)
        #iter over all the atmosphere layers
        for col in d_xsec[mol].T:
            gamma_l = lorz_param(PT[lay],avail[mol],broad_data,True)
            x = np.linspace(-lorz_width*gamma_l,lorz_width*gamma_l,res)
            x2 = x0[(x0 >= -5*gamma_l) * (x0 <= 5*gamma_l)]
            L_prof = L(x,0, gamma_l)
            L_norm = L_prof/np.sum(L_prof)
            temp[lay] = np.convolve(col,L_norm,'same')
            lay +=1
            
        b_xsec[mol] = temp.T

    return b_xsec
            

    
    
    
    