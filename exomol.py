# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 17:22:23 2017
Methods of retrieving exomol cross sections
@author: Seb.BN
"""
import requests as rq
from bs4 import BeautifulSoup as bs
import re
import os.path
import numpy as np

#setup the working directory
workdir = 'xsec_data'

#custom validation error
class ValidationError(Exception):
    pass

class OutOfMethodsError(Exception):
    pass


def checkdir():
    if not os.path.isdir(workdir):
        print('Directory not found. Creating ./{}'.format(workdir))
        os.makedirs(workdir)

#ACCEPTED PARAMETERS DEFINIIONS
#all v in cm-1, T in K
params = { '12C-16O': {'dv_min' : 0.01,
                       'dv_max' : 100,
                       'min'   : 0,
                       'max'   : 8499,
                       'Tmin'   : 296,
                       'Tmax'   : 5000
                       },
           '12C-1H4': {'dv_min' : 0.01,
                       'dv_max' : 100,
                       'min'   : 0,
                       'max'   : 12000,
                       'Tmin'   : 296,
                       'Tmax'   : 2000
                       },
           '1H2-16O': {'dv_min' : 0.01,
                       'dv_max' : 100,
                       'min'   : 0,
                       'max'   : 30000,
                       'Tmin'   : 296,
                       'Tmax'   : 3000
                       }
          }


#check if gien parameters are within the permited range from Exomol
def params_check(dv, vmin, vmax,temp, molec):
    if (
        vmin < vmax and
        dv >= params[molec]['dv_min'] and 
        dv <= params[molec]['dv_max'] and
        vmin >= params[molec]['min'] and
        vmax <= params[molec]['max'] and
        temp >= params[molec]['Tmin'] and
        temp <= params[molec]['Tmax']
        ):
            return True
    else:
        return False

#return the expected file name according to the exomol naming convention
def file_name(dv, vmin, vmax, temp, molec):
    return '{}_{}-{}_{}K_{:.6f}.sigma'.format(molec,vmin,vmax,temp,dv)

def method1(dv, vmin, vmax, temp, molec): #if file is alrady downloaded and file name is known
    file = file_name(dv, vmin, vmax, temp, molec)
    if os.path.isfile('./{}/{}'.format(workdir,file)):
        #print('File Found!')
        data = np.genfromtxt('./{}/{}'.format(workdir,file), names=['wavenum',molec],deletechars='')
        return data
    else:
        print('File Not Found.')
        raise FileNotFoundError



#This method is a much more brute approach to getting the xsec values from exomol
#as it goes through their HTTP server, thus coul be flagged as spam is abused.
#Considering it opens a new session at every calls, I would proceed with caution
#when calling over a large array of temperatures.
#Input goes as such: dnu= spacing of the wavenumbers, numin= minium wavenumber
#numax = maximum wave number, temp=temperature (296-2000K), molname=molecule name
#all wavenumbers mut be in cm^-1
def method2(dnu,numin,numax,temp, molname): #retrieve data over web using a POST method
    
    q_url = 'http://exomol.com/xsec/' + molname + '/'
    base_url = 'http://exomol.com/'
    
    payload = {'dnu' : dnu, 'numin':numin, 'numax':numax,'T':temp, 'spoon_feed':'on'}
    
    with rq.Session() as s:
        r = s.get(q_url)
        
        if r.status_code == 200:
            
            token = re.findall(r'csrftoken=(.*?);',r.headers['Set-Cookie'])
            payload['csrfmiddlewaretoken'] = token[0]
            r = s.post(q_url,data=payload)
            if r.status_code == 200:
                soup = bs(r.content,"html.parser")
                data = soup.find_all('div',attrs={'class':'well'})
                for div in data:
                    links = div.findAll('a',href=re.compile('.sigma'))
                    for a in links:
                        sigma = a['href']
                        file = a.contents[0]
                    r = rq.get(base_url + sigma)
                    if r.status_code == 200:
                        with open('./{}/{}'.format(workdir,file),"wb") as f:
                            f.write(r.content)
                            data = np.genfromtxt('./{}/{}'.format(workdir,file), names=['wavenum',molname],deletechars='')
                            #print("Success!")
                        return data
                    else:
                        print('Error retrieving data')
                        r.raise_for_status()
                        
            else:
                print('Error retrieving data')
                r.raise_for_status()
        else:
            print('Error retrieving data')
            r.raise_for_status()
    
    #print(r.status_code)

def method3(): #use Table Access Protocol (TAP) provided by Exomol
    pass

def getExoData(*args,**kwargs): #MAIN FUNCTION
    
    #check and/or create the working directory of xsec data
    checkdir()

    #first check the parameters are valid
    if params_check(*args) or not kwargs['val']:
        
        #if it pass then first check if the file was previouslydownloaded
        try:
            data = method1(*args)
        except FileNotFoundError:
            print('No existing cross-section file found. Will attempt to download.')
            #try method2 by downloading from Exomol if method1 failed
            try:
                data = method2(*args)
            except rq.exceptions.RequestException as e:
                print(e)
                raise OutOfMethodsError('Out of methods to retrieve cross-section data.')
        
        return data
        
    else:
        print(*args)
        raise ValidationError('The given parameters are outside the supported values of Exomol.')

def getMultExoData(dv,v_min,v_max,temps,*args,**kwargs):#return multiple xsec values
        
    #define dict for temperatures
    d_data = {}
    for T in temps:
        #define parameters
        params = [dv,v_min,v_max,T]
        #define molecules list
        mols = args
        data = getExoData(dv,v_min,v_max,T,(mols[0]),val=kwargs['val'])
        
        if len(mols) == 1:
            d_data[T] = data
            
        else:
            
            #generate dtype for final data set
            d = [('wavenum', '<f8')]
            for i in mols:
                d.append((i,'<f8'))
            
            #create empty array to to join all data into and assign initial data
            full_data = np.empty_like(data, dtype = d)
            full_data['wavenum'] = data['wavenum']
            full_data[mols[0]] = data[mols[0]]
    
            for mol in mols[1::]:#itterate over remaining molecules
                data_tmp = getExoData(*params,mol,val=kwargs['val'])
                full_data[mol] = data_tmp[mol]
                
            d_data[T] = full_data

    return d_data


