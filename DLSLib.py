#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 09:15:04 2022

@author: Christian Beck

christianbeck91@gmx.de

ORCID: https://orcid.org/0000-0001-7214-3447

"""
import numpy as np
import time
from lmfit import Model, Parameters
import pickle
import warnings
import h5py
import scipy as sp
warnings.filterwarnings("ignore")

#%% different functions
# save to pickle file
def save(FileName,Data):
    with open(FileName, "wb") as file:
        pickle.dump(Data, file, pickle.HIGHEST_PROTOCOL)
# load data
def load(Filename):
    with open(Filename, "rb") as file:
        loaded_dict = pickle.load(file)
    return loaded_dict
#%% save to h5py
def savehdf(name,dataset):
    f=h5py.File(name,'w')
    path='/'
    if isinstance(dataset,dict):
        savehdfdict(dataset,path,f)
    elif isinstance(dataset,list):
        savehdflist(dataset,path,f)
    f.close()
    
def savehdflist(liste,path,f):
    for hi,listentry in enumerate(liste):
        if isinstance(listentry,dict):
            savehdfdict(listentry, path + 'Dataset_' + str(hi) +'/', f)
        elif isinstance(listentry,list):
            savehdflist(listentry, path + 'Dataset_' + str(hi) +'/', f)
        else:
            # print(path + 'entry' + str(hi) +'/')
            # print('liste: ' + path)
            # print(hi)
            # print(listentry)
            # try:
            f.create_dataset(name=path + 'entry' + str(hi),data=listentry)
            # except:
            #     f.create_dataset(name='/Dataset_0/RedData/Dataset_0/angle/',data=np.array(30))
            #     print('hardcoded')
                    

def savehdfdict(dictionary,path,f):
    #print('dict: ' + path)
    for key in dictionary.keys():
        if isinstance(dictionary[key],dict):
            savehdfdict(dictionary[key],path + key + '/',f)
        elif isinstance(dictionary[key],list):
            savehdflist(dictionary[key],path + key + '/',f)
        else:
            if dictionary[key] is None:
                dictionary[key]=np.nan
            f.create_dataset(path + key,data=dictionary[key])
def loadhdf(Filename):
    with h5py.File(Filename, "r") as f:
        a_group_key = list(f.keys())[0]
        print('loading of hdf5 file is not yet implemented')
        
#%% alternative way of implementation
# def savehdf(name,dataset):
#     f=h5py.File(name,'w')
#     if isinstance(dataset,dict):
#         savehdfdict(dataset,f)
#     elif isinstance(dataset,list):
#         savehdflist(dataset,f)
#     f.close()
# def savehdflist(liste,group):
#     for hi,listentry in enumerate(liste):
#         if isinstance(listentry,dict):
#             grp=group.create_group('Dataset_' + str(hi))
#             savehdfdict(listentry,grp)
#         elif isinstance(listentry,list):
#             grp=group.create_group('Dataset_' + str(hi))
#             savehdflist(listentry, grp)
#         else:
#             print('')
#             # print(liste)
#             # try:
#             #     group['value']=liste
#             # except:
#             #     print(group.name +' already exists...')
            
# def savehdfdict(dictionary,group):
#     #print('dict: ' +  group.name)
#     for key in dictionary.keys():
#         if isinstance(dictionary[key],dict):
#             grp=group.create_group(key)
#             savehdfdict(dictionary[key],grp)
#         elif isinstance(dictionary[key],list):
#             grp=group.create_group(key)
#             savehdflist(dictionary[key],grp)
#         else:
#             if dictionary[key] is None:
#                 dictionary[key]=np.nan
#             # print('==========')
#             # print(dictionary[key])
#             # print(key)
#             # print(group.name)
#             # group[key]=dictionary[key]
        
#%%
# readin data
def readin(file):
    f=open(file,'r',encoding="ISO-8859-1")
    FC={}
    dat=0
    for line in f:
        if dat==0 and line.__contains__(":"):
            idx=line.find(":")
            FC.update({line[0:idx-1]:line[idx+1:-1]})
        elif not line.__contains__(":") and line.__contains__('"'):
            dat=1
            dattype=line.replace('"','')
            DAT=np.array([])
        elif dat==1:
            if len(DAT)==0:
                try:
                    DAT=[np.array(line.split(),dtype=float)]
                except:
                    time.sleep(0)
            elif 'Monitor' in line:
                FC.update({'monitor':float(line[15:])})
            else:
                try:
                    DAT=np.append(DAT,[np.array(line.split(),dtype=float)], axis=0)
                except:
#                    print('##############')
#                    if 'Monitor' in line:
#                        FC.update({'monitor':float(line[15:])})
#                        print('Monitor added' + str(FC['monitor']))
#                    else:
                    time.sleep(0)
                FC.update({dattype:DAT})
    f.close()
    data={"angle": float(FC['Angle [°]      ']),
          "intensities":FC['Count Rate\n'],
          "g2m1":FC['Correlation\n'][:,1:3],
          "tau":FC['Correlation\n'][:,0],#ms
          "T":float(FC['Temperature [K]']),
          "filename":file,
          "Date":str(FC['Date']).replace("\t", ''),
          "Time":str(FC['Time']).replace("\t", ''),
          "Monitor":FC['monitor']}
    wl=632.8*1e-9
    data.update({'q':np.array(4*np.pi*1.332/wl*np.sin(data['angle']/360*np.pi))})#nm^-1
    return data
def writesummary(filename,line,perm='a'):
    with open(filename,perm) as f:
        f.write(line)
        f.write('\n')
#%% define fit functions
def gaussdistexpnorm(x,w1,sig):
    y= (0.5*np.exp(0.5*x*(sig**2*x-2*w1))*(sp.special.erf((w1-sig**2*x)/(2*sig**2)**0.5)+1))**2
    y[np.isnan(y)]=0
    y=y/y[0]
    return y
def normalandstretchednorm(x,a1,w1,w2,b):
    return a1*np.exp(-2*w1*x)+(1-a1)*np.exp(-(2*w2*x)**b)
def stretchedstretchednorm(x,a1,w1,w2,b1,b2):
    return a1*np.exp(-(2*w1*x)**b1)+(1-a1)*np.exp(-(2*w2*x)**b2)
def doubleexponential(x,a1,a2,w1,w2):
    return a1*np.exp(-2*w1*x)+a2*np.exp(-2*w2*x)
def singleexponentialnorm(x,w1):
    return np.exp(-2*w1*x)
def doubleexponentialnorm(x,a1,w1,w2):
    return a1*np.exp(-2*w1*x)+(1-a1)*np.exp(-2*w2*x)
def Dq2(x,D):
    return x*D
def returnfitmodel(name):
    if name=='doubleexponential':
        model = Model(doubleexponential)
        params = Parameters()
        params.add('a1', value=0.50, min=0, max=1)
        params.add('a2', value=0.50, min=0, max=1)
        params.add('w1', value=0.50, min=0)
        params.add('delta', value=0.50, min=0)
        params.add('w2', expr='w1+delta')
    elif name=='gaussdistexpnorm':
        model = Model(gaussdistexpnorm)
        params = Parameters()
        params.add('w1', value=10, min=0)
        params.add('sig', value=1, min=0)
    elif name=='singleexponentialnorm':
        model = Model(singleexponentialnorm)
        params = Parameters()
        params.add('w1', value=0.50, min=0)
    elif name=='doubleexponentialnorm':
        model = Model(doubleexponentialnorm)
        params = Parameters()
        params.add('a1', value=0.50, min=0, max=1)
        params.add('w1', value=0.50, min=0)
        params.add('delta', value=0.50, min=0)
        params.add('w2', expr='w1+delta')
    elif name=='normalandstretchednorm':
        model = Model(normalandstretchednorm)
        params = Parameters()
        params.add('a1', value=0.50, min=0, max=1)
        params.add('w1', value=0.50, min=0)
        params.add('delta', value=0.50, min=0)
        params.add('w2', expr='w1+delta')
        params.add('b',  value=1, min=0)
    elif name=='stretchedstretchednorm':
        model = Model(stretchedstretchednorm)
        params = Parameters()
        params.add('a1', value=0.50, min=0, max=1)
        params.add('w1', value=0.50, min=0)
        params.add('delta', value=0.50, min=0)
        params.add('w2', expr='w1+delta')
        params.add('b1',  value=1, min=0)
        params.add('b2',  value=1, min=0)
    elif name=='Dq2':
        model=Model(Dq2)
        params = Parameters()
        params.add('D', value=0.50, min=0)
    return model,params
