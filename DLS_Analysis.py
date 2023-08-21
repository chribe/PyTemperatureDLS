#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version 15/12/2023 14:00

@author: Christian Beck

christianbeck91@gmx.de

ORCID: https://orcid.org/0000-0001-7214-3447

"""
import DLSLib as DL
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # <--- This is important for 3d plotting 
import numpy as np
from tqdm import tqdm
import time
import ilt
import warnings
warnings.filterwarnings("ignore")
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
#%% some constants
dT=1.25
DLS=[]
rdpath='./rawdata/'
dirfiles = os.listdir(rdpath)
#fullpaths = map(lambda name: os.path.join(dirname, name), dirfiles)
dirs = []
for folder in dirfiles:
    if os.path.isdir(rdpath + folder): dirs.append(folder)
print('\n######################\n########Files#########\n######################\n')
print(dirs)
time.sleep(0.1)
print('\n######################\n#######Loading########\n######################\n')
for i in tqdm (range(len(dirs)),desc="Loading…",ascii=False, ncols=75):
    Rawdata=dirs[i]
    dlsdata=[]
    files = os.listdir(rdpath+Rawdata)
    files.sort()
    for F in files:
        if '.ASC' in F:
            dlsdata.append(DL.readin(rdpath+Rawdata + '/' + F))
    DLS.append({'Name':Rawdata,
                'Data':dlsdata})
#%% merging different runs
colors = plt.get_cmap('jet', 13)
inttresh=5
for DLSmeas in DLS:
    currT=DLSmeas['Data'][0]['T']
    currangle=DLSmeas['Data'][0]['angle']
    currdata=[]
    reddata={}
    DLSmeas.update({'RedData':[]})
    for hi,meas in enumerate(DLSmeas['Data']):
        if meas['angle'] in np.arange(30,160,10):    
            if (((meas['T']-currT)**2<dT**2 and currangle==meas['angle']) or len(currdata)==0) \
                and hi!=len(DLSmeas['Data'])-1:
                currdata.append(meas)
            else:
                for entr in currdata[0].keys():
                    if entr not in ['g2m1','intensities','tau']:
                        reddata.update({entr:[]})
                        for curd in currdata:
                            reddata[entr].append(curd[entr])
                reddata.update({'tau':currdata[0]['tau']})
                discarted=[]
                for hi1,curd in enumerate(currdata):
                    if hi1==0:
                        g2m1=curd['g2m1']
                        intensity=curd['intensities'][:,1:3]
                    else:
                        g2m1=np.append(g2m1,curd['g2m1'],axis=1)
                        if len(intensity[:,0])==len(curd['intensities'][:,0]):
                            intensity=np.append(intensity,curd['intensities'][:,1:3],axis=1)
                        else:
                            if len(intensity[:,0])<len(curd['intensities'][:,0]):
                                intensity=np.append(intensity,curd['intensities'][0:len(curd['intensities'][:,0])-1,1:3],axis=1)
                            else:
                                intensity=np.append(intensity,intensity[:,0:2]*10**6,axis=1)
                                intensity[0:len(curd['intensities'][:,0]),-2:]=curd['intensities'][:,1:3]
                for hi2 in np.arange(len(g2m1[0,:])-1,-1,-1):
                    if (intensity[:,hi2]<inttresh).any():
                        np.delete(g2m1,np.s_[hi2],axis=1)
                        discarted.append(hi2)
                for hi2 in np.arange(len(g2m1[0,:])):
                    g2m1[:,hi2]=g2m1[:,hi2]/g2m1[reddata['tau']<0.001,hi2].mean()
                reddata.update({'g2m1':g2m1.mean(axis=1),
                                'dg2m1':g2m1.std(axis=1),
                                'discarted':np.array(discarted)})
                DLSmeas['RedData'].append(reddata)
                reddata={}
                currdata=[]
                currdata.append(meas)
                currT=meas['T']
                currangle=meas['angle']
#%% create files for summaries
if not os.path.exists('Summaries'):
    os.mkdir('Summaries')
for DLSmeas in DLS:
    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\documentclass{article}\usepackage{graphicx}\usepackage{hyperref}\begin{document}',perm='w')
    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\author{Automatically Compiled Summary}\title{' + DLSmeas['Name'].replace('_',' ') +r'}\maketitle\newpage\tableofcontents')

#%% Contin algorithm
# based on https://github.com/caizkun/pyilt
bound = np.array([0.0001, 100])
alpha = 1
print('\n######################\n#######Contin#########\n######################\n')
for DLSmeas in DLS:
    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\newpage')
    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\section{Contin Analysis}')
    for j in tqdm (range(len(DLSmeas['RedData'])),desc="Contin " + DLSmeas['Name'] +': ',ascii=False, ncols=100):
        meas=DLSmeas['RedData'][j]
        t=meas['tau']
        F=meas['g2m1']
        z_kww, f_kww, res_lsq, res_reg = ilt.ilt(t, F, bound, len(F), alpha)
        Contin={'t':1./z_kww,
                'distrib':z_kww*f_kww}
        meas.update({'Contin':Contin})
#%% determine Temperatures
Tlimit=[]
for DLSmeas in DLS:
    for meas in DLSmeas['RedData']:
        Tlimit.append((np.array(meas['T'])*2).mean().round()/2)
Tlimit=sorted(list(set(Tlimit)))
Tdifference=np.array(Tlimit)[1:] -np.array(Tlimit)[0:-1]
Tlimits=[]
hi1=0
# merge different temperatures
for hi,Td in enumerate(Tdifference):
    #print(Td)
    if Td>dT:
        #print(Tlimit[hi1:hi+1])
        Tlimits.append(np.mean(Tlimit[hi1:hi+1]))
        hi1=hi+1
Tlimits.append(np.mean(Tlimit[hi1:]))        
#%% plot Contin
#Tlimits=np.arange(285,319,2.5)
if not os.path.exists('figures'):
    os.mkdir('figures')
if not os.path.exists('figures/contin'):
    os.mkdir('figures/contin')
containscurve=False
for Tl in Tlimits:
    for DLSmeas in DLS:
        fig=plt.figure()
        ax = fig.add_subplot(projection='3d')
        for meas in DLSmeas['RedData']:
            if np.abs(np.mean(meas['T'])-Tl)<dT:
                ax.scatter(np.ones(np.shape(meas['Contin']['t']))*np.mean(meas['q'])**2,\
                           np.log10(meas['Contin']['t']),meas['Contin']['distrib'],marker='o',\
                               color=colors(int(np.arange(0,13)[np.arange(30,160,10)==np.mean(meas['angle'])])))
            if not np.shape([])==np.shape(meas['Contin']['t']):
                containscurve=True
        ax.set_title(str(Tl))
        ax.set_xlabel(r'$q^2$ [nm$^{-2}$]')
        ax.set_ylabel(r'log$_{10}\left(t\right)$')
        ax.set_zlim([0,1])
        if containscurve:
            fig.savefig('figures/contin/' +DLSmeas['Name'] +str(Tl) +'.pdf',bbox_inches='tight')
            DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex', r'\includegraphics[width=.4\textwidth]{'+'../figures/contin/' +DLSmeas['Name'] +str(Tl) +'.pdf'+'} ')
        plt.close()
#%%
print('\n######################\n########Fits##########\n######################\n')
Fittypes=['doubleexponential','doubleexponentialnorm','normalandstretchednorm','stretchedstretchednorm']
for fit in Fittypes:
    if not os.path.exists('figures/'+fit):
        os.mkdir('figures/'+fit)
for DLSmeas in DLS:
    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\newpage')
    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\section{Fits}')
    for j in tqdm (range(len(DLSmeas['RedData'])),desc="Fit " + DLSmeas['Name'] +': ',ascii=False, ncols=100):
        meas=DLSmeas['RedData'][j]
        meas.update({'Fit':[]})
        for fit in Fittypes:
            plt.errorbar(meas['tau'], meas['g2m1'],meas['dg2m1'],fmt='o',\
                         color=colors(int(np.arange(0,13)[np.arange(30,160,10)==np.mean(meas['angle'])])),alpha=0.1)
            plt.plot(meas['tau'], meas['g2m1'],'o',\
                         color=colors(int(np.arange(0,13)[np.arange(30,160,10)==np.mean(meas['angle'])])))
            model,param=DL.returnfitmodel(fit)
            param.add('q',value=np.array(meas['q']).mean(),vary=False)
            FR=model.fit(meas['g2m1'],param,x=meas['tau'],weights=1/meas['dg2m1'])
            Fit={'Name':fit}
            for p in FR.best_values.keys():
                Fit.update({'d'+p : FR.params[p].stderr})
                Fit.update({p : FR.params[p].value})
            meas['Fit'].append(Fit)
            x=np.logspace(np.log10(min(meas['tau'])),np.log10(max(meas['tau'])),100)
            plt.plot(x,FR.eval(x=x),color='k',zorder=10)#color=colors(int(np.arange(0,13)[np.arange(30,160,10)==np.mean(meas['angle'])]))
            name=DLSmeas['Name']+ ' ' +str(np.mean(meas['q']))+' ' +str(np.mean(meas['T']))
            plt.xscale('log')
            plt.xlabel(r'$\tau$ [ms]')
            plt.ylabel(r'g$_2$-1')
            plt.title(r'$q$='+str( format(np.array(meas['q']).mean(), '.2f'))+ r'nm$^{-1}$')
            plt.savefig('figures/'+ fit +'/' +  name +'.pdf',bbox_inches='tight')
            plt.close()
    for fit in Fittypes:
        DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\subsection{'+ fit +'}')
        for Tl in Tlimits:
            DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\subsubsection{'+ fit +':'+ str(Tl)+'K}')
            for meas in DLSmeas['RedData']:
                if np.abs(np.mean(meas['T'])-Tl)<1.25:
                    name=DLSmeas['Name']+ ' ' +str(np.mean(meas['q']))+' ' +str(np.mean(meas['T']))
                    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\includegraphics[width=.4\textwidth]{../figures/'+ fit +'/' +  name +'.pdf}')

#%% determine D from fit of Gamma vs q^2
for DLSmeas in DLS:
    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\section{Fit post processing}')
    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\subsection{$q$-dependence of fit parameter}')
    DLSmeas.update({'Fitanalysis':{}})
    for hi,fit in enumerate(Fittypes):
        DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\subsubsection{'+ fit +'}')
        if not os.path.exists('figures/fitanalysis_'+fit):
            os.mkdir('figures/fitanalysis_'+fit)
        w1=[]
        dw1=[]
        w2=[]
        dw2=[]
        T=[]
        q=[]
        name=[]
        DLSmeas['Fitanalysis'].update({fit:{'T':np.array([]),'D1':np.array([]),'dD1':np.array([]),'D2':np.array([]),'dD2':np.array([])}})
        for meas in DLSmeas['RedData']:
            w1.append(meas['Fit'][hi]['w1'])
            dw1.append(meas['Fit'][hi]['dw1'])
            w2.append(meas['Fit'][hi]['w2'])
            dw2.append(meas['Fit'][hi]['dw2'])
            name.append(meas['Fit'][hi]['Name'])
            q.append(np.mean(meas['q']))
            T.append(np.mean(meas['T']))
        dw1=np.array(dw1,dtype=float)
        dw2=np.array(dw2,dtype=float)
        for Tl in Tlimits:
            if Tl<max(np.array(T)+3):
                if not np.shape([])==np.shape(np.array(w1)[(T<Tl+1.25)&(T>Tl-1.25)]):
                    DLSmeas['Fitanalysis'][fit]['T']=np.append(DLSmeas['Fitanalysis'][fit]['T'],Tl)
                    fig, ax1 = plt.subplots()
                    ax1.errorbar(np.array(q)[(T<Tl+1.25)&(T>Tl-1.25)]**2,np.array(w1)[(T<Tl+1.25)&(T>Tl-1.25)],dw1[(T<Tl+1.25)&(T>Tl-1.25)],fmt='o')
                    model,param=DL.returnfitmodel('Dq2')
                    dw1[dw1==None]=1000
                    dw1=np.array(dw1,dtype=float)
                    dw1[np.isnan(dw1)]=1000
                    g1=np.array(w1)[(T<Tl+1.25)&(T>Tl-1.25)]
                    dg1=1/dw1[(T<Tl+1.25)&(T>Tl-1.25)]
                    q2=np.array(q)[(T<Tl+1.25)&(T>Tl-1.25)]**2
                    Dq2fit1=model.fit(g1[q2>250],param,x=q2[q2>250],weights=dg1[q2>250])
                    x=np.linspace(0,650)
                    ax1.plot(x,Dq2fit1.eval(x=x),color='k',zorder=10)
                    ax1.set_xlabel(r'$q^2$ [nm$^{-2}$]')
                    ax1.set_ylabel(r'$\gamma$ [ms$^{-1}$]')
                    ax2 = ax1.twinx()
                    color = 'tab:green'
                    ax2.errorbar(np.array(q)[(T<Tl+1.25)&(T>Tl-1.25)]**2,np.array(w2)[(T<Tl+1.25)&(T>Tl-1.25)],dw2[(T<Tl+1.25)&(T>Tl-1.25)],fmt='s',color = color)
                    dw2[dw2==None]=1000
                    dw2[dw2==np.inf]=1000                
                    dw2=np.array(dw2,dtype=float)
                    dw2[np.isnan(dw2)]=1000
                    g2=np.array(w2)[(T<Tl+1.25)&(T>Tl-1.25)]
                    dg2=1/dw2[(T<Tl+1.25)&(T>Tl-1.25)]
                    model2,param2=DL.returnfitmodel('Dq2')
                    Dq2fit2=model2.fit(g2[q2>250],param2,x=q2[q2>250],weights=dg2[q2>250])
                    ax2.plot(x,Dq2fit2.eval(x=x),color=color,zorder=10)
                    DLSmeas['Fitanalysis'][fit]['D1']=np.append(DLSmeas['Fitanalysis'][fit]['D1'],Dq2fit1.params['D'].value)
                    DLSmeas['Fitanalysis'][fit]['dD1']=np.append(DLSmeas['Fitanalysis'][fit]['dD1'],Dq2fit1.params['D'].stderr)
                    DLSmeas['Fitanalysis'][fit]['D2']=np.append(DLSmeas['Fitanalysis'][fit]['D2'],Dq2fit2.params['D'].value)
                    DLSmeas['Fitanalysis'][fit]['dD2']=np.append(DLSmeas['Fitanalysis'][fit]['dD2'],Dq2fit2.params['D'].stderr)
                    ax2.set_ylabel(r'$\gamma_2$ [ms$^{-1}$]', color = color)
                    ax2.tick_params(axis ='y', labelcolor = color)
                    ax2.set_ylim([0,1.1*max(np.array(w2)[(T<Tl+1.25)&(T>Tl-1.25)])])
                    ax1.set_ylim([0,1.1*max(np.array(w1)[(T<Tl+1.25)&(T>Tl-1.25)])])
                    ax1.set_title(str(Tl))
                    fig.savefig('figures/fitanalysis_'+fit +'/' +DLSmeas['Name']+ str(Tl) +'.pdf',bbox_inches='tight')
                    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\includegraphics[width=.4\textwidth]{../figures/fitanalysis_'+fit +'/' +DLSmeas['Name']+ str(Tl) +'.pdf}')
                    plt.close()
#%% plot fit results as a function of T
for DLSmeas in DLS:
    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\section{Temperature dependence}')
    for fittype in DLSmeas['Fitanalysis'].keys():
        plt.figure()
        plt.errorbar(DLSmeas['Fitanalysis'][fittype]['T'],DLSmeas['Fitanalysis'][fittype]['D1'],DLSmeas['Fitanalysis'][fittype]['dD1'],fmt='s')
        plt.errorbar(DLSmeas['Fitanalysis'][fittype]['T'],DLSmeas['Fitanalysis'][fittype]['D2'],DLSmeas['Fitanalysis'][fittype]['dD2'],fmt='o')
        plt.title(fittype)
        plt.xlabel('T [K]')
        plt.ylabel(r'D [nm$^2$/ms]')
        plt.savefig('figures/fitanalysis_'+fittype +'/' +DLSmeas['Name'] +'T_dependence.pdf',bbox_inches='tight')
        plt.close()
        DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\includegraphics[width=.4\textwidth]{../figures/fitanalysis_'+fittype +'/' +DLSmeas['Name'] +'T_dependence.pdf}')
#%% compile Summaries
print('\n######################\n###compile summaries##\n######################\n')
for DLSmeas in DLS:
    DL.writesummary('Summaries/' + DLSmeas['Name'] + '.tex',r'\end{document}')
    os.chdir('Summaries')
    try:
        os.system('xelatex ' + DLSmeas['Name'] + '.tex')
        os.system('xelatex ' + DLSmeas['Name'] + '.tex')
    except:
        print('Compiling of ' + DLSmeas['Name'] + '.tex  failed!')
    os.chdir('..')
#%% save Data Structure
DL.save("DLS_analysis.pickle",DLS)
DL.savehdf("DLS_analysis.hdf5",DLS)
