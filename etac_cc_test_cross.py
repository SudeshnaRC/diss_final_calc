# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 18:13:20 2019

@author: srcha
"""

import collections 
import gvar as gv
import numpy as np
import corrfitter as cf
import lsqfit
import scipy
import matplotlib.pyplot as plt

def make_data(file):
    return cf.read_dataset(file)


data = make_data("relabelled-etac-fullwall-testset.gpl")

#set correlator with tag 'etac_t0+2' as main eta_c correlator
def make_etac():
    etac = data['etac_t0+2']
    return(etac)

#keep all other tags in matrix cc to make cc_loop
def make_cc():
    cc = np.zeros((10,820)) # potentially might throw up an error later
    for i in range(2,11):
        tag = "etac_t0+{0}".format(2*i)  
        for j in range(820):
            for k in range(64):
                cc[i-1][j]+=data[tag][j][k]
                
    return(cc)
 
etac = make_etac()
cc = make_cc()
    
#print(cc)

print(np.shape(etac))

def correlator():
    correlator = collections.OrderedDict()
    for i in range(2,11):
        correlator['etac'] = np.zeros((820,64))
        correlator['cc{0}'.format(i)] = np.zeros(820)
        correlator['etac_cc{0}'.format(i)] = np.zeros((820,64))
        for t in range(64):
            for j in range(820):
                correlator['etac_cc{0}'.format(i)][j][t] += etac[j][t]*cc[i-1][j]
        correlator['etac'] = etac
        correlator['cc{0}'.format(i)] = cc[i-1]
    AvCorr = gv.dataset.avg_data(correlator)
        #print(AvCorr)
    return(AvCorr)
AvCorr = correlator()

    
def plot_data_21():
    T_samp = 21
    T = np.zeros(9)
    mean = np.zeros(9)
    yerr = np.zeros(9)
    plt.figure(1)
    for i in range(2,11):
        T[i-2]=i
        mean[i-2] = -1*0.6382*(48**3)*((AvCorr['etac_cc{0}'.format(i)][T_samp]-AvCorr['etac'][T_samp]*AvCorr['cc{0}'.format(i)])/(AvCorr['etac'][T_samp])).mean
        yerr[i-2] = 0.6382*(48**3)*((AvCorr['etac_cc{0}'.format(i)][T_samp]-AvCorr['etac'][T_samp]*AvCorr['cc{0}'.format(i)])/(AvCorr['etac'][T_samp])).sdev
        
        #print(AvCorr['etac'][T_samp])
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T-2, mean, yerr=[yerr,yerr], capsize=2, fmt = '-ko')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    plt.show()
    
    return()
    
def plot_data_22():
    T_samp = 22
    T = np.zeros(9)
    mean = np.zeros(9)
    yerr = np.zeros(9)
    plt.figure(1)
    for i in range(2,11):
        T[i-2]=i
        mean[i-2] = -1*0.6382*(48**3)*((AvCorr['etac_cc{0}'.format(i)][T_samp]-AvCorr['etac'][T_samp]*AvCorr['cc{0}'.format(i)])/(AvCorr['etac'][T_samp])).mean
        yerr[i-2] = 0.6382*(48**3)*((AvCorr['etac_cc{0}'.format(i)][T_samp]-AvCorr['etac'][T_samp]*AvCorr['cc{0}'.format(i)])/(AvCorr['etac'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T-2, mean, yerr=[yerr,yerr], capsize=2, fmt = '-go')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    #plt.legend(['T=21','T=23','T=25'])
    plt.grid(True)
    plt.title('$<\eta_c|c\overline{c}|\eta_C>$')
    plt.show()
    
    return()
    
def plot_data_23():
    T_samp = 23
    T = np.zeros(9)
    mean = np.zeros(9)
    yerr = np.zeros(9)
    plt.figure(1)
    for i in range(2,11):
        T[i-2]=i
        mean[i-2] = -1*0.6382*(48**3)*((AvCorr['etac_cc{0}'.format(i)][T_samp]-AvCorr['etac'][T_samp]*AvCorr['cc{0}'.format(i)])/(AvCorr['etac'][T_samp])).mean
        yerr[i-2] = 0.6382*(48**3)*((AvCorr['etac_cc{0}'.format(i)][T_samp]-AvCorr['etac'][T_samp]*AvCorr['cc{0}'.format(i)])/(AvCorr['etac'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T-2, mean, yerr=[yerr,yerr], capsize=2, fmt = '-ro')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    plt.show()
    
    return()
    
def plot_data_24():
    T_samp = 24
    T = np.zeros(9)
    mean = np.zeros(9)
    yerr = np.zeros(9)
    plt.figure(1)
    for i in range(2,11):
        T[i-2]=i
        mean[i-2] = -1*0.6382*(48**3)*((AvCorr['etac_cc{0}'.format(i)][T_samp]-AvCorr['etac'][T_samp]*AvCorr['cc{0}'.format(i)])/(AvCorr['etac'][T_samp])).mean
        yerr[i-2] = 0.6382*(48**3)*((AvCorr['etac_cc{0}'.format(i)][T_samp]-AvCorr['etac'][T_samp]*AvCorr['cc{0}'.format(i)])/(AvCorr['etac'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T-2, mean, yerr=[yerr,yerr], capsize=2, fmt = '-mo')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    #plt.legend(['T=21','T=23','T=25'])
    plt.grid(True)
    plt.title('$<\eta_c|c\overline{c}|\eta_C>$')
    plt.show()
    
    return()
    
def plot_data_25():
    T_samp = 25
    T = np.zeros(9)
    mean = np.zeros(9)
    yerr = np.zeros(9)
    plt.figure(1)
    for i in range(2,11):
        T[i-2]=i
        mean[i-2] = -1*0.6382*(48**3)*((AvCorr['etac_cc{0}'.format(i)][T_samp]-AvCorr['etac'][T_samp]*AvCorr['cc{0}'.format(i)])/(AvCorr['etac'][T_samp])).mean
        yerr[i-2] = 0.6382*(48**3)*((AvCorr['etac_cc{0}'.format(i)][T_samp]-AvCorr['etac'][T_samp]*AvCorr['cc{0}'.format(i)])/(AvCorr['etac'][T_samp])).sdev
        
        print(yerr)
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T-2, mean, yerr=[yerr,yerr], capsize=2, fmt = '-bo')
    plt.xlabel('t\'(a)')
    plt.ylabel('$C_{3pt}/C_{2pt}$')
    plt.legend(['T=21','T=22','T=23','T=24','T=25'])
    plt.grid(True)
    plt.title('$<\eta_c|\overline{c}c|\eta_c>$, 996 configurations')
    plt.show()
    
    return()
    
#plot_data_19()
plot_data_21()
plot_data_22()
plot_data_24()

plot_data_23()
plot_data_25()
