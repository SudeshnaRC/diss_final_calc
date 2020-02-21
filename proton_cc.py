# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 21:28:51 2019

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


data_etac = make_data("relabelled-etac-fullwall-testset.gpl")
data_p_fold = make_data("relabelled-proton-c-fold.gpl")
data_p_unfold = make_data("relabelled-proton-c-unfold.gpl")
#set correlator with tag 'etac_t0+2' as main eta_c correlator
    
def make_proton_fold():
    proton_fold = data_p_fold['Proton_s2_s2']
    return(proton_fold)
    
def make_proton_unfold():
    proton_unfold = data_p_unfold['Proton_s2_s2']
    return(proton_unfold)

#keep all other tags in matrix cc to make cc_loop
def make_cc():
    #cc = np.zeros((10,966)) 
    cc=np.zeros((10,820))
    for i in range(1,11):
        tag = "etac_t0+{0}".format(2*i)  
        #for j in range(966):
        for j in range(820):
            for k in range(64):
                cc[i-1][j]+=data_etac[tag][j][k]
                
    return(cc)
 

cc = make_cc()
proton_fold = make_proton_fold()
proton_unfold = make_proton_unfold()
    
print(cc)

print(np.shape(proton_fold))
print(np.shape(proton_unfold))

def correlator_fold():
    correlator = collections.OrderedDict()
    for i in range(1,11):
        correlator['proton'] = np.zeros((820,64))
        correlator['cc{0}'.format(i)] = np.zeros(820)
        correlator['proton_cc{0}'.format(i)] = np.zeros((820,64))
        #correlator['proton'] = np.zeros((433,64))
        #correlator['cc{0}'.format(i)] = np.zeros(433)
        #correlator['proton_cc{0}'.format(i)] = np.zeros((433,64))
        for t in range(64):
            for j in range(820):
            #for j in range(433):
                correlator['proton_cc{0}'.format(i)][j][t] += proton_fold[j][t]*cc[i-1][j]
        correlator['proton'] = proton_fold[0:820][:]
        #correlator['proton'] = proton_fold[0:433][:]
        correlator['cc{0}'.format(i)] = cc[i-1]
    AvCorr = gv.dataset.avg_data(correlator)
        #print(AvCorr)
    return(AvCorr)
AvCorr_Fold = correlator_fold()

"""
def correlator_unfold():
    correlator = collections.OrderedDict()
    for i in range(1,11):
        correlator['proton'] = np.zeros((966,64))
        correlator['cc{0}'.format(i)] = np.zeros(966)
        correlator['proton_cc{0}'.format(i)] = np.zeros((966,64))
        for t in range(64):
            for j in range(966):
                correlator['proton_cc{0}'.format(i)][j][t] += proton_unfold[j][t]*cc[i-1][j]
        correlator['proton'] = proton_unfold[0:966][:]
        correlator['cc{0}'.format(i)] = cc[i-1]
    AvCorr = gv.dataset.avg_data(correlator)
        #print(AvCorr)
    return(AvCorr)
    
AvCorr_Unfold = correlator_unfold()
"""

def plot_data_20():
    T_samp = 20
    T = np.zeros(10)
    mean = np.zeros(10)
    yerr = np.zeros(10)
    plt.figure(1)
    for i in range(1,11):
        T[i-1]=i
        mean[i-1] = -1*0.6382*(48**3)*((AvCorr_Fold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Fold['proton'][T_samp]*AvCorr_Fold['cc{0}'.format(i)])/(AvCorr_Fold['proton'][T_samp])).mean
        yerr[i-1] = 0.6382*(48**3)*((AvCorr_Fold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Fold['proton'][T_samp]*AvCorr_Fold['cc{0}'.format(i)])/(AvCorr_Fold['proton'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T, mean, yerr=[yerr,yerr], capsize=2, fmt = '-ko')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    plt.show()
    
    return()
    
def plot_data_22():
    T_samp = 22
    T = np.zeros(10)
    mean = np.zeros(10)
    yerr = np.zeros(10)
    plt.figure(1)
    for i in range(1,11):
        T[i-1]=i
        mean[i-1] = -1*0.6382*(48**3)*((AvCorr_Fold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Fold['proton'][T_samp]*AvCorr_Fold['cc{0}'.format(i)])/(AvCorr_Fold['proton'][T_samp])).mean
        yerr[i-1] = 0.6382*(48**3)*((AvCorr_Fold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Fold['proton'][T_samp]*AvCorr_Fold['cc{0}'.format(i)])/(AvCorr_Fold['proton'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T, mean, yerr=[yerr,yerr], capsize=2, fmt = '-go')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    plt.show()
    
    return()
    
def plot_data_23():
    T_samp = 23
    T = np.zeros(10)
    mean = np.zeros(10)
    yerr = np.zeros(10)
    plt.figure(1)
    for i in range(1,11):
        T[i-1]=i
        mean[i-1] = -1*0.6382*(48**3)*((AvCorr_Fold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Fold['proton'][T_samp]*AvCorr_Fold['cc{0}'.format(i)])/(AvCorr_Fold['proton'][T_samp])).mean
        yerr[i-1] = 0.6382*(48**3)*((AvCorr_Fold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Fold['proton'][T_samp]*AvCorr_Fold['cc{0}'.format(i)])/(AvCorr_Fold['proton'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T, mean, yerr=[yerr,yerr], capsize=2, fmt = '-mo')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    plt.show()
    
    return()
    
def plot_data_24():
    T_samp = 24
    T = np.zeros(10)
    mean = np.zeros(10)
    yerr = np.zeros(10)
    plt.figure(1)
    for i in range(1,11):
        T[i-1]=i
        mean[i-1] = -1*0.6382*(48**3)*((AvCorr_Fold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Fold['proton'][T_samp]*AvCorr_Fold['cc{0}'.format(i)])/(AvCorr_Fold['proton'][T_samp])).mean
        yerr[i-1] = 0.6382*(48**3)*((AvCorr_Fold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Fold['proton'][T_samp]*AvCorr_Fold['cc{0}'.format(i)])/(AvCorr_Fold['proton'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
        print(yerr)
    plt.errorbar(2*T, mean, yerr=[yerr,yerr], capsize=2, fmt = '-ro')
    plt.xticks(np.arange(0,22,step=2))
    plt.grid(True)
    plt.legend(['T=21','T=22','T=24'])
    plt.title('$<p|\overline{c}c|p>,\ C_{22}$, 966 configurations.')
    plt.xlim(1,21)
    plt.xlabel('t\'(a)')
    plt.ylabel('$C_{3pt}/C_{2pt}$')
    plt.show()
    
    return()
    
def plot_data_25():
    T_samp = 25
    T = np.zeros(10)
    mean = np.zeros(10)
    yerr = np.zeros(10)
    plt.figure(1)
    for i in range(1,11):
        T[i-1]=i
        mean[i-1] = -1*0.6382*(48**3)*((AvCorr_Fold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Fold['proton'][T_samp]*AvCorr_Fold['cc{0}'.format(i)])/(AvCorr_Fold['proton'][T_samp])).mean
        yerr[i-1] = 0.6382*(48**3)*((AvCorr_Fold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Fold['proton'][T_samp]*AvCorr_Fold['cc{0}'.format(i)])/(AvCorr_Fold['proton'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')
    print(mean)
    print(yerr)
    plt.errorbar(2*T, mean, yerr=[yerr,yerr], capsize=2, fmt = '-bo')
    plt.xlabel('$t-t_0$')
    plt.xticks(2,20,2)
    plt.xlim(1,21)
    plt.ylabel('$C_3/C_2$')
    plt.grid(True)
    plt.legend(['T=21','T=22','T=24'])
    plt.title('$<p|\overline{c}c|p>$, s2.s2, folded')
    plt.show()
    
    return()

"""
def plot_data_18_u():
    T_samp = 18
    T = np.zeros(10)
    mean = np.zeros(10)
    yerr = np.zeros(10)
    plt.figure(2)
    for i in range(1,11):
        T[i-1]=i
        mean[i-1] = -1*0.6382*(48**3)*((AvCorr_Unfold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Unfold['proton'][T_samp]*AvCorr_Unfold['cc{0}'.format(i)])/(AvCorr_Unfold['proton'][T_samp])).mean
        yerr[i-1] = 0.6382*(48**3)*((AvCorr_Unfold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Unfold['proton'][T_samp]*AvCorr_Unfold['cc{0}'.format(i)])/(AvCorr_Unfold['proton'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T, mean, yerr=[yerr,yerr], capsize=2, fmt = '-bo')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    plt.show()
    
def plot_data_20_u():
    T_samp = 20
    T = np.zeros(10)
    mean = np.zeros(10)
    yerr = np.zeros(10)
    plt.figure(2)
    for i in range(1,11):
        T[i-1]=i
        mean[i-1] = -1*0.6382*(48**3)*((AvCorr_Unfold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Unfold['proton'][T_samp]*AvCorr_Unfold['cc{0}'.format(i)])/(AvCorr_Unfold['proton'][T_samp])).mean
        yerr[i-1] = 0.6382*(48**3)*((AvCorr_Unfold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Unfold['proton'][T_samp]*AvCorr_Unfold['cc{0}'.format(i)])/(AvCorr_Unfold['proton'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T, mean, yerr=[yerr,yerr], capsize=2, fmt = '-ro')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    plt.show()
    
    
def plot_data_22_u():
    T_samp = 22
    T = np.zeros(10)
    mean = np.zeros(10)
    yerr = np.zeros(10)
    plt.figure(2)
    for i in range(1,11):
        T[i-1]=i
        mean[i-1] = -1*0.6382*(48**3)*((AvCorr_Unfold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Unfold['proton'][T_samp]*AvCorr_Unfold['cc{0}'.format(i)])/(AvCorr_Unfold['proton'][T_samp])).mean
        yerr[i-1] = 0.6382*(48**3)*((AvCorr_Unfold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Unfold['proton'][T_samp]*AvCorr_Unfold['cc{0}'.format(i)])/(AvCorr_Unfold['proton'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T, mean, yerr=[yerr,yerr], capsize=2, fmt = '-go')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    plt.show()
    
def plot_data_24_u():
    T_samp = 24
    T = np.zeros(10)
    mean = np.zeros(10)
    yerr = np.zeros(10)
    plt.figure(2)
    for i in range(1,11):
        T[i-1]=i
        mean[i-1] = -1*0.6382*(48**3)*((AvCorr_Unfold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Unfold['proton'][T_samp]*AvCorr_Unfold['cc{0}'.format(i)])/(AvCorr_Unfold['proton'][T_samp])).mean
        yerr[i-1] = 0.6382*(48**3)*((AvCorr_Unfold['proton_cc{0}'.format(i)][T_samp]-AvCorr_Unfold['proton'][T_samp]*AvCorr_Unfold['cc{0}'.format(i)])/(AvCorr_Unfold['proton'][T_samp])).sdev
    #plt.plot(2*T, mean, 'ro')  
    plt.errorbar(2*T, mean, yerr=[yerr,yerr], capsize=2, fmt = '-mo')
    plt.xlabel('$t-t_0$')
    plt.ylabel('$c_3(t-t_0)$')
    plt.grid(True)
    plt.legend(['T=20','T=22','T=24'])
    plt.title('$<p|c\overline{c}|p>$, s2.s2, unfolded')
    plt.show()
 """   
#plot_data_18()
plot_data_20()
plot_data_22()
#plot_data_23()
plot_data_24()
#plot_data_25()

#plot_data_18_u()
#plot_data_20_u()
#plot_data_22_u()
#plot_data_24_u()