# -*- coding: utf-8 -*-
"""
Spyder Editor
@author: mfinkel
"""

import numpy as np
from itsdangerous import NoneAlgorithm

import Modells
import matplotlib.pyplot as plt
# import mem_profile
#plt.rc('xtick', labelsize=20)
#plt.rc('ytick', labelsize=20)
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 11}
plt.rc('font', **font)

#plt.rc('ylabel', labelsize=20)
#from scipy.optimize import curve_fit


import math
#import scipy
#import os
import methods
#import DataCursor
#import lgs
import Fittingtools


def linfunk(x, m, b):
    y = m*x+b
    return y

path_to_data  =  "..\\Daten-bearbeitet\\Stahl ST37\\"
odf_name = "ST37_MTODF.txt"

Data_Iron = methods.Data(path_to_data, odf_name, 6)
Data_Iron.read_scattering_data(path_of_unstraind_data="Euler-Scans ohne Last\\",
                               path_of_straind_data="Euler-Scans unter 5kN\\")
modi = ["reus", "voigt", "hill", "eshelby"]
# for i in modi:
Data_Iron.Fit_the_data_with_texture(filename="Result_iron_", method=modi[0], number_of_datapoints=None)
# Data_Iron.plot_odf()
# Data_Iron.integral_over_total_odf
# Data_Iron.calc__deltavals()

plt.show()
quit()


'''
read Data:
filelist1 contains the unstrainddata,
filelist2 the data under strain
'''

from glob import glob
filelist1 = glob('..\\Daten-bearbeitet\\Stahl ST37\\Euler-Scans ohne Last\\*.eth')
filelist2 = glob('..\\Daten-bearbeitet\\Stahl ST37\\Euler-Scans unter 5kN\\*.eth')
filelist  = filelist1 + filelist2
filelist1.sort()
filelist2.sort()
filelist.sort()





'''
-------------------------------------------------------------------------------------------
two lists of dataset objects containing the straind and the unstraind data
-------------------------------------------------------------------------------------------
'''
#list of unstraind data
unstraind=[]
for i in filelist1:
    unstraind.append(methods.dataset(i))

#list of sraind data
straind=[]
for i in filelist2:
    straind.append(methods.dataset(i))


'''
-------------------------------------------------------------------------------------------
select the peaks and calculate 2Theta
hkl_setting is a list of lists. The elements of hkl_setting contain
- the hkl (first three indiices)
- the index of the minimal and maximal value of the respectiv peak.

with hkl_setting = unstraind[0].select_hkl() it is possible to select the peaks
with some kinde of user interface
stand 08.12.2015:
the userinterface is not vary sophisticated
-------------------------------------------------------------------------------------------
'''
#'normal' Iron
#hkl setting 3 peaks [[1.0, 1.0, 0.0, 780, 980], [2.0, 0.0, 0.0, 1180, 1400], [2.0, 1.0, 1.0, 1540, 1740]]
#hkl_setting 5 peaks [[1.0, 1.0, 0.0, 852, 904], [2.0, 0.0, 0.0, 1260, 1314], [2.0, 1.0, 1.0, 1600, 1676], [2.0, 2.0, 0.0, 1926, 2020], [3.0, 1.0, 0.0, 2266, 2380]]
hkl_setting = [[1.0, 1.0, 0.0, 852, 904], [2.0, 0.0, 0.0, 1260, 1314], [2.0, 1.0, 1.0, 1600, 1676], [2.0, 2.0, 0.0, 1926, 2020], [3.0, 1.0, 0.0, 2266, 2380]]
print hkl_setting
#do it for the rest automaticaly
for i in xrange(0,len(unstraind)):
    unstraind[i].select_hkl(auto=True, rang=hkl_setting)
    straind[i].select_hkl(auto=True, rang=hkl_setting)

print unstraind[0].hkl_TTheta
print straind[0].hkl_TTheta



'''
-------------------------------------------------------------------------------------------
shift the straind data by the an offset due to missalignement of the sample calculated
out of the unstraind data
-------------------------------------------------------------------------------------------
'''

'''
latice_param = 2.8679779
wavelength = 1.5484
HKL =[[1,1,0],[2,0,0],[2,1,1]]
filename = "shiftvals.dat"
dat = open(filename,"w")
r45=[]
r60=[]
r75=[]
chi=[0,15,30,45,60,75,90]
hkl_dTheta_unstraind = []#[[label, (dTheta of: )[110,200,211]]]
hkl_dTheta_straind = []#[[label, (dTheta of: )[110,200,211]]]
for i in xrange(len(unstraind)):
    res = []
    res1=[]
    dTheta =[]
    dTheta_straind = []
    label_unstraind = str(unstraind[i].Omega) + " " + str(unstraind[i].Chi)
    label_straind = str(straind[i].Omega) + " " + str(straind[i].Chi)
    for j in HKL:
        h=j[0]
        k=j[1]
        l=j[2]
        dTheta.append(unstraind[i].calc_dTheta(latice_param, wavelength, h, k, l))
        dTheta_straind.append(straind[i].calc_dTheta(latice_param, wavelength, h, k, l))
        #dTheta=unstraind[i].calc_dTheta(latice_param, wavelength, h, k, l)
    #dTheta = np.average(dTheta)

    hkl_dTheta_unstraind.append([label_unstraind,dTheta])
    hkl_dTheta_straind.append([label_straind,dTheta_straind])


    #shift each hkl by its own dTheta-value
    straind[i].write_shifted_theta_data(dTheta[0],dTheta[1],dTheta[2])


    dTheta = np.average(dTheta)
    dTheta = 0.054912698 + dTheta
    res.append([j,dTheta])
    res1.append(dTheta)
    if unstraind[i].Omega ==45:
        r45.append(dTheta)
    if unstraind[i].Omega ==60:
        r60.append(dTheta)
    if unstraind[i].Omega ==75:
        r75.append(dTheta)
    dat.write('Omg: %i, Chi: %i; hkl:%i%i%i  shift: %f\n'%(unstraind[i].Omega,unstraind[i].Chi,h,k,l, dTheta))
    #plt.figure("hallo")
    #plt.plot(res1)
    #plt.show()
    print res
dat.close()
print r45

plt.figure()
plt.plot(r45)

plt.figure()
plt.plot(r60)

plt.figure()
plt.plot(r75)

plt.figure("dTheta, unstraind")
for i in hkl_dTheta_unstraind:
    plt.plot(i[1], label=i[0])
plt.legend(loc='lower right', numpoints = 1)

plt.figure("dTheta, straind")
for i in hkl_dTheta_straind:
    plt.plot(i[1], label=i[0])
plt.legend(loc='lower right', numpoints = 1)
plt.show()
'''




'''
-------------------------------------------------------------------------------------------
store objekt's, containing one hkl each, the angles with respect to the labframe
and the eulerangles psi and PHI with respekct to the specimen frame.
-------------------------------------------------------------------------------------------
'''
Force =5000
diameter = 6.
nu = 0#.30
HKL_object_list=[]
for i in xrange(0, len(unstraind[0].hkl_TTheta)):#loop over all measured hkl
    a =[]
    #print unstraind[0].hkl_TTheta[i]
    for j in xrange(0,len(unstraind)):#loop over all datapoints
        sigma       = methods.sigma3(Force, diameter) #stress
        chi         = unstraind[j].Chi  #chi of j'th datapoint with respect to labsys
        chi_        = 90. #chi of scattering vector
        Omega       = unstraind[j].Omega #omega of the j'th datapoint
        TTheta_0    = unstraind[j].hkl_TTheta[i][3] #theta_0 of the j'th datapoint and the i'th hkl
        TTheta_0_weight = unstraind[j].hkl_TTheta[i][4]
        TTheta      = straind[j].hkl_TTheta[i][3]   #theta of the j'th datapoint and the i'th hkl
        TTheta_weight = straind[j].hkl_TTheta[i][4]
        phi         = 0. #phi rot of the specimen

        h = unstraind[j].hkl_TTheta[i][0]
        k = unstraind[j].hkl_TTheta[i][1]
        l = unstraind[j].hkl_TTheta[i][2]

        lam = 0.154840#wavelength of SPODI


        if math.isnan(TTheta_0) or math.isnan(TTheta):
            pass
        else:
            b=methods.one_HKL(h, k, l, TTheta_0, TTheta_0_weight, TTheta, TTheta_weight,
                              chi, chi_, Omega, phi, sigma, lam)
            if (math.isnan(b.psi)!=True and math.isnan(b.PHI)!=True):
                a.append(b)
                #print TTheta_0
    HKL_object_list.append(a)



def Plot_chi_Psi(HKL_object_list, hkl, omeg, fig=True):
    '''
        Plots psi over chi to check if the calculation of psi makes sence
    '''
    hkl_ = {'110': 0,
            '200': 1,
            '211': 2,
            '220': 3,
            '310': 4}

    HKL = hkl_[hkl]
    chi = []
    psi = []
    for i in HKL_object_list[HKL]:
        if i.Omega == omeg:
            chi.append(i.chi)
            psi.append(methods.r_t_d(i.psi))

    if fig == True:
        plt.figure('%s chi - psi plot, omeg=%d' % (hkl, omeg))
        plt.plot(chi, psi, label=str(omeg))


    else:
        plt.plot(chi, psi)
    plt.xlabel('chi')
    plt.ylabel('psi')
    plt.legend()





def plot_chi2psi(HKL_object_list, hkl):
    '''
    make an sin^2 Psi Plot for the hkl (string for example '110')
    and calculate the REK's
    '''
    hkl_ = {'110': 0,
            '200': 1,
            '211': 2,
            '220': 3,
            '310': 4}
    HKL = hkl_[hkl]

    hkl_200 = HKL_object_list[HKL]
    dikt = {0: 'b', 1: 'g', 2: 'y'}
    k = 0
    xxx = []
    yyy = []
    for j in [45, 60, 75]:

        sig = 0
        e_s = 0
        e_p = 0.
        e = []
        psi = []
        chi__ = []
        for i in hkl_200:
            if (i.Omega == j):
                sig = i.sig_L()
                sigma = i.sigma
                if i.chi == 0:
                    e_s = i.d_epsilon
                if i.chi == 90:
                    e_p = i.d_epsilon

                e.append(i.d_epsilon / sigma)
                psi.append(np.sin(np.pi / 2 - i.psi) ** 2)  # attention, som crazy trick
                chi__.append(i.chi)

        sigma = hkl_200[0].sigma
        print "hkl: {4}, omega: {0}, \nchi: \n{1}, \n psi:\n{3}, \ne:\n{2}\n\n".format(j, chi__, e, psi, hkl)
        e = e[0:-1]  # -1]
        psi = psi[0:-1]  # -1]
        xxx = xxx + psi
        yyy = yyy + e
        # print 'Sigma', hkl_200[0].sig_L(), sigma
        # print 'Omega: ', j, "E: ", sig/e_p, "\nnü/E: ", e_s/sig, " nü:  ", -e_s/e_p
        try:
            a, b = Fittingtools.lin_fit(psi, e)
            y = Fittingtools.lin(np.array(psi), a, b)

            a1 = sig / e_p
            b1 = 0
            y1 = Fittingtools.lin(np.array([0, e_p]), a1, b1)

            col1 = dikt[k] + 'o'
            col2 = dikt[k] + '-'
            plt.figure('sin^2psi-epsilon plot, hkl: %s' % (hkl))
            plt.plot(psi, e, col1, label='%s, omg=%d, s1= %E, s2= %E' % (hkl, j, b, a))
            plt.plot(psi, y, col2)
            plt.legend(loc='lower left', numpoints=1)
            plt.xlabel(r'$sin^2(\Psi)$')
            plt.ylabel(r'$\varepsilon/\sigma$')

            plt.figure('sigma-epsilon, hkl: %s' % (hkl))
            plt.plot([0, e_p], [0, sig], col1, label='%s, omg=%d, E= %.2f' % (hkl, j, a1 / np.power(10., 9)) + 'GPa')
            plt.plot([0, e_p], y1, col2)
            plt.xlabel(r'$\epsilon_{\| \|}$')
            plt.ylabel(r'$\sigma$')
            plt.legend(loc='lower right', numpoints=1)
            k = k + 1
        except TypeError:
            pass
    a,b =  Fittingtools.lin_fit(xxx,yyy)
    y = Fittingtools.lin(np.array(xxx),a,b)
    plt.figure('sin^2psi-epsilon plot, hkl: %s' %(hkl))
    plt.plot(xxx,y,'r-',label='%s, s1= %E, s2= %E'%(hkl, b, a))
    plt.legend(loc='lower left', numpoints = 1)
    plt.xlabel(r'$sin^2(\Psi)$')
    plt.ylabel(r'$\varepsilon/\sigma$')
    return hkl, a, b


'''
-------------------------------------------------------------------------------
calculate s1 and 1/2s2 with the function plot_chi2psi and do
afterwords the leastsquare fit afrer Gnäupel-Herold
-------------------------------------------------------------------------------
'''

HKLs =['110', '200', '211', '220', '310']
omegs = [45, 60, 75]
for i in HKLs:
    plt.figure('%s chi - psi plot'%(i))
    for j in omegs:
        Plot_chi_Psi(HKL_object_list, hkl=i, omeg=j, fig=False)
    plt.legend(loc='lower right', numpoints = 1)
plt.show()


print "\nhkl = 110"
hkl = []
s1 = []
s2 = []
for i in HKLs:
    temp = plot_chi2psi(HKL_object_list, i)
    HHKl = [float(temp[0][0]), float(temp[0][1]), float(temp[0][2])]
    hkl.append(HHKl)
    s1.append(temp[2])
    s2.append(temp[1])
Data_Fe = Modells.InParams()
Data_Fe.set_data(hkl    = hkl,
                 s1     = s1,
                 s2     = s2,
                 s1_err = None,
                 s2_err = None)

print Data_Fe.data['hkl']
print Data_Fe.data['s1']
print Data_Fe.data['s2']

Modells.call_all(Data_Fe, 'Iron_chi^2-plot', filename = 'My_results', mode='w')

'''
-------------------------------------------------------------------------------
#calculate s1(hkl) and 1/2s2(hkl) Gnäupel herold eq (1)
-> system of linear equations -> s1(hkl), 1/2s2(hkl)
-> average over many diverent orientations to get good statistic

afterwards do the least square fit to optain the single crystal compliances
-------------------------------------------------------------------------------
'''
result = []


for s in HKL_object_list:
    res=[]
    a45=[]
    a60=[]
    a75=[]
    epsilon45=[]
    epsilon60=[]
    epsilon75=[]
    for i in s:
        h = i.h
        k = i.k
        l = i.l
        hkl=[h,k,l]


        #calculate a i
        omega= i.Omega
        if omega==45:

            epsilon45.append(i.d_epsilon)

            s2_2_const = i.s2_2_const(nu)
            s1_const = i.s1_const(nu)

            a45.append([s1_const, s2_2_const])


        elif omega==60:

            epsilon60.append(i.d_epsilon)

            s2_2_const = i.s2_2_const(nu)
            s1_const = i.s1_const(nu)

            a60.append([s1_const, s2_2_const])


        elif omega==75:

            epsilon75.append(i.d_epsilon)

            s2_2_const = i.s2_2_const(nu)
            s1_const = i.s1_const(nu)

            a75.append([s1_const, s2_2_const])


    epsilon = [epsilon45,epsilon60,epsilon75]
    a=[a45,a60,a75]

   #calculate the solutionfor the different omega:
    for m in xrange(0,len(a)):
        for g in xrange(0,len(a[m])):
            A=[]
            E=[]
            for b in xrange(int(g-2),g):
                A.append(a[m][b])
                E.append(epsilon[m][b])

           # print "\nh,k,l:  ",h,k,l, "  " ,"Epsilon: ",E, "\na:\n", A\
           #           ,"\nsolution: ", np.linalg.solve(A,E),"\n"
            res.append([hkl,list(np.linalg.solve(A,E))])

    result.append(res)
#print np.array(res)

#calculate s1 and 1/2s2



result = np.array(result)
print 'result: \n', result
s1_110=[]
s1_200=[]
s1_211=[]
s1_220=[]
s1_310=[]
s2_110=[]
s2_200=[]
s2_211=[]
s2_220=[]
s2_310=[]
#110:
for i in result[0]:
    if (i[1][0]<0 and i[1][1]>0.):
        s1_110.append(i[1][0])
        s2_110.append(i[1][1])
        print 's1_110: ',i[1][0], i[1][1]
#200
for i in result[1]:
    if (i[1][0]<0 and i[1][1]>0.):
        s1_200.append(i[1][0])
        s2_200.append(i[1][1])
        print 's1_200: ',i[1][0], i[1][1]
#211
for i in result[2]:
    if (i[1][0]<0 and i[1][1]>0.):
        s1_211.append(i[1][0])
        s2_211.append(i[1][1])
        print 's1_211: ',i[1][0], i[1][1]
#220
for i in result[3]:
    if (i[1][0]<0and i[1][1]>0):
        s1_220.append(i[1][0])
        s2_220.append(i[1][1])
        print 's1_220: ',i[1][0], i[1][1]
#310
for i in result[4]:
    if (i[1][0]<0 and i[1][1]>0.):
        s1_310.append(i[1][0])
        s2_310.append(i[1][1])
        print 's1_310: ',i[1][0], i[1][1]
    #print i[0][1][0], i[1][1][0], i[2][1][0], i[0][1][1], i[1][1][1], i[2][1][1]
data_Fe = Modells.InParams()
data_Fe.set_data(hkl =[[1,1,0],[2,0,0],[2,1,1],[2,2,0],[3,1,0]],\
        s1  = [np.average(s1_110),
                                np.average(s1_200),
                                np.average(s1_211),
                                np.average(s1_220),
                                np.average(s1_310)],\
        s2=   [np.average(s2_110),
                                np.average(s2_200),
                                np.average(s2_211),
                                np.average(s2_220),
                                np.average(s2_310)],\
        s1_err = [np.std(s1_110),
                                np.std(s1_200),
                                np.std(s1_211),
                                np.std(s1_220),
                                np.std(s1_310)],\
        s2_err = [np.std(s2_110),
                                np.std(s2_200),
                                np.std(s2_211),
                                np.std(s2_220),
                                np.std(s2_310)])
print "s1 110:\n  "  , np.average(s1_110), np.std(s1_110),\
      "\ns1_200:\n  ", np.average(s1_200), np.std(s1_200),\
      "\ns1_211:\n  ", np.average(s1_211), np.std(s1_211),\
      "\ns1_220:\n  ", np.average(s1_220), np.std(s1_220),\
      "\ns1_310:\n  ", np.average(s1_310), np.std(s1_310),\
      "\ns2_110:\n  ", np.average(s2_110), np.std(s2_110),\
      "\ns2_200:\n  ", np.average(s2_200), np.std(s2_200),\
      "\ns2_211:\n  ", np.average(s2_211), np.std(s2_211),\
      "\ns2_220:\n  ", np.average(s2_220), np.std(s2_220),\
      "\ns2_310:\n  ", np.average(s2_310), np.std(s2_310)

Modells.call_all(data_Fe, 'Iron', filename = 'My_results', mode='a')
#print "end"



plt.show()

