# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:53:08 2015

@author: mfinkel
"""
import lmfit as lm
import numpy as np
# from numpy.linalg import inv
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.figure import Figure
from scipy.optimize import fsolve
import scipy
# import scipy.optimize
import time as tm
# import os.path
import tempfile
import itertools as IT
import os
import sys
import gui
import PyQt4.QtGui


# import multiprocessing as multi
# from multiprocessing.pool import ThreadPool
# from multiprocessing.pool import Pool
# import Plot_coordsys as cplot


class CalcError(Exception):
    def __init__(self):
        self.value = 1
        # def __str__(self):
        #    return repr(self.value)


class InParams(object):
    """
    A dictionary of all the Parameters required to calculate teh compliance's
    """

    def __init__(self):

        # if err = None, there are now uncertaties
        self.__err = False
        self.__dict = {'hkl': None,
                       's1': None,
                       's2': None,
                       's1_err': None,
                       's2_err': None}

    def set_data(self, hkl=None, s1=None, s2=None,
                 s1_err=None, s2_err=None,
                 exp=False, EXP=1e-12
                 # if exp = False, nead not to multiplye datavalues by EXP
                 # if exp = True, multiplye data by EXP
                 ):
        """
        this memberfunction of InParams defines the data inserted by the user
        per default everythif is set to None
        if exp is False, it is not nessessary to multiply by EXP.
        else the data is multiplyed by EXP, so you don't need to type it in like
        s1 = [1e-12, 2e-12, ...] but you can do now s1 = [1,2,...] and set
        exp to True and EXP = 1e-12

        you don't need to specifie the uncertencies (s1_err and s2_err)
        """
        s1 = np.array(s1)
        s2 = np.array(s2)
        if (not (s1_err is None) and not (s2_err is None)):
            print 'Hallo'
            s1_err = np.array(s1_err)
            s2_err = np.array(s2_err)
        if exp != False:
            s1 = s1 * EXP
            s2 = s2 * EXP
            if (not (s1_err is None) and not (s2_err is None)):
                s1_err = s1_err * EXP
                s2_err = s2_err * EXP
        print 'HHKKLL: ', self.__dict['hkl']
        if not (self.__dict['hkl'] is None):
            self.__dict['hkl'] = np.append(self.__dict['hkl'], hkl)
            self.__dict['s1'] = np.append(self.__dict['s1'], s1)
            self.__dict['s2'] = np.append(self.__dict['s2'], s2)
            self.__dict['s1_err'] = np.append(self.__dict['s1_err'], s1_err)
            self.__dict['s2_err'] = np.append(self.__dict['s2_err'], s2_err)
        else:
            self.__dict['hkl'] = hkl
            self.__dict['s1'] = s1
            self.__dict['s2'] = s2
            self.__dict['s1_err'] = s1_err
            self.__dict['s2_err'] = s2_err

    def __getdata(self):
        return self.__dict

    data = property(__getdata)


'''
-----------------------------------------------------------------------------
First define the functions to compute the tehoretical values

Voigt's Model
-----------------------------------------------------------------------------
'''


def Voigt(Gamma, c_11, c_12, c_44):
    '''
        Implementation of the Voigt Model
        returns s1, s2, ...
    '''
    x = 3 * c_11
    y = 3 * c_12
    z = 3 * c_44
    s_1_v = -(3 / 2.) * ((x + 4 * y - 2 * z) / ((x - y + 3 * z) * (x + 2 * y)))
    s_2_v = (15. / (2 * x - 2 * y + 6 * z))
    res = []
    try:
        for i in range(len(Gamma)):
            res.append([s_1_v, s_2_v])
    except TypeError:
        res = [s_1_v, s_2_v]
    res = np.array(res)
    return res.flatten()


'''
-----------------------------------------------------------------------------
Reus's Model
-----------------------------------------------------------------------------
'''


def Reus(Gamma, c_11, c_12, c_44=126.4 * np.power(10., 9)):
    '''
        Implimentation of the Reus Model
    '''
    a_11 = (c_11 ** 2 - c_12 ** 2) / (c_11 ** 3 + 2 * c_12 ** 3 - 3 * c_11 * (c_12 ** 2))
    a_12 = (c_12 ** 2 - c_12 * c_11) / (c_11 ** 3 + 2 * c_12 ** 3 - 3 * c_11 * (c_12 ** 2))
    a_44 = 1. / c_44

    s_1 = a_12 + (a_11 - a_12 - 0.5 * a_44) * Gamma

    s_2 = a_11 - a_12 - 3 * (a_11 - a_12 - 0.5 * a_44) * Gamma

    res = np.array([s_1, s_2])
    return res.flatten()


'''
-----------------------------------------------------------------------------
Hill's Model
-----------------------------------------------------------------------------
'''


def RV(Gamma, c_11, c_12, c_44=126.4 * np.power(10., 9)):
    '''
        Implementation of Hill's Model
    '''

    a_11 = (c_11 ** 2 - c_12 ** 2) / (c_11 ** 3 + 2 * c_12 ** 3 - 3 * c_11 * c_12 ** 2)

    a_12 = (c_12 ** 2 - c_12 * c_11) / (c_11 ** 3 + 2 * c_12 ** 3 - 3 * c_11 * c_12 ** 2)

    a_44 = 1. / c_44

    s_1_r = a_12 + (a_11 - a_12 - 0.5 * a_44) * Gamma

    s_2_r = a_11 - a_12 - 3 * (a_11 - a_12 - 0.5 * a_44) * Gamma

    x = 3 * c_11
    y = 3 * c_12
    z = 3 * c_44
    s_1_v = -(3 / 2.) * ((x + 4 * y - 2 * z) / ((x - y + 3 * z) * (x + 2 * y)))
    s_2_v = (15. / (2 * x - 2 * y + 6 * z))

    s_1 = (s_1_r + s_1_v) / 2
    s_2 = (s_2_r + s_2_v) / 2
    return np.array([s_1, s_2]).flatten()


'''
-----------------------------------------------------------------------------
calculate shearmodulus
-----------------------------------------------------------------------------
'''


def GG(G, a, b, c):
    '''
    calculate the shearmodulus
    '''
    return G ** 3 + a * G ** 2 + b * G + c


def Gama(h, k, l):
    '''
    Function to calc Gamma
    '''
    g = ((h ** 2) * (k ** 2) + (h ** 2) * (l ** 2) + (k ** 2) * (l ** 2)) / (((h ** 2) + (k ** 2) + (l ** 2)) ** 2)
    return g


'''
-----------------------------------------------------------------------------
de Wit's Model
-----------------------------------------------------------------------------
'''


def BHM_dW(Gamma, c_11=222.6 * np.power(10., 9), c_12=123.2 * np.power(10., 9), c_44=121.7 * np.power(10., 9)):
    '''
        an Implementation of the Theory of Bollenrath, Hauk and Müller
        with the parrameters corrected by de Wit
    '''
    K = (c_11 + 2 * c_12) / 3.

    mu = c_44
    etha = (c_11 - c_12) / 2.

    a = 3. / 8. * (3. * K + 4. * (mu + 3. * (etha - mu) * Gamma)) - 1. / 5. * (2. * etha + 3. * mu)
    b = 3. / 4. * K * (mu + 3. * (etha - mu) * Gamma) - 3. / 40. * (6. * K * etha + 9. * K * mu + 20. * etha * mu)
    c = (-(3 * K * mu * etha) / 4.)
    G = fsolve(GG, 79.3 * np.power(10., 9.), args=(a, b, c), maxfev=1000000)
    s_1_b = ((1 / (9. * K)) - (1 / (6. * G)))
    s_2_b = (1 / (2. * G))
    res = np.array([s_1_b, s_2_b])
    return res.flatten()


'''
-----------------------------------------------------------------------------
Model nach Bollenrath, Hauk und Müller
-----------------------------------------------------------------------------
'''


def BHM(Gamma, c_11=217.8 * np.power(10., 9), c_12=125.8 * np.power(10., 9), c_44=126.4 * np.power(10., 9)):
    '''
        an Implementation of the Theory of Bollenrath, Hauk and Müller
    '''
    K_K = (c_11 + 2 * c_12) / 3.  # compresionmodulus of the Krystal (Korn)
    K_M = (c_11 + 2 * c_12) / 3.  # Kompresionsmudul of the Matrix
    mu = c_44
    etha = (c_11 - c_12) / 2.
    '''
    Calculation of the shearmodulus
    using th formula of Bollerntath, Hauk and Müller 1966:
    '''
    a = ((9 * K_K + 4 * etha) / 8.)
    b = -(((3 * K_K + 12 * etha) * mu) / 8.)
    c = -((3 * K_K * mu * etha) / 4.)
    G = float(fsolve(GG, 79.3 * np.power(10., 9.), args=(a, b, c), maxfev=1000000))

    # calculation of the REK's:
    A = (K_M - K_K) / (K_M * (3 * K_K + 4 * G))
    B = ((G - etha) * (3 * K_M + 6 * G)) / (G * (8 * G ** 2 + G * (9 * K_M + 12 * etha) + 6 * etha * K_M))
    C = ((G - mu) * (3 * K_M + 6 * G)) / (G * (8 * G ** 2 + G * (9 * K_M + 12 * mu) + 6 * mu * K_M))
    D = 1 / (3. * K_M)
    E = 1 / (2 * G)

    '''
    calculate the REK's, using the formula of
    Bollerntath, Hauk and Müller 1966:
    '''
    s_1_b = (D - E) / 3. + (A - B) / 3. + Gamma * (B - C)
    s_2_b = E + B - 3 * Gamma * (B - C)

    res = np.array([s_1_b, s_2_b])
    return res.flatten()


def Voigt__(Gamma, c_11, c_12, c_44=126.4 * np.power(10., 9)):
    '''
        returns s1, s2, s1, s2, ...
    '''
    x = 3 * c_11
    y = 3 * c_12
    z = 3 * c_44
    s_1_v = -(3 / 2.) * ((x + 4 * y - 2 * z) / ((x - y + 3 * z) * (x + 2 * y)))
    s_2_v = (15. / (2 * x - 2 * y + 6 * z))
    res = []
    try:
        for i in range(len(Gamma)):
            res.append([s_1_v, s_2_v])
    except TypeError:
        res = [s_1_v, s_2_v]
    res = np.array(res)
    return res


'''
def fitting_Voigt(Gamma,y, w):
    Vmod = lm.Model(RV, independent_vars=['Gamma'])
    Vmod.set_param_hint('c_11', value=217.8*np.power(10.,9))
    Vmod.set_param_hint('c_12', value=125.8*np.power(10.,9))
    Vmod.set_param_hint('c_44', value=126.4*np.power(10.,9))
    pars = Vmod.make_params()
    print Vmod.param_names
    y1 = Vmod.eval(Gamma=Gamma, params=pars)
    result = Vmod.fit(y,Gamma=Gamma,params=pars, weights=w)
    print result.fit_report()
    return result.fit_report()
'''
'''
-----------------------------------------------------------------------------
fitting the Hill Moddel:
-----------------------------------------------------------------------------
'''


# Function to be minimized
def residual_Hill(params, gam, data1=None, data2=None, weight1=None, weight2=None):
    '''
        Hill Model to calculate the compliance's
    '''

    c_11 = params['c_11'].value
    c_12 = params['c_12'].value
    c_44 = params['c_44'].value
    model1 = []
    model2 = []
    try:
        for i in gam:
            model1.append(RV(i, c_11, c_12, c_44)[0])  # s1
            model2.append(RV(i, c_11, c_12, c_44)[1])  # s2
    except all:
        pass
    mde = []
    if (data1 is None and data2 is None):
        return model1, model2
    if (weight1 is None and weight2 is None):
        for i in xrange(len(data1)):
            mde.append((model1[i] - data1[i]))
            mde.append((model2[i] - data2[i]))
        mde = np.array(mde)
        return mde
    mde.append((model1 - data1) / weight1)
    mde.append((model2 - data2) / weight2)
    mde = np.array(mde)
    return mde


def fitting_Hill(gam, data1, data2, weight1, weight2):
    '''
    -----------------------------------------------------------------------------
    fitting the Hill Moddel:
    -----------------------------------------------------------------------------
    '''
    # create a set of parameters
    params = lm.Parameters()
    #          name,   value                           Min,                      Max
    params.add('c_11', value=217.8 * np.power(10., 9))  # ,  min=100.*np.power(10.,9), max=300.*np.power(10.,9))
    params.add('c_12', value=125.8 * np.power(10., 9))  # ,  min=100.*np.power(10.,9), max=300.*np.power(10.,9))
    params.add('c_44', value=126.4 * np.power(10., 9))  # ,  min=100.*np.power(10.,9), max=300.*np.power(10.,9))

    result = lm.minimize(residual_Hill, params, method='leastsq', args=(gam,), \
                         kws={'data1': data1, 'data2': data2, 'weight1': weight1, 'weight2': weight2})
    pars = result
    print '#datapoints: ', pars.ndata
    print 'redchi: ', pars.redchi
    print 'chi: ', pars.chisqr
    print 'bic: ', pars.bic
    return result


'''
-----------------------------------------------------------------------------
fitting the BHM Moddel:
-----------------------------------------------------------------------------
'''


# Function to be minimized
def residual_BHM(params, gam, data1=None, data2=None, weight1=None, weight2=None):
    '''
        BHM Model to calculate the compliances
    '''

    c_11 = params['c_11'].value
    c_12 = params['c_12'].value
    c_44 = params['c_44'].value
    model1 = []
    model2 = []
    try:
        for i in gam:
            model1.append(BHM(i, c_11, c_12, c_44)[0])  # s1
            model2.append(BHM(i, c_11, c_12, c_44)[1])  # s2
    except all:
        pass
    mde = []
    if (data1 is None and data2 is None):
        return model1, model2
    if (weight1 is None and weight2 is None):
        for i in xrange(len(data1)):
            mde.append((model1[i] - data1[i]))
            mde.append((model2[i] - data2[i]))
        mde = np.array(mde)
        return mde
    mde.append((model1 - data1) / weight1)
    mde.append((model2 - data2) / weight2)
    mde = np.array(mde)
    return mde


def fitting_BHM(gam, data1, data2, weight1, weight2):
    # create a set of parameters
    params = lm.Parameters()
    #          name,   value                           Min,                      Max
    params.add('c_11', value=217.8 * np.power(10., 9))  # ,  min=100.*np.power(10.,9), max=300.*np.power(10.,9))
    params.add('c_12', value=125.8 * np.power(10., 9))  # ,  min=100.*np.power(10.,9), max=300.*np.power(10.,9))
    params.add('c_44', value=126.4 * np.power(10., 9))  # ,  min=100.*np.power(10.,9), max=300.*np.power(10.,9))

    result = lm.minimize(residual_BHM, params, method='leastsq', args=(gam,), \
                         kws={'data1': data1, 'data2': data2, 'weight1': weight1, 'weight2': weight2})

    return result


'''
-----------------------------------------------------------------------------
fitting the deWit Moddel:
-----------------------------------------------------------------------------
'''


# Function to be minimized
def residual_deWit(params, gam, data1=None, data2=None, weight1=None, weight2=None):
    '''
        deWit Model to calculate the compliances
    '''

    c_11 = params['c_11'].value
    c_12 = params['c_12'].value
    c_44 = params['c_44'].value
    model1 = []
    model2 = []
    try:
        for i in gam:
            model1.append(BHM_dW(i, c_11, c_12, c_44)[0])  # s1
            model2.append(BHM_dW(i, c_11, c_12, c_44)[1])  # s2
    except all:
        pass
    mde = []
    if (data1 is None and data2 is None):
        return model1, model2
    if (weight1 is None and weight2 is None):
        for i in xrange(len(data1)):
            mde.append((model1[i] - data1[i]))
            mde.append((model2[i] - data2[i]))
        mde = np.array(mde)
        return mde
    mde.append((model1 - data1) / weight1)
    mde.append((model2 - data2) / weight2)
    mde = np.array(mde)
    return mde


def fitting_deWit(gam, data1, data2, weight1, weight2):
    # create a set of parameters
    params = lm.Parameters()
    #          name,   value                           Min,                      Max
    params.add('c_11', value=217.8 * np.power(10., 9))  # ,  min=100.*np.power(10.,9), max=300.*np.power(10.,9))
    params.add('c_12', value=125.8 * np.power(10., 9))  # ,  min=100.*np.power(10.,9), max=300.*np.power(10.,9))
    params.add('c_44', value=126.4 * np.power(10., 9))  # ,  min=100.*np.power(10.,9), max=300.*np.power(10.,9))

    result = lm.minimize(residual_deWit, params, method='leastsq', args=(gam,), \
                         kws={'data1': data1, 'data2': data2, 'weight1': weight1, 'weight2': weight2})
    return result


'''
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
In this section I try to calculatet the strain in the measurmentdirection
(the scateringvector q) according to the book of Herfried Behnken (his habil.).
this is not very simpel and my be improfed by some clever algorithms for
selfconsistent calculation of things like microstress, phasestress
and thow on.
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
'''

'''
-----------------------------------------------------------------------------
First define a class for the ODF (the Orientation Distribution Funktion).
This class should hold the Intensity of the ODF as data.
Methods:
- Integrate_over_al_orientations
-----------------------------------------------------------------------------
'''


class FitStrainWithTexture(object):
    def __init__(self, data_object):
        self.data_object = data_object
        # self.xvals =
        # self.force = force
        # self.diameter = diameter
        # self.__strains_data = strains_data
        # self.__weights = weights
        self.__odf_phase_1 = self.data_object.odf_phase_1
        self.__odf_phase_2 = self.data_object.odf_phase_2
        self.counter2 = 0
        self.symmetry_phase_1 = self.__odf_phase_1.crystal_symmetry
        try:
            self.symmetry_phase_2 = self.__odf_phase_2.crystal_symmetry
            self.phase_flag = True  # multi phase material
        except:
            self.phase_flag = False  # single phase material

        self.odf_integral_over_all_angles_phase_1 = self.__odf_phase_1.integral_over_all_orientations_g
        if self.phase_flag:
            self.odf_integral_over_all_angles_phase_2 = self.__odf_phase_2.integral_over_all_orientations_g

        self.__params = lm.Parameters()  # Parameters of the of both phases
        self.add_params(params=self.__params, sym=self.symmetry_phase_1)  # add params of phase 1
        if self.phase_flag:
            self.add_params_phase_2(params=self.__params, sym=self.symmetry_phase_2)  # add params of phase 2

        self.__params.add("sigma_11", value=self.stress_sigma(0, 0), vary=False)
        self.__params.add("sigma_22", value=self.stress_sigma(1, 1), vary=False)
        self.__params.add("sigma_33", value=self.stress_sigma(2, 2), vary=False)
        self.__params.add("sigma_12", value=self.stress_sigma(0, 1), vary=False)
        self.__params.add("sigma_13", value=self.stress_sigma(0, 2), vary=False)
        self.__params.add("sigma_23", value=self.stress_sigma(1, 2), vary=False)

        # self.__crystal_sym = {
        #     "triklinic"       : np.zeros(21),  # 21 independent components
        #     "monoklinic"      : np.zeros(13),  # 13 independent components
        #     "rhombic"         : np.zeros( 9),  # 9  independent components
        #     "trigonal_3"      : np.zeros( 7),  # 7  independent components
        #     "trigonal_32"     : np.zeros( 6),  # 6  independent components
        #     "tetragonal_4"    : np.zeros( 7),  # 7  independent components
        #     "tetragonal_4mm"  : np.zeros( 6),  # 6  independent components
        #     "hexagonal"       : np.zeros( 5),  # 5  independent components
        #     "cubic"           : np.zeros( 3),  # 3  independent components
        #     "isotope"         : np.zeros( 2)   # 2  independent components
        #     }

        self.__constant_c_tensor_voigt_notation_fitted_phase = np.zeros((6, 6))
        self.__constant_c_tensor_extended_notation_fitted_phase = np.zeros((3, 3, 3, 3))
        self.__complience_s_tensor_voigt_notation_fitted_phase = np.zeros((6, 6))
        self.__complience_s_tensor_extended_notation_fitted_phase = np.zeros((3, 3, 3, 3))

        self.__constant_c_tensor_voigt_notation_fixed_phase = np.zeros(
            (6, 6))  # in this case, this are the tensors fore the inclousion of the same type as the Matrix
        self.__constant_c_tensor_extended_notation_fixed_phase = np.zeros((3, 3, 3, 3))
        self.__complience_s_tensor_voigt_notation_fixed_phase = np.zeros((6, 6))
        self.__complience_s_tensor_extended_notation_fixed_phase = np.zeros((3, 3, 3, 3))

    def set_odf_matrix(self, odf_path):
        new_odf = ODF()
        new_odf.read_data(odf_path)
        self.symmetry_phase_1 = new_odf.crystal_symmetry
        self.__odf_phase_1 = new_odf

    @staticmethod
    def add_params(params, sym):
        if sym == "isotope":
            params.add('c_11_p1', value=220 * np.power(10., 9), min=0 * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_12_p1', value=126 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
        elif sym == "m-3m":
            params.add('c_11_p1', value=240 * np.power(10., 9), min=10. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_12_p1', value=120 * np.power(10., 9), min=10. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_44_p1', value=120 * np.power(10., 9), min=10. * np.power(10., 9), max=600. * np.power(10., 9))
        elif sym == "hexagonal":
            params.add('c_11_p1', value=217 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_12_p1', value=120 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_13_p1', value=120 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_33_p1', value=126 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_44_p1', value=126 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))

    @staticmethod
    def add_params_phase_2(params, sym):
        if sym == "isotope":
            params.add('c_11_p2', value=220 * np.power(10., 9), min=0 * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_12_p2', value=126 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
        elif sym == "m-3m":
            params.add('c_11_p2', value=230 * np.power(10., 9), min=10. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_12_p2', value=120 * np.power(10., 9), min=10. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_44_p2', value=115 * np.power(10., 9), min=10. * np.power(10., 9), max=600. * np.power(10., 9))
        elif sym == "hexagonal":
            params.add('c_11_p2', value=217 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_12_p2', value=120 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_13_p2', value=120 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_33_p2', value=126 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_44_p2', value=126 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))

    def set_params_phase_1(self, params):
        keys = params.keys()
        for key in keys:
            for k in self.__params.keys():
                if ("p1" in k) and (key in k):
                    self.__params[k].value = params[key].value

    def set_params_phase_2(self, params):
        keys = params.keys()
        for key in keys:
            for k in self.__params.keys():
                if ("p2" in k) and (key in k):
                    self.__params[k].value = params[key].value

    def print_params(self):
        self.__params.pretty_print()

    @staticmethod
    def __voigt_notation(i, j):
        res = int
        if i == j == 1:
            res = 1
        elif i == j == 2:
            res = 2
        elif i == j == 3:
            res = 3
        elif (i == 2 and j == 3) or (i == 3 and j == 2):
            res = 4
        elif (i == 1 and j == 3) or (i == 3 and j == 1):
            res = 5
        elif (i == 2 and j == 1) or (i == 1 and j == 2):
            res = 6
        return res

    def __conv_all_voigtnot_to_extended_not(self):
        """
        not longer used
        :return:
        """
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        m = self.__voigt_notation(i + 1, j + 1)
                        n = self.__voigt_notation(k + 1, l + 1)

                        if m <= 3 and n <= 3:
                            self.__constant_c_tensor_extended_notation_fitted_phase[i, j, k, l] = \
                                self.__constant_c_tensor_voigt_notation_fitted_phase[
                                    m - 1, n - 1]
                            self.__complience_s_tensor_extended_notation_fitted_phase[i, j, k, l] = \
                                self.__complience_s_tensor_voigt_notation_fitted_phase[
                                    m - 1, n - 1]
                        elif m <= 3 < n:
                            self.__constant_c_tensor_extended_notation_fitted_phase[i, j, k, l] = \
                                self.__constant_c_tensor_voigt_notation_fitted_phase[
                                    m - 1, n - 1]
                            self.__complience_s_tensor_extended_notation_fitted_phase[i, j, k, l] = \
                                self.__complience_s_tensor_voigt_notation_fitted_phase[
                                    m - 1, n - 1] / 2.
                        elif m > 3 and n > 3:
                            self.__constant_c_tensor_extended_notation_fitted_phase[i, j, k, l] = \
                                self.__constant_c_tensor_voigt_notation_fitted_phase[
                                    m - 1, n - 1]
                            self.__complience_s_tensor_extended_notation_fitted_phase[i, j, k, l] = \
                                self.__complience_s_tensor_voigt_notation_fitted_phase[
                                    m - 1, n - 1] / 4.

    def __conv_all_extended_not_to_voigt_not(self):
        """
        not longer used
        :return:
        """
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        m = self.__voigt_notation(i + 1, j + 1)
                        n = self.__voigt_notation(k + 1, l + 1)

                        if m <= 3 and n <= 3:
                            self.__constant_c_tensor_voigt_notation_fitted_phase[m - 1, n - 1] = \
                                self.__constant_c_tensor_extended_notation_fitted_phase[
                                    i, j, k, l]
                            self.__complience_s_tensor_voigt_notation_fitted_phase[m - 1, n - 1] = \
                                self.__complience_s_tensor_extended_notation_fitted_phase[
                                    i, j, k, l]
                        elif m <= 3 < n:
                            self.__constant_c_tensor_voigt_notation_fitted_phase[m - 1, n - 1] = \
                                self.__constant_c_tensor_extended_notation_fitted_phase[
                                    i, j, k, l]
                            self.__complience_s_tensor_voigt_notation_fitted_phase[m - 1, n - 1] = 2 * \
                                                                                                   self.__complience_s_tensor_extended_notation_fitted_phase[
                                                                                                       i, j, k, l]
                        elif m > 3 and n > 3:
                            self.__constant_c_tensor_voigt_notation_fitted_phase[m - 1, n - 1] = \
                                self.__constant_c_tensor_extended_notation_fitted_phase[
                                    i, j, k, l]
                            self.__complience_s_tensor_voigt_notation_fitted_phase[m - 1, n - 1] = 4 * \
                                                                                                   self.__complience_s_tensor_extended_notation_fitted_phase[
                                                                                                       i, j, k, l]

    def __conv_voigtnot_to_extended_not_compliences_s(self, tensor_in_voigt_not):
        """
        convert the tensor of the elasic compliences from the woigt notation (m,n) to the extended notation (i,j,k,l)
        :rtype: extendet_tensor
        """
        extendet_tensor = np.zeros((3, 3, 3, 3))
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        m = self.__voigt_notation(i + 1, j + 1)
                        n = self.__voigt_notation(k + 1, l + 1)

                        if m <= 3 and n <= 3:
                            extendet_tensor[i, j, k, l] = tensor_in_voigt_not[m - 1, n - 1]
                        elif m <= 3 < n:
                            extendet_tensor[i, j, k, l] = tensor_in_voigt_not[m - 1, n - 1] / 2.
                        elif m > 3 and n > 3:
                            extendet_tensor[i, j, k, l] = tensor_in_voigt_not[m - 1, n - 1] / 4.

        return extendet_tensor

    def __conv_extended_not_to_voigt_not_compliences_s(self, tensor_in_extendet_notation):
        """
        convert the tensor of the elasic compliences from the extended notation (i,j,k,l) to the woigt notation (m,n)
        :rtype: extendet_tensor
        """
        tennsor_in_voigt_notation = np.zeros((6, 6))
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        m = self.__voigt_notation(i + 1, j + 1)
                        n = self.__voigt_notation(k + 1, l + 1)

                        if m <= 3 and n <= 3:
                            tennsor_in_voigt_notation[m - 1, n - 1] = tensor_in_extendet_notation[i, j, k, l]
                        elif m <= 3 < n:
                            tennsor_in_voigt_notation[m - 1, n - 1] = 2 * tensor_in_extendet_notation[i, j, k, l]
                        elif m > 3 and n > 3:
                            tennsor_in_voigt_notation[m - 1, n - 1] = 4 * tensor_in_extendet_notation[i, j, k, l]

        return tennsor_in_voigt_notation

    def __conv_voigtnot_to_extended_not_constants_c(self, tensor_in_voigt_not):
        """
        convert the tensor of the elastic constants from the woigt notation (m,n) to the extended notation (i,j,k,l)
        :rtype: extended_tensor
        """
        extendet_tensor = np.zeros((3, 3, 3, 3))
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        m = self.__voigt_notation(i + 1, j + 1)
                        n = self.__voigt_notation(k + 1, l + 1)

                        if m <= 3 and n <= 3:
                            extendet_tensor[i, j, k, l] = tensor_in_voigt_not[m - 1, n - 1]
                        elif m <= 3 < n:
                            extendet_tensor[i, j, k, l] = tensor_in_voigt_not[m - 1, n - 1]
                        elif m > 3 and n > 3:
                            extendet_tensor[i, j, k, l] = tensor_in_voigt_not[m - 1, n - 1]

        return extendet_tensor

    def __conv_extended_not_to_voigt_not_constants_c(self, tensor_in_extendet_notation):
        """
        convert the tensor of the elastic constants from the extended notation (i,j,k,l) to the woigt notation (m,n)
        :rtype: extended_tensor
        """
        tennsor_in_voigt_notation = np.zeros((6, 6))
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        m = self.__voigt_notation(i + 1, j + 1)
                        n = self.__voigt_notation(k + 1, l + 1)

                        if m <= 3 and n <= 3:
                            tennsor_in_voigt_notation[m - 1, n - 1] = tensor_in_extendet_notation[i, j, k, l]
                        elif m <= 3 < n:
                            tennsor_in_voigt_notation[m - 1, n - 1] = tensor_in_extendet_notation[i, j, k, l]
                        elif m > 3 and n > 3:
                            tennsor_in_voigt_notation[m - 1, n - 1] = tensor_in_extendet_notation[i, j, k, l]
        return tennsor_in_voigt_notation

    def set_parameters_of_the_FITTED_PHASE_in_matrix_representation(self, params, symmetry):
        pars = self.__return_free_or_fixed_parameters(params=params, free=True)
        if symmetry == "isotope":
            self.__constant_c_tensor_voigt_notation_fitted_phase = \
                np.array([
                    [pars['c_11'].value, pars['c_12'].value, pars['c_11'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_11'].value, pars['c_12'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_12'].value, pars['c_11'].value, 0, 0, 0],
                    [0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value), 0, 0],
                    [0, 0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value), 0],
                    [0, 0, 0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value)]
                ])

        elif symmetry == "m-3m":  # cubic
            # print "Symmetry:", self.symmetry_phase_1
            self.__constant_c_tensor_voigt_notation_fitted_phase = \
                np.array([
                    [pars['c_11'].value, pars['c_12'].value, pars['c_12'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_11'].value, pars['c_12'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_12'].value, pars['c_11'].value, 0, 0, 0],
                    [0, 0, 0, pars['c_44'].value, 0, 0],
                    [0, 0, 0, 0, pars['c_44'].value, 0],
                    [0, 0, 0, 0, 0, pars['c_44'].value]
                ])

        elif symmetry == "hexagonal" or symmetry == "hexagonal":
            self.__constant_c_tensor_voigt_notation_fitted_phase = \
                np.array([
                    [pars['c_11'].value, pars['c_12'].value, pars['c_13'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_11'].value, pars['c_13'].value, 0, 0, 0],
                    [pars['c_13'].value, pars['c_13'].value, pars['c33'].value, 0, 0, 0],
                    [0, 0, 0, pars['c_44'].value, 0, 0],
                    [0, 0, 0, 0, pars['c_44'].value, 0],
                    [0, 0, 0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value)]
                ])

        self.__constant_c_tensor_extended_notation_fitted_phase = self.__conv_voigtnot_to_extended_not_constants_c(
            self.__constant_c_tensor_voigt_notation_fitted_phase)
        self.__complience_s_tensor_voigt_notation_fitted_phase = np.linalg.inv(
            self.__constant_c_tensor_voigt_notation_fitted_phase)
        # self.__conv_all_voigtnot_to_extended_not()
        self.__complience_s_tensor_extended_notation_fitted_phase = self.__conv_voigtnot_to_extended_not_compliences_s(
            self.__complience_s_tensor_voigt_notation_fitted_phase)

    def set_parameters_of_the_FIXED_PHASE_in_matrix_representation_multi_phase(self, params, symmetry):
        pars = self.__return_free_or_fixed_parameters(params=params, free=False)
        if symmetry == "isotope":
            self.__constant_c_tensor_voigt_notation_fixed_phase = \
                np.array([
                    [pars['c_11'].value, pars['c_12'].value, pars['c_11'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_11'].value, pars['c_12'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_12'].value, pars['c_11'].value, 0, 0, 0],
                    [0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value), 0, 0],
                    [0, 0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value), 0],
                    [0, 0, 0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value)]
                ])
        elif symmetry == "m-3m":  # cubic
            # print "Symmetry:", self.symmetry_phase_2
            self.__constant_c_tensor_voigt_notation_fixed_phase = \
                np.array([
                    [pars['c_11'].value, pars['c_12'].value, pars['c_12'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_11'].value, pars['c_12'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_12'].value, pars['c_11'].value, 0, 0, 0],
                    [0, 0, 0, pars['c_44'].value, 0, 0],
                    [0, 0, 0, 0, pars['c_44'].value, 0],
                    [0, 0, 0, 0, 0, pars['c_44'].value]
                ])
        elif symmetry == "hexagonal" or self.symmetry_phase_2 == "hexagonal":
            self.__constant_c_tensor_voigt_notation_fixed_phase = \
                np.array([
                    [pars['c_11'].value, pars['c_12'].value, pars['c_13'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_11'].value, pars['c_13'].value, 0, 0, 0],
                    [pars['c_13'].value, pars['c_13'].value, pars['c33'].value, 0, 0, 0],
                    [0, 0, 0, pars['c_44'].value, 0, 0],
                    [0, 0, 0, 0, pars['c_44'].value, 0],
                    [0, 0, 0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value)]
                ])

        self.__constant_c_tensor_extended_notation_fixed_phase = self.__conv_voigtnot_to_extended_not_constants_c(
            self.__constant_c_tensor_voigt_notation_fixed_phase)
        self.__complience_s_tensor_voigt_notation_fixed_phase = np.linalg.inv(
            self.__constant_c_tensor_voigt_notation_fixed_phase)
        self.__complience_s_tensor_extended_notation_fixed_phase = self.__conv_voigtnot_to_extended_not_compliences_s(
            self.__complience_s_tensor_voigt_notation_fixed_phase)

    def set_parameters_of_the_FIXED_PHASE_in_matrix_representation_single_phase(self, params, symmetry):
        pars = self.__return_free_or_fixed_parameters(params=params, free=True)
        if symmetry == "isotope":
            self.__constant_c_tensor_voigt_notation_fixed_phase = \
                np.array([
                    [pars['c_11'].value, pars['c_12'].value, pars['c_11'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_11'].value, pars['c_12'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_12'].value, pars['c_11'].value, 0, 0, 0],
                    [0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value), 0, 0],
                    [0, 0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value), 0],
                    [0, 0, 0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value)]
                ])
        elif symmetry == "m-3m":  # cubic
            # print "Symmetry:", self.symmetry_phase_2
            self.__constant_c_tensor_voigt_notation_fixed_phase = \
                np.array([
                    [pars['c_11'].value, pars['c_12'].value, pars['c_12'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_11'].value, pars['c_12'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_12'].value, pars['c_11'].value, 0, 0, 0],
                    [0, 0, 0, pars['c_44'].value, 0, 0],
                    [0, 0, 0, 0, pars['c_44'].value, 0],
                    [0, 0, 0, 0, 0, pars['c_44'].value]
                ])
        elif symmetry == "hexagonal" or self.symmetry_phase_2 == "hexagonal":
            self.__constant_c_tensor_voigt_notation_fixed_phase = \
                np.array([
                    [pars['c_11'].value, pars['c_12'].value, pars['c_13'].value, 0, 0, 0],
                    [pars['c_12'].value, pars['c_11'].value, pars['c_13'].value, 0, 0, 0],
                    [pars['c_13'].value, pars['c_13'].value, pars['c33'].value, 0, 0, 0],
                    [0, 0, 0, pars['c_44'].value, 0, 0],
                    [0, 0, 0, 0, pars['c_44'].value, 0],
                    [0, 0, 0, 0, 0, 2 * (pars['c_11'].value - pars['c_12'].value)]
                ])

        self.__constant_c_tensor_extended_notation_fixed_phase = self.__conv_voigtnot_to_extended_not_constants_c(
            self.__constant_c_tensor_voigt_notation_fixed_phase)
        self.__complience_s_tensor_voigt_notation_fixed_phase = np.linalg.inv(
            self.__constant_c_tensor_voigt_notation_fixed_phase)
        self.__complience_s_tensor_extended_notation_fixed_phase = self.__conv_voigtnot_to_extended_not_compliences_s(
            self.__complience_s_tensor_voigt_notation_fixed_phase)

    def __residuum_with_texture(self, params, xvals, data_phase_1=None, data_phase_2=None, fitted_phase=1,
                                method=None):
        """
        fit the params of the matrix and keep the one of the inclusion fixed.
        @ eshelby, matrix is fixed inclusion is fitted
        :type method: Theoretical model to fit the data
        :param params: lm.Parameter Object
        :param xvals: list of the xvalues [[phi, psi, h, k, l], [phi, psi, h, k, l], ...
        :return:
        """
        self.__counter += 1
        self.counter2 = 0
        symmetry_fixed_phase = 0
        symmetry_fitted_phase = 0
        print "Iteration #", self.__counter
        if fitted_phase == 1:
            symmetry_fitted_phase = self.symmetry_phase_1
            if self.phase_flag:
                symmetry_fixed_phase = self.symmetry_phase_2

        elif fitted_phase == 2:
            symmetry_fitted_phase = self.symmetry_phase_2
            symmetry_fixed_phase = self.symmetry_phase_1

        self.set_parameters_of_the_FITTED_PHASE_in_matrix_representation(params, symmetry_fitted_phase)
        print "phaseflag ", self.phase_flag
        if self.phase_flag:
            self.set_parameters_of_the_FIXED_PHASE_in_matrix_representation_multi_phase(params, symmetry_fixed_phase)
        else:
            self.set_parameters_of_the_FIXED_PHASE_in_matrix_representation_single_phase(params, symmetry_fitted_phase)

        print "parameter vals:"
        print "C_11_p1: ", params["c_11_p1"].value
        print "C_12_p1: ", params["c_12_p1"].value
        print "C_44_p1: ", params["c_44_p1"].value
        if self.phase_flag:
            print "C_11_p2: ", params["c_11_p2"].value
            print "C_12_p2: ", params["c_12_p2"].value
            print "C_44_p2: ", params["c_44_p2"].value

        # self.__constant_c_tensor_extended_notation_fitted_phase = self.__conv_voigtnot_to_extended_not_constants_c(
        #     self.__constant_c_tensor_voigt_notation_fitted_phase)
        # self.__complience_s_tensor_voigt_notation_fitted_phase = np.linalg.inv(self.__constant_c_tensor_voigt_notation_fitted_phase)
        # # self.__conv_all_voigtnot_to_extended_not()
        # self.__complience_s_tensor_extended_notation_fitted_phase = self.__conv_voigtnot_to_extended_not_compliences_s(
        #     self.__complience_s_tensor_voigt_notation_fitted_phase)

        # self.__complience_s_tensor_extended = np.linalg.inv(self.__constant_c_tensor_extended)

        # print "tensor identity calc:\n",np.tensordot(self.__constant_c_tensor_extended, self.__complience_s_tensor_extended)
        # print "tensor identity teo:\n",self.fourth_rank_identity()
        # print "tensor diference:\n", np.tensordot(self.__constant_c_tensor_extended_notation_fitted_phase,
        #                                           self.__complience_s_tensor_extended_notation_fitted_phase) - \
        #                              self.fourth_rank_identity()
        # strain_epsilon_2 = []
        # strain_epsilon = 0.
        t1 = tm.clock()
        co = 0
        if method == "voigt" or method == "hill":
            euler = (0, 0, 0)
            self.a_voigt = self.__A_voigt(euler, 0, 0, 0, 0, fitted_phase=fitted_phase)

        xvals_fitted = []
        xvals_fixed = []
        strain_stress_data_fitted = []
        strain_stress_data_fixed = []

        print data_phase_1
        applied_forces = data_phase_1.keys()

        for n in xrange(len(applied_forces)):  # Loop over all forces
            force = applied_forces[n]
            # list of the xvalues [[phi, psi, h, k, l], [phi, psi, h, k, l], ...
            xvals_fitted_, strain_stre_data_fitted_ = 0, 0
            if fitted_phase == 1:
                xvals_fitted_, strain_stre_data_fitted_ = data_phase_1[force]
            if fitted_phase == 2:
                xvals_fitted_, strain_stre_data_fitted_ = data_phase_2[force]
            strain_stress_data_fitted.append(strain_stre_data_fitted_)  # list containing the data for
            #  the different applied forces
            xvals_fitted.append(xvals_fitted_)

            if self.phase_flag:
                xvals_fixed_, strain_stre_data_fixed = 0, 0
                if fitted_phase == 1:
                    xvals_fixed_, strain_stre_data_fixed = data_phase_2[force]
                if fitted_phase == 2:
                    xvals_fixed_, strain_stre_data_fixed = data_phase_1[force]

                strain_stress_data_fixed.append(strain_stre_data_fixed)
                xvals_fixed.append(xvals_fixed_)

        data_fix = []
        data_fit = []
        data_fit_err = []
        data_fix_err = []
        theory_fit = []
        theory_fix = []

        for n in xrange(len(applied_forces)):  # Loop over all forces
            xvals_fit = xvals_fitted[n]
            # print xvals_mat
            # print "applied_force: ", applied_forces
            # xvals_inc = xvals_inclusion[n]
            # print len(xvals_mat)
            for m in xrange(len(xvals_fit)):
                phi, psi, h, k, l = xvals_fit[m]
                # print "phi: ", phi
                strain_fit, strain_fit_err, stress_fit, stress_fit_err = strain_stress_data_fitted[n][m]
                data_fit.append(strain_fit / stress_fit)
                # data_mat.append(strain_mat)
                data_fit_err.append(abs(strain_fit / stress_fit) * (abs(strain_fit_err / strain_fit) +
                                                                    abs(stress_fit_err / stress_fit)))
                # data_mat_err.append(strain_mat_err)

                theory_val = 0.
                if method == "eshelby":
                    for i in xrange(3):
                        for j in xrange(3):
                            theory_val += self.F(phi, psi, h, k, l, i, j, method) * self.force_factor(i, j)
                else:
                    theory_val += self.F(phi, psi, h, k, l, 2, 2, method)

                theory_fit.append(theory_val)
                cli_progress_test(co, len(xvals_fitted[n]))
                co += 1

        t2 = tm.clock()
        dt = t2 - t1
        print "\nresult: ", (np.array(theory_fit) - np.array(data_fit)) / (np.array(data_fit_err))
        print "theory: ", theory_fit
        print "data: ", data_fit
        print "err: ", data_fit_err
        print "shape: ", np.shape((np.array(theory_fit) - np.array(data_fit)) / (np.array(data_fit_err)))
        print "time for iteration #%i: %i min %i sec" % (self.__counter, int(dt / 60), int(dt % 60))

        return (np.array(theory_fit) - np.array(data_fit)) / (np.array(data_fit_err))

    def __residuum_without_texture_single_phase(self, params, xvals, data_phase_1=None, data_phase_2=None,
                                                fitted_phase=1, method=None):
        """
        :param params: lm.Parameter Object
        :param xvals: list of the xvalues [[phi, psi, h, k, l], [phi, psi, h, k, l], ...
        :param data: the data
        :param weight: the error's of the data
        :return:
        """
        self.__counter += 1
        self.counter2 = 0
        symmetry_fixed_phase = 0
        symmetry_fitted_phase = 0
        print "Iteration #", self.__counter

        if fitted_phase == 1:
            symmetry_fitted_phase = self.symmetry_phase_1
            if self.phase_flag:
                symmetry_fixed_phase = self.symmetry_phase_2

        elif fitted_phase == 2:
            symmetry_fitted_phase = self.symmetry_phase_2
            symmetry_fixed_phase = self.symmetry_phase_1

        self.set_parameters_of_the_FITTED_PHASE_in_matrix_representation(params, symmetry_fitted_phase)

        if self.phase_flag:
            self.set_parameters_of_the_FIXED_PHASE_in_matrix_representation_multi_phase(params, symmetry_fixed_phase)
        else:
            self.set_parameters_of_the_FIXED_PHASE_in_matrix_representation_single_phase(params, symmetry_fitted_phase)

        print "parameter vals:"
        params.pretty_print()

        xvals_fitted = []
        xvals_fixed = []
        strain_stress_data_fitted = []
        strain_stress_data_fixed = []

        applied_forces = data_phase_1.keys()

        for n in xrange(len(applied_forces)):  # Loop over all forces
            force = applied_forces[n]
            # list of the xvalues [[phi, psi, h, k, l], [phi, psi, h, k, l], ...
            xvals_fitted_, strain_stre_data_fitted_ = 0, 0
            if fitted_phase == 1:
                xvals_fitted_, strain_stre_data_fitted_ = data_phase_1[force]
            if fitted_phase == 2:
                xvals_fitted_, strain_stre_data_fitted_ = data_phase_2[force]
            strain_stress_data_fitted.append(strain_stre_data_fitted_)  # list containing the data for
            #  the different applied forces
            xvals_fitted.append(xvals_fitted_)
            if self.phase_flag:
                xvals_fixed_, strain_stre_data_fixed = 0, 0
                if fitted_phase == 1:
                    xvals_fixed_, strain_stre_data_fixed = data_phase_2[force]
                if fitted_phase == 2:
                    xvals_fixed_, strain_stre_data_fixed = data_phase_1[force]

                strain_stress_data_fixed.append(strain_stre_data_fixed)
                xvals_fixed.append(xvals_fixed_)
        data_fix = []
        data_fit = []
        data_fit_err = []
        data_fix_err = []
        theory_fix = []

        theory_fit = []
        t1 = tm.clock()
        pars = self.__return_free_or_fixed_parameters(params=params, free=True)
        for n in xrange(len(applied_forces)):  # Loop over all forces
            xvals_fit = xvals_fitted[n]
            for m in xrange(len(xvals_fit)):
                phi, psi, h, k, l = xvals_fit[m]
                strain_fit, strain_fit_err, stress_fit, stress_fit_err = strain_stress_data_fitted[n][m]
                data_fit.append(strain_fit)
                data_fit_err.append(strain_fit_err)
                sigma_11 = params["sigma_11"].value * stress_fit
                sigma_22 = params["sigma_22"].value * stress_fit
                sigma_33 = params["sigma_33"].value * stress_fit
                sigma_12 = params["sigma_12"].value * stress_fit
                sigma_13 = params["sigma_13"].value * stress_fit
                sigma_23 = params["sigma_23"].value * stress_fit

                if method == "hill":
                    s1, s2 = RV(Gamma=Gama(h, k, l), c_11=pars["c_11"].value, c_12=pars["c_12"].value,
                                c_44=pars["c_44"].value)
                if method == "voigt":
                    s1, s2 = Voigt__(Gamma=Gama(h, k, l), c_11=pars["c_11"].value, c_12=pars["c_12"].value,
                                     c_44=pars["c_44"].value)
                if method == "eshelby":
                    s1, s2 = BHM_dW(Gamma=Gama(h, k, l), c_11=pars["c_11"].value, c_12=pars["c_12"].value,
                                    c_44=pars["c_44"].value)
                if method == "reus":
                    s1, s2 = Reus(Gamma=Gama(h, k, l), c_11=pars["c_11"].value, c_12=pars["c_12"].value,
                                  c_44=pars["c_44"].value)

                eps = s1 * (sigma_11 + sigma_22 + sigma_33) \
                      + s2 * (sigma_11 * np.cos(phi) ** 2 * np.sin(psi) ** 2 +
                              sigma_22 * np.sin(phi) ** 2 * np.sin(psi) ** 2 +
                              sigma_33 * np.cos(psi) ** 2) \
                      + s2 * (sigma_12 * np.sin(2 * phi) * np.sin(psi) ** 2 +
                              sigma_13 * np.cos(phi) * np.sin(2 * psi) +
                              sigma_23 * np.sin(phi) * np.sin(2 * psi))

                theory_fit.append(eps)

        t2 = tm.clock()
        dt = t2 - t1
        # print "time for iteration #%i: %i min %i sec" % (self.__counter, int(dt / 60), int(dt % 60))

        # if data is None and weight is None:
        #     return theory_fit
        #
        # if weight is None:
        #     return np.array(data) - np.array(theory_fit)

        return (np.array(data_fit) - np.array(theory_fit)) / (np.array(data_fit_err))

    @staticmethod
    def __return_free_or_fixed_parameters(params, free=True):
        """
        just returning a lm.Parameters() object containing all free/fixed parameters
        :param params:
        :return:
        """
        key_words = ["c_11", "c_12", "c_13", "c_14", "c_15", "c_16",
                     "c_22", "c_23", "c_24", "c_25", "c_26",
                     "c_33", "c_34", "c_35", "c_36",
                     "c_44", "c_45", "c_46",
                     "c_55", "c_56",
                     "c_66"]
        pars = lm.Parameters()
        for key in params.keys():
            condition = params[key].vary
            if not free:
                condition = not params[key].vary
            if condition:
                for k in key_words:
                    if k in key:
                        pars.add(k, params[key].value)
        return pars

    def do_the_fitting(self, filename, material, method="reus", path=".\\results\\", texture=False, phase=1,
                       phase_name=""):
        """
        :param phase: phase to fit
        :param filename: name of the outputfile
        :param material: name of the material
        :param method: name of the used fitting method
        :param path: path to the output file (default: .\\results\\)
        :param texture: use the texture for the fit: Yes-> True, No-> False
        :return: result of the fit
        """
        self.__counter = 0
        params = self.__params
        data = []  # self.__strains_data
        weight = []  # self.__weights
        xvals = []  # self.xvals
        t1 = tm.clock()
        date = tm.localtime()
        fit_method = 'leastsq'  # the optons are:
        # leastsq, nelder, lbfgsb, powell, cg, newton, cobyla, tnc, dogleg, slsqp,
        # differential_evolution
        data_phase_1 = self.data_object.fitted_data.get_force_dict_phase_1()
        data_phase_2 = self.data_object.fitted_data.get_force_dict_phase_2()
        print "Texture: ", texture
        params_keys = params.keys()
        for key in params_keys:
            if phase == 1:
                if "p1" in key:
                    params[key].vary = True
                if "p2" in key:
                    params[key].vary = False
            if phase == 2:
                if "p1" in key:
                    params[key].vary = False
                if "p2" in key:
                    params[key].vary = True
        if texture:
            result = lm.minimize(self.__residuum_without_texture_single_phase, params, method=fit_method, args=(xvals,),
                                 kws={'data_phase_1': data_phase_1, 'data_phase_2': data_phase_2, 'method': method,
                                      'fitted_phase': phase})
            params = result.params
            print "Params withot texture: "
            print lm.fit_report(result.params)
            self.__counter = 0
            result = lm.minimize(self.__residuum_with_texture, params, method=fit_method, args=(xvals,),
                                 kws={'data_phase_1': data_phase_1, 'data_phase_2': data_phase_2, 'method': method,
                                      'fitted_phase': phase})
        else:
            result = lm.minimize(self.__residuum_without_texture_single_phase, params, method=fit_method, args=(xvals,),
                                 kws={'data_phase_1': data_phase_1, 'data_phase_2': data_phase_2, 'method': method,
                                      'fitted_phase': phase})
        t2 = tm.clock()
        dt = t2 - t1
        print "time for fit: ", dt
        nice_result = self.__print_result_nicely(result, fitting_time=dt, date_of_fit=date, method=fit_method)

        filename = path + filename
        filename = self.__save_data(filename, material, phase_name, nice_result)
        print "over all compliance Matrix:"
        print self.calc_the_over_all_constant_c_tensor(self.__constant_c_tensor_extended_notation_fitted_phase,
                                                       fitted_phase=phase)
        return result, filename

    def do_the_fitting_self_consistent_sigma_and_el_const(self, filename, material, method="reus", path=".\\results\\",
                                                          texture=False, phase=1, phase_name=""):
        self.__counter = 0
        params = self.__params
        data = self.__strains_data
        weight = self.__weights
        xvals = self.xvals
        data_phase_1 = self.data_object.fitted_data.get_force_dict_phase_1()
        data_phase_2 = self.data_object.fitted_data.get_force_dict_phase_2()
        t1 = tm.clock()
        date = tm.localtime()
        fit_method = 'leastsq'  # the optons are:
        # leastsq, nelder, lbfgsb, powell, cg, newton, cobyla, tnc, dogleg, slsqp,
        # differential_evolution
        # calc elastic constants with given sigma as a start value
        result = lm.minimize(self.__residuum_without_texture_single_phase, params, method=fit_method, args=(xvals,),
                             kws={'data_phase_1': data_phase_1, 'data_phase_2': data_phase_2, 'method': method,
                                  'fitted_phase': phase})
        params_new = result.params
        # if this token is true, fit sigma else fit the elastic constants
        sigma_constant_token = False
        delta = []
        upper_limit = np.power(10., 2)
        while True:
            if sigma_constant_token:
                params["c_11"].vary = False
                params["c_12"].vary = False
                params["c_44"].vary = False
                params["sigma_11"].vary = True
                params["sigma_22"].vary = True
                params["sigma_33"].vary = True
                params["sigma_12"].vary = True
                params["sigma_13"].vary = True
                params["sigma_23"].vary = True
                sigma_constant_token = False
                result = lm.minimize(self.__residuum_without_texture_single_phase, params, method=fit_method,
                                     args=(xvals,),
                                     kws={'data': data, 'weight': weight, 'method': "eshelby"})
            else:
                params["c_11"].vary = True
                params["c_12"].vary = True
                params["c_44"].vary = True
                params["sigma_11"].vary = False
                params["sigma_22"].vary = False
                params["sigma_33"].vary = False
                params["sigma_12"].vary = False
                params["sigma_13"].vary = False
                params["sigma_23"].vary = False
                sigma_constant_token = True
                result = lm.minimize(self.__residuum_without_texture_single_phase, params, method=fit_method,
                                     args=(xvals,),
                                     kws={'data': data, 'weight': weight, 'method': "eshelby"})
            params_new = result.params
            # now compare the old and the new result:
            delta1 = abs(params_new["c_11"].value - params["c_11"].value)
            delta2 = abs(params_new["c_12"].value - params["c_12"].value)
            delta3 = abs(params_new["c_44"].value - params["c_44"].value)
            delta4 = abs(params_new["sigma_11"].value - params["sigma_11"].value)
            delta5 = abs(params_new["sigma_22"].value - params["sigma_22"].value)
            delta6 = abs(params_new["sigma_33"].value - params["sigma_33"].value)
            delta7 = abs(params_new["sigma_12"].value - params["sigma_12"].value)
            delta8 = abs(params_new["sigma_13"].value - params["sigma_13"].value)
            delta9 = abs(params_new["sigma_23"].value - params["sigma_23"].value)
            delta_sum = delta1 + delta2 + delta3 + delta4 + delta5 + delta6 + delta7 + delta8 + delta9

            delta.append(delta_sum)
            print delta_sum, len(delta)

            if delta_sum < upper_limit:
                break

            if len(delta) > 1000000:
                upper_limit = np.min(np.array(delta))
            params = params_new
        params["c_11"].vary = True
        params["c_12"].vary = True
        params["c_44"].vary = True
        params["sigma_11"].vary = False
        params["sigma_22"].vary = False
        params["sigma_33"].vary = False
        params["sigma_12"].vary = False
        params["sigma_13"].vary = False
        params["sigma_23"].vary = False
        result = lm.minimize(self.__residuum_without_texture_single_phase, params, method=fit_method, args=(xvals,),
                             kws={'data': data, 'weight': weight, 'method': "eshelby"})
        # if texture:
        #     result = lm.minimize(self.__residuum_without_texture, params, method=fit_method, args=(xvals,),
        #                          kws={'data': data, 'weight': weight, 'method': "eshelby"})
        #     params = result.params
        #     self.__counter = 0
        #     result = lm.minimize(self.__residuum_with_texture, params, method=fit_method, args=(xvals,),
        #                          kws={'data': data, 'weight': weight, 'method': method})
        # else:
        #     result = lm.minimize(self.__residuum_without_texture, params, method=fit_method, args=(xvals,),
        #                          kws={'data': data, 'weight': weight, 'method': method})
        t2 = tm.clock()
        dt = t2 - t1
        print "time for fit: ", dt
        print "delta_val: ", delta_sum
        nice_result = self.__print_result_nicely(result, fitting_time=dt, date_of_fit=date, method=fit_method)
        Material = material
        filename = path + filename
        self.__save_data(filename, Material, nice_result)
        return result

    def __print_result_nicely(self, res, **kwargs):
        """
        returns a string with the results and some other information about the fit.
        :param res:
        :param kwargs:
        :return:
        """
        fit_time = kwargs["fitting_time"]
        h = int(fit_time / 3600)
        m = int((fit_time % 3600) / 60)
        s = int(fit_time % 60)
        comment = "selfconsistent calculation of sigma tensor and elastic constants"
        # "Using an other definition for the reus model (wreite s(g) =(g_im * g_jn * g_ko * g_lp * c^0_mnop)^-1"  # add some comment here
        time = "%ih %i min %i sec" % (h, m, s)
        date = kwargs["date_of_fit"]
        da = tm.strftime("%d.%m.%Y, %H:%M", date)
        pars = lm.fit_report(res.params)
        sym_matrix = self.symmetry_phase_1 + " " + self.__odf_phase_1.crystal_symmetry
        sym_inclusion = "None"
        if self.phase_flag:
            sym_inclusion = self.symmetry_phase_2 + " " + self.__odf_phase_2.crystal_symmetry
        out = \
            "\nMethod:         %s\
           \nSymmetry:       %s\
           \nDate:           %s\
           \nFitting time    %s\
           \nComment:        %s\
           \nMessage:        %s\
           \nCovar matrix:   \n%s\
           \nSuccess:        %s\
           \nChisqr:         %f\
           \nReducec Chisqr: %f\
           \ndeg of freedom: %i\
           \n# of datapoint: %i\
           \nnfev:           %i\
           \nParameters:\n%s " \
            % (kwargs["method"], sym_matrix, da, time, comment, res.message, str(res.covar), res.success, res.chisqr,
               res.redchi, res.nfree, res.ndata, res.nfev, pars)
        return out

    def __test_if_file_exists(self, filename):
        filename = filename
        if os.path.isfile(filename):
            try:
                index = int(filename[-5])
                filename = filename[0:-5] + str(index + 1) + filename[-4:]
            except ValueError:
                filename = filename[0:-4] + str(1) + filename[-4:]
            return self.__test_if_file_exists(filename)
        else:
            return filename

    @staticmethod
    def __uniquify(path, sep=''):
        def name_sequence():
            count = IT.count()
            yield ''
            while True:
                yield '{s}{n:d}'.format(s=sep, n=next(count))

        orig = tempfile._name_sequence
        with tempfile._once_lock:
            tempfile._name_sequence = name_sequence()
            path = os.path.normpath(path)
            dirname, basename = os.path.split(path)
            filename, ext = os.path.splitext(basename)
            fd, filename = tempfile.mkstemp(dir=dirname, prefix=filename, suffix=ext)
            tempfile._name_sequence = orig
        return filename

    def __save_data(self, filename, material, phase, data):
        filename = u'{0:s}.txt'.format(filename)
        # filename = self.__test_if_file_exists(filename)
        filename = self.__uniquify(filename)
        result = open(filename, "w")
        string_to_write = "Material: %s\nFitted phase: %s\n" % (material, phase) + data
        result.write(string_to_write)
        result.close()
        return filename

    def invert_four_rank_c_tensor(self, c_tensor):
        c_voigt = self.__conv_extended_not_to_voigt_not_constants_c(c_tensor)
        s_voigt = np.linalg.inv(c_voigt)
        s_extend = self.__conv_voigtnot_to_extended_not_compliences_s(s_voigt)
        return s_extend

    def invert_four_rank_s_tensor(self, s_tensor):
        s_voigt = self.__conv_voigtnot_to_extended_not_compliences_s(s_tensor)
        c_voigt = np.linalg.inv(s_voigt)
        c_extend = self.__conv_extended_not_to_voigt_not_constants_c(c_voigt)
        return c_extend

    @staticmethod
    def tensor_dotprouct(a, b):
        return np.tensordot(a, b)

    @property
    def get_compliance_tensor_voigt_not(self):
        return self.__constant_c_tensor_voigt_notation_fitted_phase

    def kronneker_deltha(self, i, j):
        if i == j:
            return 1
        else:
            return 0

    def fourth_rank_identity(self):
        res = np.zeros((3, 3, 3, 3))
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        res[i, j, k, l] = 0.5 * (self.kronneker_deltha(i, k) * self.kronneker_deltha(j, l) + \
                                                 self.kronneker_deltha(i, l) * self.kronneker_deltha(j, k))
        return res

    def F(self, phi, psi, h, k, l, i, j, method, fitted_phase=1, use_in_fit=True):
        """
        calculates the force factors (dt: Spannungsfaktoren) for the selected phase using the model defined in method
        :param method:
        :param phi:
        :param psi:
        :param h:
        :param k:
        :param l:
        :param i:
        :param j:
        :rtype: object
        """
        ODF = 0
        if fitted_phase == 1:
            ODF = self.__odf_phase_1
        elif fitted_phase == 2:
            ODF = self.__odf_phase_2
        if not use_in_fit:
            self.set_parameters_of_the_FITTED_PHASE_in_matrix_representation(self.__params, ODF.crystal_symmetry)

        self.counter2 += 1
        res = 0
        if method == "reus":
            for u in xrange(3):
                for w in xrange(3):
                    res += ODF.integrate(self.A_reus, phi, psi, h, k, l, u, w, i, j) / \
                           ODF.integrate_(phi, psi, h, k, l)

        elif method == "voigt":
            for u in xrange(3):
                for w in xrange(3):
                    res += ODF.integrate(self.__A_voigt_call, phi, psi, h, k, l, u, w, i, j) / \
                           ODF.integrate_(phi, psi, h, k, l)

        elif method == "hill":
            for u in xrange(3):
                for w in xrange(3):
                    res += ((ODF.integrate(self.A_reus, phi, psi, h, k, l, u, w, i, j) +
                             ODF.integrate(self.__A_voigt_call, phi, psi, h, k, l, u, w, i, j)) / 2) / \
                           ODF.integrate_(phi, psi, h, k, l)

        elif method == "eshelby":
            for u in xrange(3):
                for w in xrange(3):
                    res += ODF.integrate(self.A_eshelby, phi, psi, h, k, l, u, w, i, j) / \
                           ODF.integrate_(phi, psi, h, k, l)
        # print "Func. calls: ", self.counter2, "Spannungsfaktor: ", res, " phi: ", rad_to_deg(phi), " psi: ", \
        #     rad_to_deg(psi), "hkl: ", h, k, l
        return res

    def stress_sigma(self, i, j):
        """
            seting up the stress tensor, the applied stress is draged out of the tensor
            This function determines the 33-component of the macro-straintensor due


            to the Force in this (and only in this) direction with respect to
            the specimen reference frame.
            The component is:
                sigma33 = F/A
            where A is the crosssection of the Specimen
            :param j:
            :param i:
        """
        sig = np.zeros((3, 3))
        nu = 0.3
        sig[2, 2] = 1
        # sig[1, 1] = - nu * sigma3
        # sig[0, 0] = - nu * sigma3
        # sig = add_rot().dot(sig.dot(add_rot().transpose()))
        # print sig
        # print "sigma_33: ", sigma3, self.force, self.diameter
        return sig[i, j]

    def stress_sigma_in_L_frame(self, i, j, psi, phi, phi2=0):
        """
            This function determines the 33-component of the macro-straintensor due
            to the Force in this (and only in this) direction with respekt to
            the specimen reference frame.
            The component is:
                sigma33 = F/A
            where A is the crosssection of the Specimen
            :param j:
            :param i:
        """
        A = (self.diameter * np.power(10., -3.)) ** 2 / 4 * np.pi
        sigma3 = self.force * np.power(10., 3) / A  # *np.power(10.,-9)
        sig = np.zeros((3, 3))
        # sig[0,0]=np.sin(psi)**2 * np.sin(phi2)**2
        # sig[2, 2] = sigma3
        g = self.__odf_phase_1.g1(phi, psi, phi2)
        for i in xrange(3):
            for j in xrange(3):
                sig[i, j] = g[i, 2] * g[j, 2] * 1
        # print "sigma_33: ", sigma3, self.force, self.diameter
        return sig[i, j]

    def __over_all_inner_sum(self, phi1, phi, phi2, a, b, i, j, c_Matrix_tensor):
        res = 0
        g = self.__odf_phase_1.g(phi1, phi, phi2).transpose()  # this is the transformation tensor it is identical
        # for both phases
        g = g.tolist()
        c_Matrix_tensor = c_Matrix_tensor.tolist()
        for m in xrange(3):
            for n in xrange(3):
                for o in xrange(3):
                    for p in xrange(3):
                        res += g[a][m] * g[b][n] * g[i][o] * g[j][p] \
                               * c_Matrix_tensor[m][n][o][p]
        return res

    def calc_the_over_all_constant_c_tensor(self, compliance_tensor, fitted_phase):
        """
        calculate the over all constant tensor of the fitted phase
        :param compliance_tensor:
        :param fitted_phase:
        :return:
        """
        over_all_tensor = np.zeros((3, 3, 3, 3))  # compliance tensor g independent
        over_all_tensor = over_all_tensor.tolist()
        cout = 0
        integral = 0
        for a in xrange(3):
            for b in xrange(3):
                for f in xrange(3):
                    for d in xrange(3):
                        cli_progress_test_voigt(cout, 81, (a, b, f, d))
                        if fitted_phase == 1:
                            over_all_tensor[a][b][f][d] = self.__odf_phase_1.integrate_a_over_all_orientations_g(
                                self.__over_all_inner_sum, a, b, f, d, compliance_tensor)
                            integral = self.odf_integral_over_all_angles_phase_1
                        elif fitted_phase == 2:
                            over_all_tensor[a][b][f][d] = self.__odf_phase_2.integrate_a_over_all_orientations_g(
                                self.__over_all_inner_sum, a, b, f, d, compliance_tensor)
                            integral = self.odf_integral_over_all_angles_phase_2
                        cout += 1
        over_all_tensor = np.array(over_all_tensor)
        over_all_tensor /= integral
        over_all_tensor = self.__conv_extended_not_to_voigt_not_constants_c(over_all_tensor)

        return over_all_tensor

    def A_reus(self, g2, g1, u, w, i, j):
        """
        A(g) = s(g)
        :param euler:
        :param u:
        :param w:
        :param i:
        :param j:
        :return:
        """
        # phi1, phi, phi2 = euler
        g2 = g2.transpose()

        res = 0.
        for m in xrange(3):
            for n in xrange(3):
                for o in xrange(3):
                    for p in xrange(3):
                        res += g2[i, m] * g2[j, n] * g2[u, o] * g2[w, p] * \
                               self.__complience_s_tensor_extended_notation_fitted_phase[m, n, o, p]

        # res = np.zeros((3, 3, 3, 3))
        # for a in xrange(3):
        #     for b in xrange(3):
        #         for c in xrange(3):
        #             for d in xrange(3):
        #                 for m in xrange(3):
        #                     for n in xrange(3):
        #                         for o in xrange(3):
        #                             for p in xrange(3):
        #                                 res[a, b, c, d] += g[a, m] * g[b, n] * g[c, o] * g[d, p] * \
        #                                                    self.__constant_c_tensor_extended_notation_fitted_phase[m, n, o, p]
        # print res
        # print self.__constant_c_tensor_extended_notation_fitted_phase
        # print self.__constant_c_tensor_voigt_notation_fitted_phase

        # res = self.invert_four_rank_c_tensor(res)
        return res  # [2, 2, u, w]  # [u, w, i, j]

    def __voigt_inner_sum(self, phi1, phi, phi2, a, b, i, j):
        res = 0
        g = self.__odf_phase_1.g(phi1, phi, phi2).transpose()  # this is just the transformation tensor it is identical
        # for both phases
        g = g.tolist()
        c_Matrix_tensor = self.__constant_c_tensor_extended_notation_fitted_phase.tolist()
        for m in xrange(3):
            for n in xrange(3):
                for o in xrange(3):
                    for p in xrange(3):
                        res += g[a][m] * g[b][n] * g[i][o] * g[j][p] \
                               * c_Matrix_tensor[m][n][o][p]
        return res

    # def Queue(self, q, a, b, i, j):
    #     q.put([a, b, i, j, self.odf_phase_1.integrate_a_over_all_orientations_g(self.__voigt_inner_sum, a, b, i, j)])

    def __calc(self, d_list):
        res = np.zeros((3, 3, 3, 3))

        # print "a: ", d_list["a"], "b: ", d_list["b"], "f: ", d_list["f"], "d: ", d_list["d"]
        res[d_list["a"], d_list["b"], d_list["f"], d_list["d"]] = self.__odf_phase_1. \
            integrate_a_over_all_orientations_g(self.__voigt_inner_sum,
                                                d_list["a"], d_list["b"], d_list["f"], d_list["d"])
        return res

    def __A_voigt(self, g2, u, w, i, j, fitted_phase):
        """
        Model of voigt
        A(g) = <c(g)>^-1
        :param g2:
        :param u:
        :param w:
        :param i:
        :param j:
        :return:
        """
        # phi1, phi, phi2 = g2
        # g = self.odf_phase_1.g(phi1, phi, phi2)
        c = np.zeros((3, 3, 3, 3))  # compliance tensor g independent
        c = c.tolist()
        cout = 0
        integral = 0
        for a in xrange(3):
            for b in xrange(3):
                for f in xrange(3):
                    for d in xrange(3):
                        cli_progress_test_voigt(cout, 81, (a, b, f, d))
                        if fitted_phase == 1:
                            c[a][b][f][d] = self.__odf_phase_1.integrate_a_over_all_orientations_g(
                                self.__voigt_inner_sum, a, b, f, d)
                            integral = self.odf_integral_over_all_angles_phase_1
                        elif fitted_phase == 2:
                            c[a][b][f][d] = self.__odf_phase_2.integrate_a_over_all_orientations_g(
                                self.__voigt_inner_sum, a, b, f, d)
                            integral = self.odf_integral_over_all_angles_phase_2
                        cout += 1
        c = np.array(c)
        c /= integral  # self.odf_integral_over_all_angles_phase_1
        s = self.invert_four_rank_c_tensor(c)
        # s = self.invert_tensor(c)
        return s  # [u, w, i, j]

    def __A_voigt_call(self, g2, g1, u, w, i, j):
        """
        A(g) = <c(g)>^-1
        :param g1:
        :param u:
        :param w:
        :param i:
        :param j:
        :return:
        """
        # S_L = 0.
        # for m in xrange(3):
        #     for n in xrange(3):
        #         for o in xrange(3):
        #             for p in xrange(3):
        #                 S_L += g1[2, m] * g1[2, n] * g1[u, o] * g1[w, p] * \
        #                        self.a_voigt[m, n, o, p]

        return self.a_voigt[2, 2, u, w]  # S_L  # self.a_voigt[u, w, i, j]

    def __A_hill(self, g2, g1, u, w, i, j):
        """
        A(g) = (s_voigt(g) + s_reus(g))/2
        :param u:
        :param w:
        :param i:
        :param j:
        :return:
        """
        hill = (self.a_voigt[u, w, i, j] + self.A_reus(g2, g1, u, w, i, j)) / 2
        return hill

    def c(self, g):
        c = np.zeros((3, 3, 3, 3))
        for u in xrange(3):
            for w in xrange(3):
                for i in xrange(3):
                    for j in xrange(3):
                        for m in xrange(3):
                            for n in xrange(3):
                                for o in xrange(3):
                                    for p in xrange(3):
                                        c[u, w, i, j] += g[u, m] * g[w, n] * g[i, o] * g[j, p] * \
                                                         self.__constant_c_tensor_extended_notation_fitted_phase[
                                                             m, n, o, p]
        return c

    @staticmethod
    def __k(alpha, beta):
        res = np.array([np.sin(alpha) * np.cos(beta),
                        np.sin(alpha) * np.sin(beta),
                        np.cos(alpha)
                        ]
                       )
        return res

    def __D(self, alpha, beta, C_Matrix):
        k = self.__k(alpha, beta)
        res = np.array([[C_Matrix[0, 0, 0, 0] * k[0] ** 2 +
                         C_Matrix[1, 0, 0, 1] * k[1] ** 2 +
                         C_Matrix[2, 0, 0, 2] * k[2] ** 2,

                         C_Matrix[0, 0, 1, 1] * k[0] * k[1] +
                         C_Matrix[1, 0, 1, 0] * k[0] * k[1],

                         C_Matrix[0, 0, 2, 2] * k[0] * k[2] +
                         C_Matrix[2, 0, 2, 0] * k[0] * k[2]
                         ],
                        [C_Matrix[0, 1, 0, 1] * k[0] * k[1] +
                         C_Matrix[1, 1, 0, 0] * k[0] * k[1],

                         C_Matrix[0, 1, 1, 0] * k[0] ** 2 +
                         C_Matrix[1, 1, 1, 1] * k[1] ** 2 +
                         C_Matrix[2, 1, 1, 2] * k[2] ** 2,

                         C_Matrix[1, 1, 2, 2] * k[1] * k[2] +
                         C_Matrix[2, 1, 2, 1] * k[1] * k[2]
                         ],
                        [C_Matrix[0, 2, 0, 2] * k[0] * k[2] +
                         C_Matrix[2, 2, 0, 0] * k[0] * k[2],

                         C_Matrix[1, 2, 1, 2] * k[1] * k[2] +
                         C_Matrix[2, 2, 1, 1] * k[1] * k[2],

                         C_Matrix[0, 2, 2, 0] * k[0] ** 2 +
                         C_Matrix[1, 2, 2, 1] * k[1] ** 2 +
                         C_Matrix[3, 3, 3, 3] * k[2] ** 2
                         ]
                        ]
                       )
        return np.linalg.inv(res)

    def __integrand(self, alpha, beta, *args):
        n, g, j, i, C_Matrix = args
        return np.sin(alpha) * self.__D(alpha, beta, C_Matrix)[n, j] \
               * self.__k(alpha, beta)[g] * self.__k(alpha, beta)[i]

    def __calc_E(self, C_Matrix):
        E = np.zeros((3, 3, 3, 3))
        for n in xrange(3):
            for g in xrange(3):
                for j in xrange(3):
                    for i in xrange(3):
                        E[n, g, j, i] = 1 / (4 * np.pi) * scipy.integrate.nquad(self.__integrand,
                                                                                [[0, np.pi], [0, 2 * np.pi]],
                                                                                args=(n, g, j, i, C_Matrix))
        return E

    def __u(self, phi1, phi, phi2, C):
        E = self.__calc_E(C)
        g = self.__odf_phase_1.g(phi1, phi, phi2).transpose()
        c = self.c(g)  # constant tensor of the inclusion depending on the orientation g
        E_inv = self.invert_four_rank_c_tensor(E)
        h = self.invert_four_rank_c_tensor(c + C + E_inv)
        res = np.zeros((3, 3, 3, 3))
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        for m in xrange(3):
                            for n in xrange(3):
                                res[i, j, k, l] += h[i, j, m, n] * E_inv[m, n, k, l]
        # return w.tensordot(np.tensordot(self.invert_four_rank_c_tensor(c - C), w) + C).tensordot(c - C)
        return res - self.fourth_rank_identity()

    def __integrand_int_u(self, phi1, phi, phi2, C):
        return self.__u(phi1, phi, phi2, C) * self.__odf_phase_1.f(phi1, phi, phi2) * np.sin(phi)

    def __int_u(self, C):
        res = np.zeros((3, 3, 3, 3))
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        for phi1 in xrange(0, 360, 5):
                            for phi in xrange(0, 180, 5):
                                for phi2 in xrange(0, 360, 5):
                                    res[i, j, k, l] += self.__integrand_int_u(phi1, phi, phi2, C)
        return res

    def A_eshelby(self, euler, u, w, i, j):
        """
        A(g) = S + t(g)
        :param euler:
        :param u:
        :param w:
        :param i:
        :param j:
        :return:
        """

        phi1, phi, phi2 = euler
        g = self.__odf_phase_1.g(phi1, phi, phi2).transpose()
        S = np.zeros((3, 3, 3, 3))  # averaged constant tensor of the matrix
        S = self.__complience_s_tensor_extended_notation_fitted_phase  # use as a start value
        C = self.invert_four_rank_s_tensor(S)  # averaged compliance tensor of the matrix
        c = self.c(g)  # constant tensor of the inclusion depending on the orientation g

        E = self.__calc_E(C)
        E_inv = self.invert_four_rank_c_tensor(E)

        # t_g
        return S[u, w, i, j] + t(c, C)[u, w, i, j]

    def set_params(self, params_Matrix, params_Inclusion=None):
        self.__params = params_Matrix

    def force_factor(self, i, j):
        res = 1
        return res


def cli_progress_test_voigt(i, end_val, tuple, bar_length=20):
    """
    just a status Bar
    :param i:
    :param end_val:
    :param bar_length:
    :return:
    """
    # for i in xrange(0, end_val):
    a, b, c, d = tuple
    percent = float(i) / end_val
    hashes = '=' * int(round(percent * bar_length))
    spaces = '_' * (bar_length - len(hashes))
    sys.stdout.write(
        "\rc: ({},{},{},{}), Percent: [{}] {}%".format(a, b, c, d, hashes + spaces, int(round(percent * 100))))
    sys.stdout.flush()


def cli_progress_test(i, end_val, bar_length=20):
    """
    just a status Bar
    :param i:
    :param end_val:
    :param bar_length:
    :return:
    """
    percent = float(i) / end_val
    hashes = '=' * int(round(percent * bar_length))
    spaces = '_' * (bar_length - len(hashes))
    sys.stdout.write(
        "\rPercent: [{}] {}%".format(hashes + spaces, int(round(percent * 100))))
    sys.stdout.flush()


'''
Fiting class, that inherits from the class "Fit_strain_with_texture_single_phase"
'''


class Fit_strain_with_texture_two_phase_material(FitStrainWithTexture):
    def __init__(self, odf_Matrix, odf_Inclusion, force, diameter, strains_Matrix_data, xvals_Matrix,
                 strains_Inclusion_data, xvals_Inclusion, weights_Inclusion=None, weights_Matrix=None):
        FitStrainWithTexture.__init__(odf_Matrix=odf_Matrix, force=force, diameter=diameter,
                                      strains_data=strains_Matrix_data, xvals=xvals_Matrix,
                                      weights=weights_Matrix)

        self.__odf_phase_2 = odf_Inclusion
        self.__strains_data_Inclusion = strains_Inclusion_data
        self.__weights_Inclusion = weights_Inclusion
        self.xvals_Inclusion = xvals_Inclusion
        self.__sym_Inclusion = self.__odf_phase_2.crystal_symmetry
        self.__params_Inclusion = lm.Parameters()
        self.add_params(self.__params_Inclusion, self.__sym_Inclusion)


'''
This clas inherits from "Fit_strain_with_texture_single_phase" class and makes nice plots
as input it is nessessary to define the single krystal elastic constants
'''


class make_some_nice_plots(FitStrainWithTexture):
    def __init__(self, odf_Matrix, odf_Inclusion, force, diameter):
        FitStrainWithTexture.__init__(self, odf_Matrix=odf_Matrix, force=force, diameter=diameter,
                                      strains_data=None, xvals=None, weights=None)
        self.params_matrix = lm.Parameters()
        self.params_inclusion = lm.Parameters()
        params = self.insert_constant_params()
        self.print_params()
        for i in xrange(0, 361, 45):
            self.plot_F_sin_2_psi(h=1, k=1, l=0, phi=deg_to_rad(i), method="reus")
        plt.show()

    def recive_data(self, params_matrix, params_inclusion=None):
        self.params_matrix = params_matrix
        if params_inclusion is not None:
            self.params_inclusion = params_inclusion

    def recive_odf_matrix_path(self, path):
        self.set_odf_matrix(path)

    def insert_constant_params(self):
        # gui.main()
        params_matrix = self.params_matrix
        app = PyQt4.QtGui.QApplication(sys.argv)
        mygui = gui.gui("Instert the elastic constants", params_matrix=params_matrix)
        mygui.show()
        # mygui.connect(self, PyQt4.QtGui.SIGNAL("ODF_matrix_path"), self.set_odf_matrix)
        # mygui.connect(PyQt4.QtGui.SIGNAL("params"), self.recive_data)
        app.exec_()

    def print_params(self):
        print "C_11: ", self.params_matrix['c_11_m'].value
        print "C_12: ", self.params_matrix['c_12_m'].value
        print "C_44: ", self.params_matrix['c_44_m'].value
        self.set_params(params_Matrix=self.params_matrix, params_Inclusion=self.params_inclusion)

    def plot_F_sin_2_psi(self, h, k, l, phi, method="reus"):
        Psi = np.arange(0, np.pi / 2, np.pi / 100.)
        F_11_list = []
        F_22_list = []
        F_33_list = []
        F_12_list = []
        for b in xrange(len(Psi)):
            a = Psi[b]
            F_11_list.append(self.F(phi=phi, psi=a, h=h, k=k, l=l, i=0, j=0, method=method, use_in_fit=False))
            F_22_list.append(self.F(phi=phi, psi=a, h=h, k=k, l=l, i=1, j=1, method=method, use_in_fit=False))
            F_33_list.append(self.F(phi=phi, psi=a, h=h, k=k, l=l, i=2, j=2, method=method, use_in_fit=False))
            F_12_list.append(self.F(phi=phi, psi=a, h=h, k=k, l=l, i=0, j=1, method=method, use_in_fit=False))
            cli_progress_test(b, len(Psi))
        F_11_list = np.array(F_11_list) * np.power(10., 12.)
        F_22_list = np.array(F_22_list) * np.power(10., 12.)
        F_33_list = np.array(F_33_list) * np.power(10., 12.)
        F_12_list = np.array(F_12_list) * np.power(10., 12.)
        print F_11_list

        fig = plt.figure("f_sin^2Psi plot phi=%.1f" % (rad_to_deg(phi)))
        ax1 = fig.add_subplot(1, 1, 1)
        # ax2 = fig.add_subplot(2,1,1)
        ax1.plot(np.sin(Psi) ** 2, F_11_list, 'b-', label='F_11, hkl{}{}{}'.format(h, k, l))
        ax1.plot(np.sin(Psi) ** 2, F_22_list, 'r-', label='F_22, hkl{}{}{}'.format(h, k, l))
        ax1.plot(np.sin(Psi) ** 2, F_33_list, 'g-', label='F_33, hkl{}{}{}'.format(h, k, l))
        ax1.plot(np.sin(Psi) ** 2, F_12_list, 'y-', label='F_12, hkl{}{}{}'.format(h, k, l))
        # ax2.plot(np.sin(Psi) ** 2, F_11_list, 'r-', label='fit')
        ax1.set_title('F_11')

        # ax1.ylim((min(F_11_list), max(F_11_list)))

        # axarr[1].scatter(x_list, av_count_l - y_)
        ax1.legend(loc='upper right', numpoints=1)
        # plt.show()


'''
Odf classdefinition
'''


class ODF(object):
    def __init__(self):
        self.__data = np.zeros((0, 4))
        self.__phi1_max = 0
        self.__Phi_max = 0
        self.__phi2_max = 0
        self.__imax = 0
        self.__imin = 0
        self.__path_to_data = ""
        self.__odf_is_interpolated = False
        # dictionary containing some angles
        # phi, psi = angles of m in specimen frame
        # phi_2 = rotation around m in measurement frame
        # phi_b, betha_b = acimut, polar angle, respectively, of m in the crystal frame
        # phi2_ = rotation around m in the crystal frame
        self.__params = {"psi": 0, "phi": 0, "phi_2": deg_to_rad(0), 'phi_b': 0, 'betha_b': 0,
                         'phi2_': 0}  # "phi_2": deg_to_rad(-90)

    def read_data(self, filename):
        """
        This Function reads the data stored in filename
        """
        filename = os.path.normpath(filename)
        path, filename = os.path.split(filename)
        path = "{0}\\".format(path)
        self.name = filename.split(".")[0]
        self.__path_to_data = path
        filename = path + filename
        data = open(filename, "r")
        lines = data.readlines()
        print "filename: ", filename

        phi1 = []
        Phi = []
        phi2 = []
        Intens = []
        mach = r"\"-*[a-zA-Z0-9]*-*[a-zA-Z0-9]*\""
        p = re.compile(mach)
        counter = 0
        for i in xrange(len(lines)):
            if "%" not in lines[i]:
                break
            counter += 1
        print "counter: ", counter
        self.__data = np.zeros((len(lines) - counter, 4))
        for i in xrange(0, len(lines)):
            line = lines[i]
            if "%" in line:
                if "crystal symmetry:" in line:
                    self.crystal_symmetry = re.findall(p, line)[0]
                    self.crystal_symmetry = self.crystal_symmetry.replace('\"', '')
                if "specimen symmetry:" in line:
                    self.specimen_symmetry = re.findall(p, line)[0]
                    self.specimen_symmetry = self.specimen_symmetry.replace('\"', '')
            else:
                line = lines[i].strip()  # removes spaces at the frond and the end
                l = re.split(r'\s*', line)
                self.__data[i - counter] = np.array(l)
                phi1.append(float(l[0]))
                Phi.append(float(l[1]))
                phi2.append(float(l[2]))
                Intens.append(float(l[3]))
        print self.__data
        if abs(int(phi1[1] - phi1[0])) > 1e-13:
            self.__stepwidth = abs(int(phi1[1] - phi1[0]))
        if abs(int(Phi[1] - Phi[0])) > 1e-13:
            self.__stepwidth = abs(int(Phi[1] - Phi[0]))
        if abs(int(phi2[1] - phi2[0])) > 1e-13:
            self.__stepwidth = abs(int(phi2[1] - phi2[0]))

        self.__phi1_max = int(max(phi1))
        self.__Phi_max = int(max(Phi))
        self.__phi2_max = int(max(phi2))
        self.__imax = max(Intens)
        self.__imin = min(Intens)
        data.close()
        print self.specimen_symmetry, self.crystal_symmetry
        print self.__phi1_max, self.__Phi_max, self.__phi2_max
        print self.__stepwidth, phi1[1]
        self.__transform_to_on_3_dim_np_array()

    def print_data(self):
        print self.__data
        print self.__phi1_max, "  ", self.__Phi_max, "  ", self.__phi2_max

    def f(self, phi1, phi, phi2):
        """
        :param phi1: phi1 in deg
        :param phi: Phi in deg
        :param phi2: phi2 in deg
        :return: value of the ODF at (phi1, Phi, phi2)
        """
        phi1 = int(phi1 / self.__stepwidth)
        phi = int(phi / self.__stepwidth)
        phi2 = int(phi2 / self.__stepwidth)
        res = self.__ODF_data[phi2][phi][phi1]
        return res

    @staticmethod
    def calc_phi_b(h, k, l):
        """
        :param h: Millerindex
        :param k: Millerindex
        :param l: Millerindex
        :return: Polar angle of [h,k,l] in the crystall system
        """
        h = float(h)
        angle = np.arccos(float(l) / np.sqrt(h ** 2 + k ** 2 + l ** 2))
        return angle

    @staticmethod
    def calc_betha_b(h, k, l):
        """
        :param h: Millerindex
        :param k: Millerindex
        :param l: Millerindex
        :return: acimut angle of [h,k,l] in the crystall system
        """
        if h > 0:
            return np.arctan(k / float(h))
        elif h == 0:
            return np.sign(k) * np.pi / 2
        elif h < 0 <= k:
            return np.arctan(k / float(h)) + np.pi
        elif h < 0 and k < 0:
            return np.arctan(k / float(h)) - np.pi
        else:
            raise CalcError

    # implementation of the rotations (after gneupel-Herold)
    @staticmethod
    def g1(phi, psi, phi2):
        """
        rotation from the specimen frame to the measurement frame
        g1.transpose().dot(vector) displays the vector in the measurement frame.
        therefore it is necessary to returne the transposed matrix
        :param phi: Azimuth angel of q
        :param psi: polar angle of q
        :param phi2: arbitrary angle
        :return: rotation matrix
        """
        phi += np.pi / 2

        res = np.array([[-np.cos(psi) * np.cos(phi) * np.sin(phi2) - np.sin(phi) * np.cos(phi2),
                         -np.cos(psi) * np.sin(phi) * np.sin(phi2) + np.cos(phi) * np.cos(phi2),
                         np.sin(psi) * np.sin(phi2)],

                        [-np.cos(psi) * np.cos(phi) * np.cos(phi2) + np.sin(phi) * np.sin(phi2),
                         -np.cos(psi) * np.sin(phi) * np.cos(phi2) - np.cos(phi) * np.sin(phi2),
                         np.sin(psi) * np.cos(phi2)],

                        [np.sin(psi) * np.cos(phi),
                         np.sin(psi) * np.sin(phi),
                         np.cos(psi)]
                        ])
        return res  # np.dot(add_rot(), np.dot(res, add_rot().transpose()))

    @staticmethod
    def omega(phi, psi, phi2):
        """
        rotation from the spesimenframe to the measurmentframe
        :param phi: Azimut angel of q
        :param psi: polar angle of q
        :param phi2: arbitary angle
        :return: raotation matrix
        """
        res = np.array([[np.cos(phi) * np.cos(psi),
                         np.sin(phi) * np.cos(psi),
                         -np.sin(psi)],

                        [-np.sin(phi),
                         np.cos(phi),
                         0.],

                        [np.cos(phi) * np.sin(psi),
                         np.sin(phi) * np.sin(psi),
                         np.cos(psi)]
                        ])
        return res

    @staticmethod
    def g2(phi2_, phi_b, beta_b):
        """
        rotation from the measurement frame to the crystal system
        :param phi2_: rotation angle
        :param phi_b: polar angel of the hkl direction
        :param beta_b: azimuth angle of the hkl direction
        :return:
        """
        beta_b = np.pi / 2 - beta_b
        res = np.array([[np.cos(phi2_) * np.sin(beta_b) - np.sin(phi2_) * np.cos(beta_b) * np.cos(phi_b),
                         np.sin(phi2_) * np.sin(beta_b) + np.cos(phi2_) * np.cos(beta_b) * np.cos(phi_b),
                         np.cos(beta_b) * np.sin(phi_b)],
                        [-np.cos(phi2_) * np.cos(beta_b) - np.sin(phi2_) * np.sin(beta_b) * np.cos(phi_b),
                         -np.sin(phi2_) * np.cos(beta_b) + np.cos(phi2_) * np.sin(beta_b) * np.cos(phi_b),
                         np.sin(beta_b) * np.sin(phi_b)],
                        [np.sin(phi2_) * np.sin(phi_b),
                         -np.cos(phi2_) * np.sin(phi_b),
                         np.cos(phi_b)]
                        ])
        return res

    @staticmethod
    def g(phi1, phi, phi2):
        """
        rotation from specimen frame directly into crystal frame using the Euler angles phi1, phi, phi2
        using the definition from Jan Pospiech
        :param phi1: Euler angles
        :param phi: Euler angles
        :param phi2: Euler angles
        :return: rotation matrix
        """
        res = np.array([[np.cos(phi1) * np.cos(phi2) - np.sin(phi1) * np.sin(phi2) * np.cos(phi),
                         np.sin(phi1) * np.cos(phi2) + np.cos(phi1) * np.sin(phi2) * np.cos(phi),
                         np.sin(phi2) * np.sin(phi)],
                        [-np.cos(phi1) * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.cos(phi),
                         -np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(phi),
                         np.cos(phi2) * np.sin(phi)],
                        [np.sin(phi1) * np.sin(phi),
                         -np.cos(phi1) * np.sin(phi),
                         np.cos(phi)]
                        ])
        return res

    @staticmethod
    def m(phi, psi):
        """
        :param phi:
        :param psi:
        :return: measurement direction in the specimen frame
        """
        return np.array([np.cos(phi) * np.sin(psi), np.sin(phi) * np.sin(psi), np.cos(psi)]).transpose()

    def calc_eulerangles(self, r, h, k, l):
        """
        need the angles phi, psi in rad
        :param r:
        :param h:
        :param k:
        :param l:
        :return: eulerangles in rad
        """
        self.__params['phi2_'] = r
        phi_ = self.__params["phi"]
        psi = self.__params["psi"]
        phi_2 = self.__params['phi_2']

        phi2_ = self.__params['phi2_']
        phi_b = self.__params['phi_b']
        betha_b = self.__params['betha_b']
        # print "::::::::phi2_: ", rad_to_deg(phi2_), " phi_2: ", rad_to_deg(phi_2)
        # rotation = self.xi(phi2_, phi_b, betha_b)#xi(phi2_, phi_b, betha_b)
        # g_ = np.dot(self.g2(phi2_, phi_b, betha_b), self.g1(phi_, psi, phi_2))  # rotation, self.omega(phi_, psi, phi_2)
        # g_ = np.dot(rotation.transpose(),  self.omega(phi_, psi, phi_2))#rotation, self.omega(phi_, psi, phi_2)
        # g_ = g_.transpose()
        g_ = np.dot(self.g1(phi_, psi, phi_2).transpose(), self.g2(phi2_, phi_b, betha_b).transpose())
        # print g_
        phi = np.arccos(g_[2, 2])  # this calculation seams to be correct

        # definition of the specimen frame
        P1 = np.array([[1], [0], [0]])
        P2 = np.array([[0], [1], [0]])
        P3 = np.array([[0], [0], [1]])

        # definition of the unit vectors of the crystal frame represented in the specimen frame
        C1 = np.dot(g_, np.array([[1], [0], [0]]))
        C2 = np.dot(g_, np.array([[0], [1], [0]]))
        C3 = np.dot(g_, np.array([[0], [0], [1]]))

        # print "............................................................\n" \
        #       "C1: \n", C1, "\n............................................................\n" \
        #       "C2: \n", C2, "\n............................................................\n" \
        #       "C3: \n", C3
        # a is the projection of C3 to the 1,2 - plane of the specimen frame
        a = np.array([[C3[0]], [C3[1]], [0]])
        # set default vals for phi1 and phi2
        phi1 = 0
        phi2 = 0

        # calculate phi1:
        if abs(phi) < 1e-13 or abs(phi - np.pi) < 1e-13:
            phi1_temp = np.arccos((P1.transpose().dot(C1)[0, 0]) / (np.linalg.norm(P1) * np.linalg.norm(C1)))
            # print rad_to_deg(phi1_temp), P1.transpose().dot(C1), np.linalg.norm(P1)*np.linalg.norm(C1)
            if abs(C1[1]) < 1e-13 and C1[0] > 0:
                # print "fall 1"
                phi1 = 0
            elif abs(C1[1]) < 1e-13 and C1[0] < 0:
                # print "fall 2"
                phi1 = np.pi
            elif C1[1] >= 1e-13:
                # print "fall 3"
                phi1 = phi1_temp
            elif C1[1] <= -1e-13:
                # print "fall 4"
                phi1 = 2 * np.pi - phi1_temp

        else:
            # print "hallo: ", -P2.transpose().dot(a)[0, 0], np.arccos(-P2.transpose().dot(a)[0, 0])
            phi1_temp = np.arccos((-P2.transpose().dot(a)[0, 0]) / (np.linalg.norm(P2) * np.linalg.norm(a)))
            if abs(C3[0]) < 1e-13 and C3[1] < 0:
                phi1 = 0
            elif abs(C3[0]) < 1e-13 and C3[1] > 0:
                phi1 = np.pi
            elif C3[0] >= 1e-13:
                phi1 = phi1_temp
            elif C3[0] <= -1e-13:
                phi1 = 2 * np.pi - phi1_temp

        # calculation of phi2
        I = np.array([[np.cos(phi1)], [np.sin(phi1)], [0]])
        if abs(phi) < 1e-13 or abs(phi - np.pi) < 1e-13:
            # print "bamm"
            phi2 = 0
            I = C1
        else:
            phi2_temp = np.arccos((I.transpose().dot(C1)[0, 0]) / (np.linalg.norm(I) * np.linalg.norm(C1)))
            # print I.transpose().dot(C1)[0, 0], (np.linalg.norm(I)*np.linalg.norm(C1)), \
            #     (I.transpose().dot(C1)[0, 0])/(np.linalg.norm(I)*np.linalg.norm(C1))
            # print phi2_temp
            if C1[2] >= 1e-13:
                # print "fall 1"
                phi2 = phi2_temp
            elif C1[2] <= -1e-13:
                # print "fall 2"
                phi2 = 2 * np.pi - phi2_temp
            elif abs(C1[2]) < 1e-13:
                # print "fall 3"
                if abs(I[0] - C1[0]) < 1e-13 and abs(I[1] - C1[1]) < 1e-13 and abs(I[2] - C1[2]) < 1e-13:
                    phi2 = 0
                elif abs(I[0] + C1[0]) < 1e-13 and abs(I[1] + C1[1]) < 1e-13 and abs(I[2] + C1[2]) < 1e-13:
                    phi2 = np.pi
        try:
            # print "phi1: ", float(phi1)
            phi1 = float(phi1)
        except ValueError:
            phi1 = phi1[0]

        try:
            phi = float(phi)
        except ValueError:
            phi = phi[0]

        try:
            phi2 = float(phi2)
        except ValueError:
            phi2 = phi2[0]

        # L_1 = np.dot(self.g1(phi_, psi, phi_2).transpose(), np.array([[1], [0], [0]]))
        # L_2 = np.dot(self.g1(phi_, psi, phi_2).transpose(), np.array([[0], [1], [0]]))
        # L_3 = np.dot(self.g1(phi_, psi, phi_2).transpose(), np.array([[0], [0], [1]]))

        # crystalframe = (np.dot(g_, np.array([[1], [0], [0]])),
        #                 np.dot(g_, np.array([[0], [1], [0]])),
        #                 np.dot(g_, np.array([[0], [0], [1]])),
        #                 rad_to_deg(self.__params['phi2_']))
        #
        # g = self.g(phi1, phi, phi2).transpose()
        # crystalframe2 = (np.dot(g, np.array([[1], [0], [0]])),
        #                  np.dot(g, np.array([[0], [1], [0]])),
        #                  np.dot(g, np.array([[0], [0], [1]])),
        #                  rad_to_deg(self.__params['phi2_']))
        #
        # title = "phi: %i, psi: %i, hkl: %i%i%i" % (rad_to_deg(self.__params["phi"]), rad_to_deg(self.__params["psi"]),
        #                                            h, k, l)
        # Q = self.g1(phi_, psi, phi_2).transpose().dot(np.array([[0], [0], [1]]))
        # # cplot.plot_coordinatframe(L_1, L_2, L_3, Q, titel=title, crystalframe=crystalframe,
        # #                           crystalframe2=crystalframe2, I=I, rotation=rad_to_deg(self.__params["phi2_"]))
        #
        # Q_C = g.dot(
        #     np.array([[np.sin(phi_b) * np.cos(betha_b)], [np.sin(phi_b) * np.sin(betha_b)], [np.cos(phi_b)]]))
        # # print "angle between c_3 and q: ", rad_to_deg(float(np.arccos(crystalframe2[2].transpose().dot(Q)[0, 0])))
        # delta_angle = rad_to_deg(
        #     float(np.arccos(Q_C.transpose().dot(Q)[0, 0] / (np.linalg.norm(Q_C) * np.linalg.norm(Q)))))
        # if abs(float(Q_C.transpose().dot(Q)[0, 0] / (np.linalg.norm(Q_C) * np.linalg.norm(Q))) - 1) < 1e-13:
        #     delta_angle = 0.0
        # print "angle between q_c and q: ", delta_angle, Q_C.transpose().dot(Q)[0, 0]/(np.linalg.norm(Q_C)*np.linalg.norm(Q))
        return phi1, phi, phi2  # , delta_angle

    def calc_delta_vals(self, h, k, l, psi, phi, r1, step):
        """
        psi, phi, r1, step in grad
        :param h:
        :param k:
        :param l:
        :param psi:
        :param phi:
        :param r1:
        :param step:
        :return: deulerangles in rad
        """
        self.__params["psi"] = deg_to_rad(psi)
        self.__params["phi"] = deg_to_rad(phi)
        self.__params['phi2_'] = deg_to_rad(r1)  # rotation in crystal frame
        self.__params['phi_b'] = self.calc_phi_b(h, k, l)
        self.__params['betha_b'] = self.calc_betha_b(h, k, l)
        r1 = deg_to_rad(r1)
        r2 = r1 + deg_to_rad(step)
        phi1, phi_, phi2 = self.calc_eulerangles(r1 % (2 * np.pi), h, k, l)
        phi1_1, phi__1, phi2_1 = self.calc_eulerangles(r2 % (2 * np.pi), h, k, l)

        dphi1 = abs(phi1_1 - phi1)
        if dphi1 > deg_to_rad(300):
            dphi1 = abs(dphi1 - 2 * np.pi)

        dphi_ = abs(phi__1 - phi_)
        if dphi_ > deg_to_rad(300):
            dphi_ = abs(dphi1 - 2 * np.pi)

        dphi2 = abs(phi2_1 - phi2)
        if dphi2 > deg_to_rad(300):
            dphi2 = abs(dphi1 - 2 * np.pi)

        dphis = [dphi1, dphi_, dphi2]
        phis = [rad_to_deg(phi1), rad_to_deg(phi_), rad_to_deg(phi2)]
        return dphi1, dphi_, dphi2

    def f_interpolate(self, phi1, phi, phi2):
        """
        this function interpolates between the neighboring values.
        It is a simple linear interpolation.
        just using a one rank taylor approximation
        It should increase the precision of the calculations.
        :param phi1: Eulerangle 1
        :param phi: Eulerangle 2
        :param phi2: Eulerangle 3
        :return: Value of the ODF
        """
        phi1_0 = int(phi1 / self.__stepwidth) * int(self.__stepwidth)
        phi_0 = int(phi / self.__stepwidth) * int(self.__stepwidth)
        phi2_0 = int(phi2 / self.__stepwidth) * int(self.__stepwidth)

        f_0 = self.f(phi1_0, phi_0, phi2_0)  # 0. rank approx
        f_1_phi1 = ((self.f(phi1_0 + int(self.__stepwidth), phi_0, phi2_0) - f_0) / int(self.__stepwidth)) * \
                   (phi1 - phi1_0)

        f_1_phi = ((self.f(phi1_0, phi_0 + int(self.__stepwidth), phi2_0) - f_0) / int(self.__stepwidth)) * \
                  (phi - phi_0)

        f_1_phi2 = ((self.f(phi1_0, phi_0, phi2_0 + int(self.__stepwidth)) - f_0) / int(self.__stepwidth)) * \
                   (phi2 - phi2_0)

        f = f_0 + f_1_phi + f_1_phi1 + f_1_phi2
        return f

    def f_interpol(self, phi1, phi, phi2):
        phi1 = int(phi1 / self.__interpolated_stebwidth)
        phi = int(phi / self.__interpolated_stebwidth)
        phi2 = int(phi2 / self.__interpolated_stebwidth)
        res = self.__ODF_interpolated[phi2][phi][phi1]
        return res

    def interpolate_whole_odf(self):
        if not self.__odf_is_interpolated:
            step = 1
            number_of_phi1vals = int(self.__phi1_max / step + 1)
            number_of_Phivals = int(self.__Phi_max / step + 1)
            number_of_phi2vals = int(self.__phi2_max / step + 1)

            interpolated_odf = np.zeros((number_of_phi2vals, number_of_Phivals, number_of_phi1vals))
            try:
                interpolated_odf = np.load(self.__path_to_data + 'ODF_data_interpolated' + self.name)
            except IOError:
                for phi1 in xrange(0, 360, step):
                    for phi in xrange(0, 360, step):
                        for phi2 in xrange(0, 360, step):
                            f = self.f_interpolate(phi1, phi, phi2)
                            interpolated_odf[phi2, phi, phi1] = f
                np.save(self.__path_to_data + 'ODF_data_interpolated' + self.name, interpolated_odf)
            self.__interpolated_stebwidth = step
            self.__ODF_interpolated = interpolated_odf.tolist()
            self.__odf_is_interpolated = True
        else:
            pass

    def plot_polfigure(self, h, k, l):
        """
        this function should plot the hkl pole figure (once it is implemented :D)
        :param h: Millerindex
        :param k: Millerindex
        :param l: Millerindex
        :return: None
        """
        # self.interpolate_whole_odf()
        step = deg_to_rad(5)
        PHI = np.arange(0, 2 * np.pi + step, step)
        PSI = np.arange(0, np.pi / 2 + step, step)
        VAL = np.zeros((PHI.size, PSI.size))
        endval = len(PHI) * len(PSI)
        i = 0
        t1 = tm.clock()
        print self.calc_al_symetrical_identical_hkl(h, k, l)
        for m in range(len(PHI)):
            for n in range(len(PSI)):
                phi = PHI[m]
                psi = PSI[n]
                val = self.integrate_(phi, psi, h, k, l)
                VAL[m, n] = val
                cli_progress_test(i=i, end_val=endval)
                i += 1
        t2 = tm.clock()

        dt = t2 - t1
        print "time: h: {}, m: {}, s: {}".format(int(dt / 3600), int((dt % 3600) / 60), (dt % 3600) % 60)
        # r, theta = np.meshgrid(PSI, PHI)
        # print r, PSI
        # print theta, PHI

        # ----------
        # the plot
        # ----------
        r, theta = np.meshgrid(PSI, PHI)
        # plt.figure("pole figure {}{}{}\n".format(h, k, l))
        # fig = Figure()
        # axs = fig.add_subplot(111, projection='polar')  # 1, 1, figsize=(12,5),subplot_kw=dict(projection='polar'))
        # pp = plt.axes(polar=True)
        # axs.set_title("pole figure {}{}{}".format(h, k, l))
        # pp.gird(linestyle='-.')
        # pp.pcolormesh(VAL, cmap = cm.gist_rainbow)  # theta, r,
        fig, axs = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
        p1 = axs.contourf(theta, np.degrees(r), VAL, 100)
        # p2 = axs[1].contourf(theta, r, values2, 100)
        cbar = plt.colorbar(p1, ax=axs)
        axs.set_title("pole figure {}{}{}\n".format(h, k, l))
        # pl = axs[0].contourf(theta, np.degrees(r), VAL)
        # CB = plt.colorbar(pl, ax=axs)
        # plt.colorbar(pl, ax=pp)
        # plt.colorbar(pp)
        # pp.scatter(theta, r, s=27, c=VAL, marker='o')  # this is working

        # plt.colorbar(mappable=VAL, cax=pp)
        # fig.set_title("pole figure {}{}{}".format(h, k, l))

        # ax.contourf(theta, r, VAL)
        # fig.colorbar(ax)
        plt.show()

    def integrate(self, A, phi, psi, h, k, l, *args):
        # performe the integration around q//h
        self.__params["phi"] = phi
        self.__params["psi"] = psi
        self.__params['phi_b'] = self.calc_phi_b(h, k, l)
        self.__params['betha_b'] = self.calc_betha_b(h, k, l)
        u = args[0]
        w = args[1]
        step = 5
        res = 0
        counter = 0
        m = self.m(phi, psi)
        g1 = self.g1(phi, psi, phi2=0)
        # print "g1", g1
        for i in range(0, 360, step):
            counter += 1
            r = deg_to_rad(i)
            phi1, phi, phi2 = self.calc_eulerangles(r % (2 * np.pi), h, k, l)
            # dphi1, dphi, dphi2 = self.calc_delta_vals(h, k, l, rad_to_deg(psi), rad_to_deg(phi), r, step)  # [0]
            euler = (rad_to_deg(phi1), rad_to_deg(phi), rad_to_deg(phi2))
            phi1, phi, phi2 = euler
            # print self.f(phi1, phi, phi2), counter, phi1, phi, phi2
            phi2_ = deg_to_rad(i)
            phi_b = self.calc_phi_b(h, k, l)
            beta_b = self.calc_betha_b(h, k, l)
            g2 = self.g2(phi2_=phi2_, phi_b=phi_b, beta_b=beta_b)
            # print "A_ijkl: ", A(g2, g1, *args)
            res += A(g2, g1, *args) * g1[2, u] * g1[2, w] * self.f(phi1, phi, phi2) * deg_to_rad(step)  # \
            # print res
            # * np.sin(deg_to_rad(phi)) * dphi1 * dphi * dphi2 # * (2 * np.pi)  # ** 2
            #  * \

        return res / (2 * np.pi)

    def integrate_(self, phi, psi, h, k, l, *args):
        """
        integrate around scattering vector using the original ODF (not interpolated)
        :param phi:
        :param psi:
        :param h:
        :param k:
        :param l:
        :param args:
        :return:
        """
        # performe the integration around q//h
        self.__params["phi"] = phi
        self.__params["psi"] = psi

        self.__params['phi_b'] = self.calc_phi_b(h, k, l)
        self.__params['betha_b'] = self.calc_betha_b(h, k, l)
        step = 5
        res = 0
        # hkls = self.calc_al_symetrical_identical_hkl(h, k, l)
        # for h, k, l in hkls:
        for i in range(0, 360, step):
            r = deg_to_rad(i)
            phi1, phi, phi2 = self.calc_eulerangles(r % (2 * np.pi), h, k, l)
            euler = (rad_to_deg(phi1), rad_to_deg(phi), rad_to_deg(phi2))
            phi1, phi, phi2 = euler
            # dphi1, dphi, dphi2 = self.calc_delta_vals(h, k, l, rad_to_deg(psi), rad_to_deg(phi), i, step)  # [0]

            res += self.f(phi1, phi, phi2) * deg_to_rad(step)
            # \ * np.sin(deg_to_rad(phi)) * dphi1 * dphi * dphi2  # / (2 * np.pi)  # ** 2

        return res / (2 * np.pi)

    def calc_al_symetrical_identical_hkl(self, h, k, l):
        """
        calc all symmetrical equivalent hkl plains
        :param h:
        :param k:
        :param l:
        :return:
        """
        res = []
        if self.crystal_symmetry == "m-3m":
            for i in IT.permutations([h, k, l]):
                res.append(list(i))
                res.append((-np.array(i)).tolist())
            for i in IT.permutations([-h, k, l]):
                res.append(list(i))
                res.append((-np.array(i)).tolist())
            for i in IT.permutations([h, -k, l]):
                res.append(list(i))
                res.append((-np.array(i)).tolist())
            for i in IT.permutations([h, k, -l]):
                res.append(list(i))
                res.append((-np.array(i)).tolist())

        # remove redundand etries
        result = []
        for j in res:
            if j not in result:
                result.append(j)
        return result

    def integrate_interpol(self, phi, psi, h, k, l, *args):
        """
        integrate around scattering vector using the interpolated ODF
        :param phi:
        :param psi:
        :param h:
        :param k:
        :param l:
        :param args:
        :return:
        """
        # performe the integration around q//h
        self.__params["phi"] = phi
        self.__params["psi"] = psi
        step = 5
        res = 0
        for i in range(0, 360, step):
            r = deg_to_rad(i)
            phi1, phi, phi2 = self.calc_eulerangles(r % (2 * np.pi), h, k, l)
            euler = (rad_to_deg(phi1), rad_to_deg(phi), rad_to_deg(phi2))
            phi1, phi, phi2 = euler
            # dphi1, dphi, dphi2 = self.calc_delta_vals(h, k, l, rad_to_deg(psi), rad_to_deg(phi), i, step)  # [0]

            res += self.f(phi1, phi, phi2) * deg_to_rad(step)
            # \ * np.sin(deg_to_rad(phi)) * dphi1 * dphi * dphi2  # / (2 * np.pi)  # ** 2

        return res / (2 * np.pi)

    @property
    def integral_over_all_orientations_g(self):
        sum_total = 0

        t1 = tm.clock()
        for k in xrange(0, self.__phi2_max + self.__stepwidth, self.__stepwidth):  # sum over all phi2 vals
            for j in xrange(0, self.__Phi_max + self.__stepwidth, self.__stepwidth):  # sum over all Phi vals
                sum_1 = 0
                sum_2 = 0
                for i in xrange(0, self.__phi1_max + self.__stepwidth,
                                self.__stepwidth):  # sum over all phi1 vals (with Phi_j)
                    sum_1 += self.f((i + self.__stepwidth) % 360, (j + self.__stepwidth) % 180, k % 360) \
                             + self.f(i % 360, (j + self.__stepwidth) % 180, k % 360) \
                             + self.f((i + self.__stepwidth) % 360, (j + self.__stepwidth) % 180,
                                      (k + self.__stepwidth) % 360) \
                             + self.f(i % 360, (j + self.__stepwidth) % 180, (k + self.__stepwidth) % 360)

                    # for i in xrange(0, self.__phi1_max, self.__stepwidth):  # sum over all phi1 vals (with Phi_j+1)
                    sum_2 += self.f((i + self.__stepwidth) % 360, j % 180, k % 360) \
                             + self.f(i % 360, j % 360, k % 360) \
                             + self.f((i + self.__stepwidth) % 360, j % 180, (k + self.__stepwidth) % 360) \
                             + self.f(i % 360, j % 180, (k + self.__stepwidth) % 360)

                sum_total += np.sin(deg_to_rad((j + self.__stepwidth) % 180)) * sum_1 \
                             + np.sin(deg_to_rad(j % 180)) * sum_2
        sum_total = sum_total * deg_to_rad(self.__stepwidth) ** 3 / (8 * 8 * np.pi ** 2)
        sum_total = sum_total  # * \
        # 360 / self.__phi1_max * \
        # 180 / self.__Phi_max * \
        # 360 / self.__phi2_max
        t2 = tm.clock()
        dt = t2 - t1
        print "Integral: ", sum_total, "Time to calc it: ", dt
        return sum_total

    def integrate_a_over_all_orientations_g(self, inner_sum, *args):
        sum_total = 0
        # t1 = tm.clock()
        for phi2 in xrange(0, self.__phi2_max + self.__stepwidth, self.__stepwidth):  # sum over all phi2 vals
            for phi in xrange(0, self.__Phi_max + self.__stepwidth, self.__stepwidth):  # sum over all Phi vals
                sum_1 = 0
                for phi1 in xrange(0, self.__phi1_max + self.__stepwidth, self.__stepwidth):  # sum over all phi1 vals
                    sum_1 += inner_sum(phi1, phi, phi2, *args) * self.f(phi1, phi, phi2)
                sum_total += np.sin(deg_to_rad(phi)) * sum_1
        sum_total = sum_total * deg_to_rad(self.__stepwidth) ** 3 / (8 * np.pi ** 2)
        # t2 = tm.clock()
        # dt = t2 - t1

        # print "Integral: ", sum_total, "Time to calc it: ", dt
        return sum_total

    def __transform_to_on_3_dim_np_array(self):  # similar to the GKSS file
        """
        sets the odf data in an other representation:
        array[i,j,k]
        i = phi2 val
        j = phi  val
        k = phi1 val
        """
        number_of_phi1vals = int(self.__phi1_max / self.__stepwidth + 1)
        number_of_Phivals = int(self.__Phi_max / self.__stepwidth + 1)
        number_of_phi2vals = int(self.__phi2_max / self.__stepwidth + 1)

        data_phi2 = np.zeros((number_of_phi2vals, number_of_Phivals, number_of_phi1vals))
        # i = 0
        # k = 0
        # l = 0
        # m = 0
        t1 = tm.clock()
        self.__ODF_data = data_phi2

        print "Name: ", self.name
        try:
            self.__ODF_data = np.load(self.__path_to_data + 'ODF_data_' + self.name + '.npy')
            print "use saved ODF"
        except IOError:
            print "no saved ODF found"
            for i in self.__data:
                phi2 = int(i[2] / self.__stepwidth)
                phi = int(i[1] / self.__stepwidth)
                phi1 = int(i[0] / self.__stepwidth)
                self.__ODF_data[phi2, phi, phi1] = i[-1]
            # for i in xrange(0, self.__phi2_max + self.__stepwidth, self.__stepwidth):
            #     print "phi2: ", i
            #     for j in xrange(0, self.__Phi_max + self.__stepwidth, self.__stepwidth):
            #         print "phi: ", j
            #         for k in xrange(0, self.__phi1_max + self.__stepwidth, self.__stepwidth):
            #             for m in self.__data:
            #                 if m[2] == i and m[1] == j and m[0] == k:
            #                     data_phi2[i / self.__stepwidth, j / self.__stepwidth, k / self.__stepwidth] = m[3]
            #                     break
            # self.__ODF_data = data_phi2
            print self.__ODF_data
            np.save(self.__path_to_data + 'ODF_data_' + self.name, self.__ODF_data)
        self.__ODF_data = self.__ODF_data.tolist()
        t2 = tm.clock()
        dt = t2 - t1
        print "Time necessary for conversion: ", dt, "sec."
        ''' ...old code needs a special m-tex format (the code abov can deal each m-tex format)...
        while i <len(self.__data):  #loop over the howl data
            while m<number_of_phi2vals:  #loop over all Phi2 values
                while self.__data[i,2]==m*5:  #loop over all Phi2 (from 0 to 90°)
                    for l in xrange(number_of_phi1vals):
                        data_phi2[m,k,l] = self.__data[i,3]
                        i+=1
                    if k==number_of_Phivals-1 :
                        k=0
                        l=0
                        break
                    else:
                        k+=1
                m+=1
            break
        '''

    def plot_odf(self):
        number_of_phi2vals = int(90 / 5 + 1)
        origin = 'lower'
        phi1 = np.arange(0, self.__phi1_max + 5, 5)
        Phi = np.arange(0, self.__Phi_max + 5, 5)
        X, Y = np.meshgrid(phi1, Phi)
        # Now make a contour plot with the levels specified,
        # and with the colormap generated automatically from a list
        # of colors.
        levels = np.arange(self.__imin, self.__imax, (self.__imax - self.__imin) / 8.)
        cmap = plt.cm.get_cmap()
        cmap.set_under("magenta")
        cmap.set_over("yellow")
        # Note: contouring simply excludes masked or nan regions, so
        # instead of using the "bad" colormap value for them, it draws
        # nothing at all in them.  Therefore the following would have
        # no effect:
        # cmap.set_bad("red")

        fig, axs = plt.subplots(number_of_phi2vals / 2 + 1, 2)
        # print zip(axs.ravel())

        for i in xrange(number_of_phi2vals / 2 + 1):
            cs = axs[i, 0].contourf(X, Y, self.__ODF_data[2 * i, :, :], levels, cmap=cmap, extend=str(2 * i * 5),
                                    origin=origin)
            try:
                ax2 = axs[i, 1].contourf(X, Y, self.__ODF_data[2 * i + 1, :, :], levels, cmap=cmap,
                                         extend=str(2 * i * 5 + 5), origin=origin)
            except (IndexError):
                pass
            axs[i, 0].set_title(str(2 * 5 * i))
            axs[i, 1].set_title(str(2 * 5 * i + 5))

        for i in xrange(len(axs) - 1):
            plt.setp([a.get_xticklabels() for a in axs[i, :]], visible=False)
        plt.setp([a.get_yticklabels() for a in axs[:, 1]], visible=False)
        plt.show()

    @staticmethod
    def add_rot():
        """
        rotates the ND perp to axis of the stick and the LD parallel to the stick axis
        :return:
        """
        psi, phi, phi2 = deg_to_rad(0), -deg_to_rad(90), deg_to_rad(0)
        # phi += np.pi / 2
        # O = np.array([[np.cos(phi), np.sin(phi), 0.],
        #               [-np.sin(phi), np.cos(phi), 0.],
        #               [0., 0., 1.]
        #               ]
        #              )
        res = np.array([[np.cos(phi), 0., -np.sin(phi)],
                        [0., 1., 0.],
                        [np.sin(phi), 0., np.cos(phi)]
                        ]
                       )
        return res


'''
Printing the result and call the fittingfunctions
'''


def print_result_nicely(res):
    pars = lm.fit_report(res.params)
    out = "Message:        %s\
    \nCovar matrix:   \n%s\
    \nSuccess:        %s\
    \nChisqr:         %f\
    \nReducec Chisqr: %f\
    \ndeg of freedom: %i\
    \n# of datapoint: %i\
    \nnfev:           %i\
    \nParameters:\n%s" \
          % (res.message, str(res.covar), res.success, res.chisqr, res.redchi, res.nfree, res.ndata, res.nfev, pars)
    return out


def call_all(data, Material, filename='result', mode='a'):
    '''
    this function gets the data, which is an object of type InParam, the Material name, the name of the file where the result should be stored
    and the writing mode. If mode = 'a' the result is just appandet to the file, if mode = 'w', the outputfile will be overwriten.
    '''
    Gamma = []
    hkl = data.data['hkl']

    for i in hkl:
        Gamma.append(Gama(float(i[0]), float(i[1]), float(i[2])))

    s1 = data.data['s1']
    s2 = data.data['s2']
    s1_err = data.data['s1_err']
    s2_err = data.data['s2_err']
    print "call all: s1, s2, s1_err, s2_err:"
    print s1
    print s2
    print s1_err
    print s2_err

    result = open('%s.txt' % (filename), mode)  # file to write the results
    # print "\n----------------------------------------------------\n"
    # print "Material: %s \n" %(Material)
    # print "de Wit:\n"
    # print "----------------------------------------------------\n"
    if mode == 'a':
        result.write(
            "\n\n\n\n\n\
        \nMaterial: %s \n----------------------------------------------------\
        \nde Wit:\
        \n----------------------------------------------------\n" % (Material))
    else:
        result.write(
            "Material: %s \n----------------------------------------------------\
        \nde Wit:\
        \n----------------------------------------------------\n" % (Material))
    deWit = fitting_deWit(Gamma, s1, s2, s1_err, s2_err)
    deWit_ = print_result_nicely(deWit)
    result.write(str(deWit_))
    # print 'dewit: ', deWit


    # print "----------------------------------------------------\n"
    # print "Material: %s\n" %(Material)
    # print "BHM:\n"
    # print "----------------------------------------------------\n"

    result.write(
        "\n\n----------------------------------------------------\
    \nBHM:\
    \n----------------------------------------------------\n")
    BHM = fitting_BHM(Gamma, s1, s2, s1_err, s2_err)
    BHM_ = print_result_nicely(BHM)
    result.write(str(BHM_))
    # print BHM

    # print "\n\n----------------------------------------------------\n"
    # print "Material: %s\n" %(Material)
    # print "Hill:\n"
    # print "----------------------------------------------------\n"
    result.write(
        "\n\n----------------------------------------------------\
    \nHill:\
    \n----------------------------------------------------\n")
    Hill = fitting_Hill(Gamma, s1, s2, s1_err, s2_err)
    Hill_ = print_result_nicely(Hill)
    result.write(str(Hill_))


'''
------------------------------------------------------------------------------
Some helpefull Funktions
------------------------------------------------------------------------------
'''


def deg_to_rad(deg):
    return np.pi * deg / 180.


def rad_to_deg(rad):
    return 180. * rad / np.pi


'''
Below there is an Implementation for a "Brute force" lesast square method
and a code snipet to plot the REK's -s_1 over 3Gamma to check for consitency
'''

'''
# Create toy data for curve_fit.
Gamma=[]
hkl     = [[0.,0.,2.],[0.,2.,2.],[2.,2.,2.],[3.,1.,1.]]
for i in hkl:
    Gamma.append(Fittingtools.Gamma(i[0],i[1],i[2]))
    print "hkl %d%d%d -> Gamma: " %(i[0],i[1],i[2]), Fittingtools.Gamma(i[0],i[1],i[2])
xdata = np.array(Gamma)
y=np.array([[-1.758,7.522],[-1.121,5.499],[-1.010,5.342],[-1.692,6.852]])
ydata = y*np.power(10.,-12)
#wights
sigma=np.array([[.084,.150],[.062,.101],[.103,.172],[.099,.166]])*np.power(10.,-12)

for i in range(len(Gamma)):
    print "hkl: %d%d%d "%(hkl[i][0],hkl[i][1],hkl[i][2]),"s_1:\t   s_2:",\
    "\nReus:     ", Reus(Gamma[i],
                         c_11=114*np.power(10.,9),
                         c_12=65.3*np.power(10.,9),
                         c_44=28.5*np.power(10.,9)),\
          "\nRV:     ", RV(Gamma[i],
                         c_11=114*np.power(10.,9),
                         c_12=65.3*np.power(10.,9),
                         c_44=28.5*np.power(10.,9)),\
          "\nVoigt:    ", Voigt__(Gamma[i],
                         c_11=114*np.power(10.,9),
                         c_12=65.3*np.power(10.,9),
                         c_44=28.5*np.power(10.,9)),\
          "\nBHM:      ", BHM(Gamma[i],
                              c_11=114*np.power(10.,9),
                         c_12=65.3*np.power(10.,9),
                         c_44=28.5*np.power(10.,9)),\
          "\nBHM_dW:   ", BHM_dW(Gamma[i],
                         c_11=114*np.power(10.,9),
                         c_12=65.3*np.power(10.,9),
                         c_44=28.5*np.power(10.,9))
'''
'''
print Gamma
dW=[]
BHM_ = []
g=[]
for i in np.linspace(0,1/3.,50):
    dW.append(-BHM_dW(i)[0])
    BHM_.append(-BHM(i)[0])
    g.append(i)
print BHM_dW(np.array(Gamma))
plt.figure('rv')
plt.plot(np.array(g)*3,-Reus(np.array(g),
                         c_11=224.9*np.power(10.,9),
                         c_12=122.2*np.power(10.,9),
                         c_44=120.7*np.power(10.,9))[0], 'r-', label='Reus')
plt.plot(np.array(g)*3,-Voigt__(np.array(g),
                         c_11=224.9*np.power(10.,9),
                         c_12=122.2*np.power(10.,9),
                         c_44=120.7*np.power(10.,9))[:,0], 'b--', label='Voigt')
g = np.array(g)*3
plt.plot(g, dW, 'g', label='de Wit')
plt.plot(g, BHM_, 'k-.', label='BHM')
plt.xlabel(r'$3\Gamma$')
plt.ylabel(r'$s_1$')
plt.legend(loc='upper right', numpoints = 1)
plt.show()
'''

# Compute chi-square manifold.
'''
#Values for calculation with Steps = 501
Steps = 501  # grid size
Chi2Manifold = numpy.zeros([Steps,Steps,Steps])  # ,Stepsallocate grid
c_11min = (205.*np.power(10.,9))  # minimal value of a covered by grid
c_11max = (245.*np.power(10.,9))  # maximal value of a covered by grid
c_12min = (105.*np.power(10.,9))  # minimal value of b covered by grid
c_12max = (140.*np.power(10.,9))  # maximal value of b covered by grid
c_44min = (105.*np.power(10.,9))  # minimal value of c covered by grid
c_44max = (140.*np.power(10.,9))  # maximal value of c covered by grid

Chi2Manifold = np.load('Chi2_Values_Modell_dW.npy')
'''

'''
Steps = 51  # grid size
Chi2Manifold = numpy.zeros([Steps,Steps,Steps])  # ,Stepsallocate grid
c_11min = (220.*np.power(10.,9))  # minimal value of a covered by grid
c_11max = (230.*np.power(10.,9))  # maximal value of a covered by grid
c_12min = (115.*np.power(10.,9))  # minimal value of b covered by grid
c_12max = (130.*np.power(10.,9))  # maximal value of b covered by grid
c_44min = (115.*np.power(10.,9))  # minimal value of c covered by grid
c_44max = (130.*np.power(10.,9))  # maximal value of c covered by grid



for s1 in xrange(Steps):
    print "s1: ",s1
    for s2 in xrange(Steps):
        for s3 in xrange(Steps):
            # Current values of (c_11,c_12,c_44) at grid position (s1,s2,s3).
            c_11 = c_11min + (c_11max - c_11min)*float(s1)/(Steps-1)
            c_12 = c_12min + (c_12max - c_12min)*float(s2)/(Steps-1)
            c_44 = c_44min + (c_44max - c_44min)*float(s3)/(Steps-1)
            # Evaluate chi-squared.
            chi2 = 0.0
            for n in range(len(xdata)):
                for i in range(len(ydata[n])):
                    residual_ = (ydata[n][i] - BHM(xdata[n], c_11, c_12, c_44)[i])/sigma[n][i]
                    chi2 = chi2 + residual_*residual_
            Chi2Manifold[s3,Steps-1-s2,s1] = chi2  # write result to grid.Steps-1-



np.save('Chi2_Values_Modell_Voigt',Chi2Manifold)

#Chi2Manifold = np.load('Chi2_Values_Modell_BHM.npy')
h = np.where(Chi2Manifold == Chi2Manifold.min())
print h
#amtot = np.argmax(h)
#
chi = Chi2Manifold[h[0][0]][h[1][0]][h[2][0]]
print   "c_11:  ", c_11min + (c_11max - c_11min)*float(h[2][0])/(Steps-1),\
      "\nc_12:  ", c_12min + (c_12max - c_12min)*float(Steps-1-h[1][0])/(Steps-1), \
      "\nc_44:  ", c_44min + (c_44max - c_44min)*(float(h[0][0]))/(Steps-1),\
      "\nchi:   ", Chi2Manifold[h[0][0]][h[1][0]][h[2][0]]

# Plot grid.
#, figsize=(8,4.5)
#plt.subplots_adjust(left=0.09, bottom=0.09, top=0.97, right=0.99)
# Plot chi-square manifold.
#for i in range(Steps):
plt.figure('Reus_600')
image = plt.imshow(Chi2Manifold[h[0][0]], vmax = chi + chi*0.01,#
              extent=[c_11min, c_11max, c_12min, c_12max])
              #extent=[c_11min, c_11max, c_12min, c_12max])
# Plot where curve-fit is going to for a couple of initial guesses.
#for a_initial in -6.0, -4.0, -2.0, 0.0, 2.0, 4.0:
    # Initial guess.
#    x0   = numpy.array([a_initial, -3.5])
#    xFit = optimization.curve_fit(func, xdata, ydata, x0, sigma)[0]
#    plt.plot([x0[0], xFit[0]], [x0[1], xFit[1]], 'o-', ms=4,
#                 markeredgewidth=0, lw=2, color='orange')
plt.colorbar(image)  # make colorbar
plt.xlim(c_11min, c_11max)
plt.ylim(c_12min, c_12max)
#plt.xlim(c_11min, c_11max)
#plt.ylim(c_12min, c_12max)
plt.xlabel(r'$c_{11}$', fontsize=24)
plt.ylabel(r'$c_{12}$', fontsize=24)
#plt.savefig('chi^2, c_11 und c_12, Modell_BHM.pdf')
plt.show()
'''
