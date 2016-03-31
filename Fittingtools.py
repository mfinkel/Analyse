# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 18:28:42 2015

@author: mfinkel
"""
import numpy as np
from scipy.optimize import curve_fit
import scipy.odr.odrpack as odr
import matplotlib.pyplot as plt
from lmfit.models import PseudoVoigtModel


def breite(x, y):
    help1 = 0.
    help2 = 0.
    rechts_l = 0.
    rechts_h = 0.
    links_l = 0.
    links_h = 0.
    lambdamax = 0.
    mitte = 0.
    maximum = max(y) / 2.
    for i in xrange(-1, len(y) + 1):
        if i == len(y) - 2:
            break
        if i >= 1 and i <= len(y):
            help2 = y[i]
            help1 = y[i + 1]
        if (help1 < maximum and help2 > maximum):
            rechts_h = x[i + 1]
            rechts_l = x[i]
        elif (help2 < maximum and help1 > maximum):
            links_l = x[i]
            links_h = x[i + 1]
        elif max(y) == y[i]:
            lambdamax = float(x[i])
    mitte = ((rechts_l - links_h) + (rechts_h - links_l)) / 2.
    return [float(mitte), lambdamax]


def lin(x, a, b):
    return a * x + b


def lin_fit(x, y):
    x = np.array(x)
    y = np.array(y)
    p_final, pcov = curve_fit(lin, x, y, maxfev=100000000)
    return p_final


def gaussian(x, amplitude, xcenter, sigma, offset):
    g = amplitude * np.exp(-0.5 * ((x - xcenter) / (sigma)) ** 2) + offset
    return g


def gaussian_lin(x, amplitude, xcenter, sigma, offset, a):
    g = amplitude * np.exp(-0.5 * ((x - xcenter) / (sigma)) ** 2) + offset + a * np.array(x)
    return g


def gaus_psoido_vioight(x, amplitude, xcenter, sigma, offset, a, g_):
    g = amplitude * np.exp(-0.5 * ((x - xcenter) / sigma) ** 2 - g_ * np.abs(x - xcenter)) + offset + a * np.array(x)
    return g


def gauss_fit(x, y):
    linear = odr.Model(gaussian)
    mydata = odr.RealData(x, y)
    myodr = odr.ODR(mydata, linear, beta0=guesspara(x, y))
    myoutput = myodr.run()
    myoutput.pprint()


def gauss_fitting(x_list, av_count_l):
    p_guess = guesspara(x_list, av_count_l)
    av_count_l_y = []
    av_count_l_y = av_count_l
    p_final, pcov = curve_fit(gaussian, x_list, av_count_l_y, p_guess, maxfev=100000000)
    # y = gaussian(x_list, p_final[0], p_final[1], p_final[2], p_final[3])
    # plt.figure("fit")
    # Two subplots, the axes array is 1-d
    '''
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(x_list, av_count_l,'b.',label='Data')
    axarr[0].plot(x_list,y, 'r-',label='fit')
    axarr[0].set_title('Sharing X axis')
    axarr[1].scatter(x_list, av_count_l-y)
    plt.show()
    '''
    return [p_final[1], p_final[2]]


def gauss_lin_fitting(x_list, av_count_l, plot=False):
    '''
    fitting the psoido voigt funktion against the data
    '''
    p_guess = guesspara(x_list, av_count_l)
    p_guess.append(0)
    p_guess.append(0)
    av_count_l_y = []
    av_count_l_y = av_count_l  # gaussian_lin
    p_final, covar = curve_fit(gaus_psoido_vioight, x_list, av_count_l_y, p_guess, maxfev=100000000)

    y_ = gaus_psoido_vioight(x_list, p_final[0], p_final[1], p_final[2], p_final[3], p_final[4], p_final[5])
    # print covar
    # Two subplots, the axes array is 1-d
    if max(av_count_l) < 150. or max(av_count_l - y_) / max(av_count_l) > 0.2:
        res = ['nan', 'nan']
    else:
        res = [p_final[1], np.sqrt(covar[1, 1])] # p_final[2]
    if plot:
        x = np.linspace(x_list[0], x_list[-1], 300)
        y = gaus_psoido_vioight(x, p_final[0], p_final[1], p_final[2], p_final[3], p_final[4], p_final[5])
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(x_list, av_count_l, 'b.', label='Theta: %s' % (str(res[0])))
        axarr[0].plot(x, y, 'r-', label='fit')
        axarr[0].set_title('Sharing X axis')
        axarr[1].scatter(x_list, av_count_l - y_)
        axarr[0].legend(loc='upper right', numpoints=1)
        plt.show()

    # print 'Theta: ',p_final[1], ' ', res[0]
    return res

def gauss_lin_fitting_2(x_list, av_count_l, plot=False):
    mod = PseudoVoigtModel()
    p_guess = guesspara(x_list, av_count_l)
    pars = mod.make_params(amplitude=p_guess[0], center=p_guess[1], sigma=p_guess[2])
    weights = []
    for i in av_count_l:
        weights.append(1/np.sqrt(i))
    out = mod.fit(av_count_l, pars, x=x_list, weights=weights)  #
    y_ = out.best_fit
    res = 0
    n = len(x_list)
    if max(av_count_l) < 150. or max(av_count_l - y_) / max(av_count_l) > 0.2:
        res = ['nan', 'nan']
    else:
        # print out.params["center"].value, out.params["center"].stderr
        res = [out.params["center"].value, out.params["center"].stderr] # p_final[2]
    # print "hallo", res
    if plot:
        plt.plot(x_list, av_count_l, "bo")
        plt.plot(x_list, out.init_fit, "k--")
        plt.plot(x_list, out.best_fit, 'r-', label = "T: %.3f\nerr: %.5f" % (res[0], res[1]))
        plt.legend(loc='upper right', numpoints=1)
        plt.show()
    # print (out.fit_report(min_correl=0.75))
    return res


def gauss_fitting_neu(x_list, av_count_l):
    p_guess = guesspara(x_list, av_count_l)
    av_count_l_y = []
    av_count_l_y = av_count_l
    p_final, pcov = curve_fit(gaussian, x_list, av_count_l_y, p_guess, maxfev=10000)
    return [p_final[0], p_final[1], p_final[2], p_final[3]]


def guesspara(x_list, y_list):
    amplitude = max(y_list) - min(y_list)
    x_center = x_list[y_list.index(np.amax(y_list))]
    offset = min(y_list)
    temp_sigma = breite(x_list, y_list)
    sigma = temp_sigma[0] / 2.3548
    p_guess = [amplitude, x_center, sigma, offset]
    return p_guess


# least_square Fitt to determine the compliences

def Gamma(h, k, l):
    g = ((h ** 2) * (k ** 2) + (h ** 2) * (l ** 2) + (k ** 2) * (l ** 2)) / (((h ** 2) + (k ** 2) + (l ** 2)) ** 2)
    return g
