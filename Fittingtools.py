# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 18:28:42 2015

@author: mfinkel
"""
import numpy as np
from scipy.optimize import curve_fit
import scipy.odr.odrpack as odr
import matplotlib.pyplot as plt
from lmfit.models import PseudoVoigtConstModel, PseudoVoigtDoublePeakModel, PseudoVoigtModel, LinearModel, GaussianModel


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

def pseudo_voigt_single_peak_fit(x_list, y_list, weights=None, plot=False, dataset=False, force=False, Chi=False, phase=False):
    """
    fitting one peak with a psoido voigt function
    :param x_list:
    :param y_list:
    :param weights:
    :param plot:
    :return:
    """
    mod = PseudoVoigtConstModel()
    p_guess = guesspara(x_list, y_list)
    pars = mod.make_params(amplitude=p_guess[0], center=p_guess[1], sigma=p_guess[2], const=p_guess[3])
    pars['const'].set(min=0, max=max(y_list))
    pars['center'].set(min=min(x_list), max=max(x_list))
    out = mod.fit(y_list, pars, x=x_list, weights=weights)  #
    y_ = out.best_fit
    chisqr = out.redchi
    chis = out.redchi/y_.max()**2
    chiss = out.chisqr

    compare = max(y_list)-background(x_list, y_list, out.params['center'], out.params['fwhm'])
    if max(y_list) < 150. or max(y_list - y_) / max(y_list) > 0.2 or compare < 500 or chis > 0.8:
        res = ['nan', 'nan']
    else:
        # print out.params["center"].value, out.params["center"].stderr
        res = [out.params["center"].value, out.params["center"].stderr] # p_final[2]
    # print "hallo", res
    if plot:
        plt.figure("Dataset: {}, Phase: {}, force: {}, Chi: {}".format(dataset, phase, force, Chi))
        plt.plot(x_list, y_list, "bo")
        plt.plot(x_list, out.init_fit, "k--")
        try:
            plt.plot(x_list, out.best_fit, 'r-', label = "T: %.3f\nerr: %.5f, chisqr: %.2f\nchiout: %.1f, chiss: %.1f" % (float(res[0]), float(res[1]), chisqr, chis, chiss))
        except TypeError:
            plt.plot(x_list, out.best_fit, 'r-', label = "T: %.3f\nerr: %.5f, chisqr: %.2f\nchiout: %.1f, chiss: %.1f" % (float(res[0]), float(res[1]), chisqr, -1, chiss))
        plt.legend(loc='upper right', numpoints=1)
        plt.show()
    # print (out.fit_report(min_correl=0.75))
    return res

def pseudo_voigt_double_peak_fit(x_list, y_list, weights=None, plot=False, dataset=False, force=False, Chi=False, phase=False):
    """
    fitting one peak with a psoido voigt function
    :param x_list:
    :param y_list:
    :param weights:
    :param plot:
    :return:
    """
    pvmod1 = GaussianModel(prefix='pv1_')
    pvmod2 = GaussianModel(prefix='pv2_')
    linmod = LinearModel(prefix='lm_')
    mod = PseudoVoigtDoublePeakModel()
    # mod = pvmod1 + pvmod2
    p_guess = guesspara_double_peak(x_list, y_list)
    # params = mod.make_params(pv1_amplitud=p_guess[0], pv1_center=p_guess[1], pv1_sigma=p_guess[2])
    # pars = mod.make_params(pv1_amplitud=p_guess[0], pv1_center=p_guess[1], pv1_sigma=p_guess[2],
    #                           pv2_amplitud=p_guess[3], pv2_center=p_guess[4], pv2_sigma=p_guess[5])  # ,
                              # lm_intercept=p_guess[6], lm_slope=0)
    pars = mod.make_params(amplitude1=p_guess[0], center1=p_guess[1], sigma1=p_guess[2],
                           amplitude2=p_guess[3], center2=p_guess[4], sigma2=p_guess[5],
                           const=p_guess[6], a=0)
    # pars['lm_slope'].vary=False
    min1 = p_guess[1] - p_guess[2]
    max1 = p_guess[1] + p_guess[2]
    if abs(min1 - max1) < 0.0001:
        min1 = p_guess[1] - p_guess[1] * 0.1
        max1 = p_guess[1] + p_guess[1] * 0.1
    min2 = p_guess[4] - p_guess[5]
    max2 = p_guess[4] + p_guess[5]
    if abs(min2 - max2) < 0.0001:
        min2 = p_guess[4] - p_guess[4] * 0.1
        max2 = p_guess[4] + p_guess[4] * 0.1

    pars['center1'].set(value=p_guess[1], min=min1, max=max1)
    pars['center2'].set(value=p_guess[4], min=min2, max=max2)
    # pars['pv2_center'].set(value=p_guess[1], min=p_guess[1]-p_guess[2], max=p_guess[1]+p_guess[2])
    # pars['pv2_center'].set(value=p_guess[4], min=p_guess[4]-p_guess[5], max=p_guess[4]+p_guess[5])
    pars['a'].vary=False
    # parss = pvmod2.make_params(pv2_amplitud=p_guess[3], pv2_center=p_guess[4], pv2_sigma=p_guess[5])
    # pars.add_many(parss)
    # pars.update(parss)
    # pars.update(linmod.make_params())
    # pars.update(pvmod2.make_params(pv2_amplitud=p_guess[3], pv2_center=p_guess[4], pv2_sigma=p_guess[5]))
    # pars.update(linmod.make_params(lm_intercept=p_guess[6], lm_slope=0))


    # pars = mod.make_params(amplitude1=p_guess[0], center1=p_guess[1], sigma1=p_guess[2],
    #                        amplitude2=p_guess[3], center2=p_guess[4], sigma2=p_guess[5], const=p_guess[6])

    out = mod.fit(y_list, pars, x=x_list, weights=weights)  #
    y_ = out.best_fit
    chisqr = out.redchi

    if max(y_list) < 150. or max(y_list - y_) / max(y_list) > 0.2:
        res = ['nan', 'nan']
    else:
        # print out.params["center"].value, out.params["center"].stderr
        # res = [out.params["pv1_center"].value, out.params["pv1_center"].stderr,
        #        out.params["pv2_center"].value, out.params["pv2_center"].stderr] # p_final[2]
        res = [out.params["center1"].value, out.params["center1"].stderr,
               out.params["center2"].value, out.params["center2"].stderr] # p_final[2]
    # print "hallo", res
    if plot:
        plt.figure("Dataset: {}, Phase: {}, force: {}, Chi: {}".format(dataset, phase, force, Chi))
        plt.plot(x_list, y_list, "bo")
        plt.plot(x_list, out.init_fit, "k--")
        plt.plot(x_list, out.best_fit, 'r-', label = "T: %.3f\nerr: %.5f, chisqr: %.2f" % (float(res[0]), float(res[1]), chisqr))
        plt.legend(loc='upper right', numpoints=1)
        plt.show()
    # print (out.fit_report(min_correl=0.75))
    return res

def background(x_list, y_list, center, fwhm):
    back = 0
    count = 0
    for i in xrange(len(x_list)):
        if not center-fwhm < x_list[i] < center + fwhm:
            back += y_list[i]
            count += 1
    return back/count

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

def guesspara_double_peak(x_list, y_list):
    # calc amlitudes:
    amplitudes = []
    def df_dx(x1, f1, x2, f2):
        return (f2-f1)/(x2-x1)
    for i in xrange(8, len(y_list)-8):
        dfdx_left=[]
        dfdx_right=[]
        for ii in xrange(1, 4):
            x1 = x_list[i-ii]
            x2 = x_list[i]
            y1 = y_list[i-ii]
            y2 = y_list[i]
            dfdx_left.append(df_dx(x1, y1, x2, y2))
            x1 = x_list[i]
            x2 = x_list[i+ii]
            y1 = y_list[i]
            y2 = y_list[i+ii]
            dfdx_right.append(df_dx(x1, y1, x2, y2))

        dfdx_left = np.sum(np.array(dfdx_left))
        dfdx_right = np.sum(np.array(dfdx_right))
        condition1 = max(y_list[i-5:i+5])<=y_list[i]
        condition2 = (y_list[i] > y_list[i + 1] > y_list[i + 2])
        condition3 = (y_list[i] > y_list[i - 1] > y_list[i - 2])
        condition4 = dfdx_left > 0 and dfdx_right < 0
        if condition1 and y_list[i]-200>min(y_list) and condition4:
            amplitudes.append([i, y_list[i]])

    first, second = amplitudes[0][1], amplitudes[-1][1]
    first_i, second_i = amplitudes[0][0], amplitudes[-1][0]
    amplitude1 = first - min(y_list)
    amplitude2 = second - min(y_list)
    x_center1 = x_list[first_i]
    x_center2 = x_list[second_i]
    index_betwen_peaks = int((second_i-first_i)/2 + first_i)
    offset = min(y_list)
    temp_sigma1 = breite(x_list[0:index_betwen_peaks], y_list[0:index_betwen_peaks])
    temp_sigma2 = breite(x_list[index_betwen_peaks:], y_list[index_betwen_peaks:])
    sigma1 = temp_sigma1[0] / 2.3548
    sigma2 = temp_sigma2[0] / 2.3548
    p_guess = [amplitude1, x_center1, sigma1, amplitude2, x_center2, sigma2, offset]
    return p_guess
# least_square Fitt to determine the compliences

def Gamma(h, k, l):
    g = ((h ** 2) * (k ** 2) + (h ** 2) * (l ** 2) + (k ** 2) * (l ** 2)) / (((h ** 2) + (k ** 2) + (l ** 2)) ** 2)
    return g
