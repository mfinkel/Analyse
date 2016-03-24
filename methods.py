# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 12:59:47 2015

@author: mfinkel
"""
import re
import numpy as np
import matplotlib.pyplot as plt
import math
import DataCursor
import Fittingtools
import Modells
from glob import glob


def get_index(l, value):
    for i in xrange(1, len(l)):
        try:
            if l[i] <= value and l[i + 1] > value:
                return i
        except (IndexError):
            print "out of range"
            return i


class dataset:
    '''
    this class stores the data of one inputfile and parses the filename.
    it stors the applyed Force, the angles chi and omegathe data itselfe and
    the hkl with the respective angel as np.array
    '''

    def __init__(self, filename):
        self.filename = filename
        # Angles in the Labsystem
        self.Force = self.parseFilename()[0]
        self.Omega = self.parseFilename()[1]
        self.Chi = self.parseFilename()[2]
        self.data = self.readData()
        self.hkl_TTheta = np.zeros((3, 4))  # [[h, k, l, TTheta], ...]

    def parseFilename(self):
        '''
        splits the filename and returns a list with the Angles and the Force
        Returnvalue:
        Parameters = [Force, Omega, Chi]
        '''
        name_split = re.split('_', self.filename)
        force = name_split[1][0:-2]
        force = float(force)
        Omega = 0.
        Chi = 0.
        for i in xrange(1, len(name_split)):
            if (re.match('Ome', name_split[i]) != None):
                Omega = float(re.split('[a-z]+', name_split[i])[1])
            elif (re.match('Chi', name_split[i]) != None):
                Chi = float(re.split('[a-z]*', name_split[i])[1])
        Parameters = [force, Omega, Chi]
        return Parameters

    def readData(self):
        '''
        This Function reads the data stord in filename and returns the 2Thetavalue
        and the intensity in the list [TTheta, Intens]
        '''
        data = open(self.filename, "r")
        lines = data.readlines()
        TTheta = []
        Intens = []
        err = []
        for i in xrange(1, len(lines)):  #
            line = lines[i].strip()  # removes withspaces at the frond and the end
            l = re.split(r'\s*', line)
            TTheta.append(float(l[0]))
            Intens.append(float(l[1]))
            err.append(float(l[2]))
        data.close()
        return [TTheta, Intens, err]

    def print_a(self):
        return self.Omega, self.Chi, self.Force

    def print_d(self):
        print self.data

    def plot_data(self):
        plt.figure(re.split(r'\\', self.filename)[-1])
        plt.plot(self.data[0], self.data[1], 'b.', label='roh-Data')
        plt.show()

    def select_peaks(self):
        self.plot_data()
        number_o_peaks = input("select the number of peaks: ")
        while (True):
            try:
                number_o_peaks = int(number_o_peaks)
                return number_o_peaks
                break
            except (ValueError):
                print "type in only int Values"

    def select_hkl(self, auto=False, rang=[]):
        '''
            With this function the different peaks can be choosen and fittet.
            After all, the h,k,l and the coresponding 2Theta values are stored
            in the np.array "self.hkl_TTheta".
            The number_o_peaks is to type in the number of relevant peaks.
            In the second while loop you choos the datarange for the different
            h,k,l. The first loop is to do this for all peaks you want.
        '''
        number_o_peaks = 0
        if auto == False:
            number_o_peaks = self.select_peaks()
        elif auto == True:
            number_o_peaks = len(rang)

        self.hkl_TTheta = np.zeros((number_o_peaks, 5))
        # print self.hkl_TTheta
        j = 0
        while (j < number_o_peaks):
            if auto == False:
                print "%d. peak" % j
            beg = 0
            end = 10
            i = 0
            x = 0
            y = 10
            if auto == False:
                while i < 2:
                    fig = plt.figure(re.split(r'\\', self.filename)[-1])
                    ax = fig.add_subplot(111)
                    ax.plot(self.data[0], self.data[1], 'b.', label='roh-Data')
                    ax.plot(self.data[0][x:y], self.data[1][x:y], 'r-', label='roh-Data')
                    print "beg: ", beg, "end: ", end
                    cursor = DataCursor.DataCursor(ax)
                    plt.show()
                    print "select the range of the %d. peak" % j
                    beg = input("begin: ")
                    end = input("end:  ")
                    x = get_index(self.data[0], beg)
                    y = get_index(self.data[0], end)
                    i = i + 1
                print "type in the millerindices of the %d. peak" % j
                self.hkl_TTheta[j][0] = input("h: ")  # select the millerindizes
                self.hkl_TTheta[j][1] = input("k: ")
                self.hkl_TTheta[j][2] = input("l: ")
                rang.append([self.hkl_TTheta[j][0], self.hkl_TTheta[j][1], self.hkl_TTheta[j][2], x, y])
            elif auto == True:
                self.hkl_TTheta[j][0] = rang[j][0]
                self.hkl_TTheta[j][1] = rang[j][1]
                self.hkl_TTheta[j][2] = rang[j][2]
                x = rang[j][3]
                y = rang[j][4]
            dax = self.data[0][x:y]
            day = self.data[1][x:y]
            gauss = Fittingtools.gauss_lin_fitting_2(dax, day)
            # print "courffit: ", gauss
            # Fittingtools.gauss_lin_fitting_2(dax, day)
            # print "\n\n"
            # if max(day)<150.:
            #    self.hkl_TTheta[j][3]='nan'

            # else:
            # gauss = Fittingtools.gauss_lin_fitting(dax, day)

            # print Fittingtools.gauss_fit(dax, day)
            self.hkl_TTheta[j][3] = gauss[0]
            self.hkl_TTheta[j][4] = gauss[1]

            # print 'Theta: ', gauss[0],' ',self.hkl_TTheta[j][3]
            #  print "gauss", gauss
            j = j + 1
        # print self.hkl_TTheta
        if not auto:
            return rang

    # def calc_dTheta(self, latice_param, wavelength, h, k, l):
    #     '''
    #         calculate the offset of 2Theta due to misalignement of the sample
    #     '''
    #     theta_theo = r_t_d(2 * BRAGG_Eq(latice_param, wavelength, h, k, l))
    #     theta_exp = 0
    #
    #     for i in self.hkl_TTheta:
    #         if (i[0] == h and i[1] == k and i[2] == l):
    #             theta_exp = i[3]
    #
    #     if theta_exp == 0:
    #         print 'ERROR'
    #
    #     return theta_theo - theta_exp
    #
    # def write_shifted_theta_data(self, dTheta1, dTheta2, dTheta3):
    #     '''
    #         shifts the data by dTheta and stores it under a new folder and name
    #     '''
    #     filename = self.filename[0:-4]  # filename without '.eth'
    #     fil = re.split(r"\\", filename)
    #     filename = fil[0] + "\\" + fil[1] + "\\data_shifted\\" + fil[2]
    #
    #     dat = open(filename + '_shifted.eth', "w")
    #
    #     for i in xrange(len(self.data[0])):
    #         if (self.data[0][i] > 35 and self.data[0][i] < 55):
    #             dat.write(' %.3f       %i   %i\n' % (self.data[0][i] + dTheta1, self.data[1][i], self.data[2][i]))
    #         elif (self.data[0][i] > 55 and self.data[0][i] < 75):
    #             dat.write(' %.3f       %i   %i\n' % (self.data[0][i] + dTheta2, self.data[1][i], self.data[2][i]))
    #         elif (self.data[0][i] > 75 and self.data[0][i] < 90):
    #             dat.write(' %.3f       %i   %i\n' % (self.data[0][i] + dTheta3, self.data[1][i], self.data[2][i]))
    #         else:
    #             dat.write(' %.3f       %i   %i\n' % (self.data[0][i], self.data[1][i], self.data[2][i]))
    #     dat.close()
    #
    # def write_data(self):
    #     dat = open(self.filename + '_new', "w")
    #     # write header line:
    #     for i in xrange(len(self.data[0])):
    #         dat.write('%.3f,%f\n' % (self.data[0][i], self.data[1][i]))
    #     dat.close()


class one_HKL:
    """
        contains on hkl reflection
    """

    def __init__(self, h, k, l, TTheta_0, TTheta_0_weight, TTheta, TTheta_weight, chi, chi_, Omega, phi, sigma, lam):
        self.h = h
        self.k = k
        self.l = l

        # angles in the lab frame
        self.Theta_0 = TTheta_0 / 2.  # Theta unstrained
        self.Theta_0_weight = TTheta_0_weight  # uncertainty in theta_0
        self.Theta = TTheta / 2.  # Theta strained
        self.Theta_weight = TTheta_weight  # uncertainty in theta
        self.chi = chi  # chi of the sample with respect to the instrument frame
        self.chi_ = chi_  # chi of scattering vector
        self.Omega = Omega
        self.phi = phi
        self.phi_ = -90.  # 0.

        # Eulerangles in the spesimen frame
        self.psi = self.psii()  # attention th np.pi/2 - is possibly wrong
        self.PHI = self.PHII()

        # couple of other usefule constants
        self.d_epsilon = self.delta_epsilon()[0]
        self.d_epsilon_weight = self.delta_epsilon()[1]
        self.sigma = sigma
        self.d_0 = self.d_hkl(lam, self.Theta_0)
        self.d = self.d_hkl(lam, self.Theta)

    # define transformation from lab.sys to prob.sys.
    def transformation_L(self):
        """
           define transformation from lab.frame to specimen frame
           acording to Grässlin
        """
        chi = -deg_to_rad(self.chi)
        omega = deg_to_rad(self.Omega)
        phi = deg_to_rad(self.phi)

        # rot around z_L'
        O = np.array([[np.cos(phi), -np.sin(phi), 0.],
                       [np.sin(phi), np.cos(phi), 0.],
                       [0., 0., 1.]
                       ]
                      )
        # rotation around y_L' axis (lefthanded if chi<0)
        X = np.array([[np.cos(chi), 0., np.sin(chi)],
                       [0., 1., 0.],
                       [-np.sin(chi), 0., np.cos(chi)]
                       ]
                      )
        # rot around z_L
        W = np.array([[np.cos(omega), -np.sin(omega), 0.],
                      [np.sin(omega), np.cos(omega), 0.],
                      [0., 0., 1.]
                      ]
                     )
        # print self.chi, self.Omega, W*X*O
        return O.dot(X.dot(W)).transpose()

    def z_I(self):
        '''
            z of probsys in the lab frame (signe of chi already set negative)
        '''
        chi = -deg_to_rad(self.chi)
        omega = deg_to_rad(self.Omega)
        x = np.sin(chi) * np.cos(omega)
        y = np.sin(chi) * np.sin(omega)
        z = np.cos(chi)
        return np.array([x, y, z])

    def x_I(self):
        '''
            x of probsys in the lab frame
        '''
        chi = deg_to_rad(90. - self.chi)
        omega = deg_to_rad(self.Omega)
        x = np.sin(chi) * np.cos(omega)
        y = np.sin(chi) * np.sin(omega)
        z = np.cos(chi)
        return np.array([x, y, z])

    def R(self):
        '''
            transformation from Specimen frame P to
            the Measurmentframe, using the rotations {Phi +pi/2, Psi, phi}
            Gnäupel-Herold
        '''

        psi = self.psi
        Phi = self.PHI
        phi = deg_to_rad(self.phi_)
        res = np.zeros((3, 3))
        res[0, 0] = -np.cos(psi) * np.cos(Phi) * np.sin(phi) \
                    - np.sin(Phi) * np.cos(phi)
        res[0, 1] = -np.cos(psi) * np.sin(Phi) * np.sin(phi) \
                    + np.cos(Phi) * np.cos(phi)
        res[0, 2] = np.sin(psi) * np.sin(phi)

        res[1, 0] = -np.cos(psi) * np.cos(Phi) * np.cos(phi) \
                    + np.sin(Phi) * np.sin(phi)
        res[1, 1] = -np.cos(psi) * np.sin(Phi) * np.cos(phi) \
                    - np.cos(Phi) * np.sin(phi)
        res[1, 2] = np.sin(psi) * np.cos(phi)

        res[2, 0] = np.sin(psi) * np.cos(Phi)
        res[2, 1] = np.sin(psi) * np.sin(Phi)
        res[2, 2] = np.cos(psi)
        return res

    # Direction of the scatteringvector in the Lab.sys.
    def q(self):
        '''Direction of the scatteringvector in the Lab.sys.'''
        chi_ = deg_to_rad(self.chi_)
        Theta = deg_to_rad(-(90. + (self.Theta + self.Theta_0) / 2.))
        return np.array([[np.sin(chi_) * np.cos(Theta)],
                         [np.sin(chi_) * np.sin(Theta)],
                         [np.cos(chi_)]])

    def q_(self):
        '''Direction of the scatteringvector in the Lab.sys.'''
        chi_ = deg_to_rad(self.chi_)
        Theta = deg_to_rad(-(90. + (self.Theta + self.Theta_0) / 2.))
        return np.array([np.sin(chi_) * np.cos(Theta),
                         np.sin(chi_) * np.sin(Theta),
                         np.cos(chi_)])

    def LQ(self):
        '''q in specimen frame'''
        L = self.transformation_L()
        # Q = np.array(self.q())
        # print Q
        return L.dot(self.q())

    def PHII(self):
        '''Acimut angle of q in the Specimen frame'''
        lq = self.LQ()
        phi = 'nan'
        if (lq[0] > 0):
            phi = np.arctan((lq[1] / lq[0]))
        elif (lq[0] == 0):
            phi = np.sign(lq[1]) * np.pi / 2
        elif (lq[0] < 0 and lq[1] >= 0):
            phi = np.arctan((lq[1] / lq[0])) + np.pi
        elif (lq[0] < 0 and lq[1] < 0):
            phi = np.arctan((lq[1] / lq[0])) - np.pi
            # print "PHI: ",float(phi)

            # try an other way to compute this
        # phi = np.arccos(np.dot(self.x_I(), self.q_()))
        return float(phi)

    def psii(self):
        '''polar angle of q in the Specimen frame'''
        r = 1.
        psi = np.arccos(self.LQ()[2]/r)
        # psi = np.arccos(np.dot(self.z_I(), self.q_()))
        return float(psi)

    def delta_epsilon(self):
        '''
            This funktion determines the strain. epsilon = (sin(Theta0)/sin(Theta))-1
        '''
        Theta = deg_to_rad(self.Theta)
        Theta0 = deg_to_rad(self.Theta_0)
        Theta_weight = self.Theta_weight/180.*np.pi
        Theta_0_weight = self.Theta_0_weight*np.pi/180.

        e1 = np.sin(Theta0) / np.sin(Theta) - 1.
        weight = np.sqrt((np.cos(Theta0) / np.sin(Theta)*Theta_0_weight)**2 +
                         (np.sin(Theta0) / (np.sin(Theta)**2)*np.cos(Theta)*Theta_weight)**2)
        # e2=(Theta0-Theta)/np.tan(Theta)
        # print 'Epsilon: ',e1, "weight: ", weight, self.Theta_0_weight, Theta_0_weight
        # print "strain: ", e1, weight
        return e1, weight

    def d_hkl(self, lam, Theta):
        '''
            determine the laticespacing of the hkl plain
        '''
        theta = deg_to_rad(Theta)
        d = (lam * np.sqrt(self.h ** 2 + self.k ** 2 + self.l ** 2)) / (2 * np.sin(theta))
        return d

    def s1_const(self, nu):
        # assume sigma1=sigma2=0
        # res = -2*Gamma_*sigma+6*Gamma_*sigma*(np.cos(psi)**2)
        # asume sigma1=sigma2!=0
        # res = (1-2*nu)*self.sigma
        res = self.sigma
        return res

    def s2_2_const(self, nu):
        # asume sigma1=sigma2!=0
        res = np.cos(self.psi) ** 2 * self.sigma  # (np.cos(self.psi)**2-nu*np.sin(self.psi)**2)*self.sigma
        return res

    def sig_s(self):
        '''stresstensor with respect to the specimen frame
        '''
        res = np.zeros((3, 3))
        res[2, 2] = self.sigma
        return np.matrix(res)

    def sig_q(self, nu=None):
        sigma = np.zeros((3, 3))
        sigma[2, 2] = self.sigma
        sigma[1, 1] = -self.sigma * nu
        sigma[0, 0] = -self.sigma * nu
        if nu is None:
            return self.LQ()[2] + self.sigma

    def sig_L(self):
        '''calculates the stresstensor with respect to the
           mesurmentframe
        '''
        s = np.matrix(np.array([[np.sin(self.psi) ** 2, 0., -np.sin(self.psi) * np.cos(self.psi)],
                                [0., 0., 0.],
                                [-np.sin(self.psi) * np.cos(self.psi), 0., np.cos(self.psi) ** 2]]) * self.sigma)
        s = self.R() * self.R() * self.sig_s()
        return self.sigma * np.cos(self.psi) ** 2


class Data:
    def __init__(self, path, odf_name, diameter):
        self.__hkl_object_list = []
        self.__path_to_data = path
        self.__odf = Modells.ODF()
        self.__set_odf(odf_name)
        self.__sample_spezifikations = {"diameter": diameter}
        self.__crystal_sym = {"cubic": "cubic"}
        self.__epsilon_list = []
        self.__epsilon_weight_list = []
        self.__hkl_phi_psi_list = []  # hkl, phi, psi

    """
    Deal with the ODF
    """

    def __set_odf(self, filename):
        self.__odf.read_data(self.__path_to_data, filename)

    def get_path_to_data(self):
        return self.__path_to_data

    def plot_odf(self):
        self.__odf.plot_odf()

    def f_odf(self, phi1, Phi, phi2):
        """
        :param phi1: in deg
        :param Phi:  in deg
        :param phi2: in deg
        :return: intensity of the ODF at phi1, Phi, phi2
        """
        return self.__odf.f(phi1, Phi, phi2)

    @property
    def integral_over_total_odf(self):
        return self.__odf.integral_over_all_orientations_g

    def calc__deltavals(self):
        lis = []
        step = 15
        for i in np.arange(0, 361, step):
            h, k, l = 1, 1, 0
            psi = 0
            phi = 0
            rot = i
            # print "phi: %i, psi: %i"%(phi, psi)
            lis.append(self.__odf.calc_delta_vals(h, k, l, psi, phi, rot, step))
        # print lis
        # for i in lis:
        #     print "rot: %i, dphi1: %.1f, dphi: %.1f, dphi2: %.1f, phi1: %.1f, phi: %.1f, phi2: %.1f"%(i[1], i[0][0], i[0][1], i[0][2], i[2][0], i[2][1], i[2][2])
        # plt.show()
        '''
        res = []
        e_angles = []
        step = 5
        ran = np.arange(0, 361, step)
        for i in ran:
            res.append(self.__odf.calc_delta_vals(h, k, l, psi, phi, i, i+step))


        for i in res:
            print "dvals: %.1f, %.1f, %.1f \t\t rot-angle: %.1f \t euler-angles: %.1f, %.1f, %.1f, %.1f" % (i[0][0], i[0][1], i[0][2], i[1], i[2][0], i[2][1], i[2][2], i[2][0]+i[2][2])
        print len(e_angles), "  ", len(res)
        # print "dvals: ", dvals
        '''

    """
    Fit the data
    """

    def Fit_the_data_with_texture(self, method, filename, number_of_datapoints = None):
        print "Number of datapoints: ", len(self.__epsilon_list), method
        fit = Modells.Fit_strain_with_texture(odf1=self.__odf, odf2=None,
                                              force=self.__sample_spezifikations["force"],
                                              diameter=self.__sample_spezifikations["diameter"],
                                              strains_data=self.__epsilon_list[0:number_of_datapoints],
                                              xvals=self.__hkl_phi_psi_list[0:number_of_datapoints],
                                              weights=self.__epsilon_weight_list[0:number_of_datapoints])  #
        filename = filename + method
        for i in range(len(self.__epsilon_list)):
            print self.__hkl_phi_psi_list[i], self.__epsilon_list[i]
        fit.do_the_fitting(filename=filename, material="Iron", method=method)

    """
    Read the scattering data and process it
    """

    def read_scattering_data(self, path_of_straind_data, path_of_unstraind_data):
        '''
        :param path_of_straind_data: (e.g. 'Euler-Scans unter 5kN\\')
        :param path_of_unstraind_data: (e.g. 'Euler-Scans ohne Last\\')
        :return: void
        '''
        unstraind = self.__path_to_data + path_of_unstraind_data + "*.eth"
        straind = self.__path_to_data + path_of_straind_data + "*.eth"
        filelist_unstraind = glob(unstraind)
        filelist_straind = glob(straind)
        filelist_unstraind.sort()
        filelist_straind.sort()
        unstraind_data_object_list = self.__creat_data_object_list(filelist_unstraind)
        straind_data_object_list = self.__creat_data_object_list(filelist_straind)
        self.__select_peaks(unstraind_data_object_list, straind_data_object_list)
        self.__hkl_object_list = self.__create_hkl_object_list(unstraind_data_object_list, straind_data_object_list)
        # self.__create_epsilon_list()

    def __creat_data_object_list(self, filelist):
        # list of dataset-objekts
        list = []
        for i in filelist:
            list.append(dataset(i))
        return list

    def __select_peaks(self, unstraind, straind):
        '''
        -------------------------------------------------------------------------------------------
        select the peaks and calculate 2Theta
        hkl_setting is a list of lists. The elements of hkl_setting contain
        - the hkl (first three indiices)
        - the index of the minimal and maximal value of the respectiv peak.

        with hkl_setting = unstraind[0].select_hkl() it is possible to select the peaks
        with some kinde of user interface
        stand 08.12.2015:
        the user interface is not vary sophisticated
        -------------------------------------------------------------------------------------------
        '''
        # 'normal' Iron
        # hkl setting 3 peaks [[1.0, 1.0, 0.0, 780, 980], [2.0, 0.0, 0.0, 1180, 1400], [2.0, 1.0, 1.0, 1540, 1740]]
        # hkl_setting 5 peaks [[1.0, 1.0, 0.0, 852, 904], [2.0, 0.0, 0.0, 1260, 1314], [2.0, 1.0, 1.0, 1600, 1676], [2.0, 2.0, 0.0, 1926, 2020], [3.0, 1.0, 0.0, 2266, 2380]]
        hkl_setting = [[1.0, 1.0, 0.0, 852, 904], [2.0, 0.0, 0.0, 1260, 1314], [2.0, 1.0, 1.0, 1600, 1676],
                       [2.0, 2.0, 0.0, 1926, 2020], [3.0, 1.0, 0.0, 2266, 2380]]
        print hkl_setting
        # do it manually:
        # hkl_setting = unstraind[0].select_hkl(auto=False, rang=[])

        # do it for the rest automatically:
        for i in xrange(0, len(unstraind)):
            unstraind[i].select_hkl(auto=True, rang=hkl_setting)
            straind[i].select_hkl(auto=True, rang=hkl_setting)

    def __create_hkl_object_list(self, unstraind_data_object_list, straind_data_object_list):
        '''
        -------------------------------------------------------------------------------------------
        store objekt's, containing one hkl each, the angles with respect to the labframe
        and the eulerangles psi and PHI with respekct to the specimen frame.
        -------------------------------------------------------------------------------------------
        '''
        self.__sample_spezifikations["force"] = straind_data_object_list[0].Force
        # print straind_data_object_list[0]
        # print "FORCE: ", self.__sample_spezifikations["force"]
        Force = self.__sample_spezifikations["force"]
        diameter = self.__sample_spezifikations["diameter"]
        nu = 0  # .30
        HKL_object_list = []
        for i in xrange(0, len(unstraind_data_object_list[0].hkl_TTheta)):  # loop over all measured hkl
            a = []
            # print unstraind[0].hkl_TTheta[i]
            for j in xrange(0, len(unstraind_data_object_list)):  # loop over all datapoints
                sigma = sigma3(Force, diameter)  # stress
                chi = unstraind_data_object_list[j].Chi  # chi of j'th datapoint with respect to labsys
                chi_ = 90.  # chi of scattering vector
                Omega = unstraind_data_object_list[j].Omega  # omega of the j'th datapoint
                TTheta_0 = unstraind_data_object_list[j].hkl_TTheta[i][3]
                TTheta_0_weight = unstraind_data_object_list[j].hkl_TTheta[i][4]
                      # theta_0 of the j'th datapoint and the i'th hkl
                TTheta = straind_data_object_list[j].hkl_TTheta[i][3]  # theta of the j'th datapoint and the i'th hkl
                TTheta_weight = straind_data_object_list[j].hkl_TTheta[i][4]
                phi = -90.  # phi rot of the specimen

                h = unstraind_data_object_list[j].hkl_TTheta[i][0]
                k = unstraind_data_object_list[j].hkl_TTheta[i][1]
                l = unstraind_data_object_list[j].hkl_TTheta[i][2]

                lam = 0.154840  # wavelength of SPODI
                # print "hkl: ", h, k, l, TTheta_0, TTheta

                if math.isnan(TTheta_0) or math.isnan(TTheta):
                    pass
                else:
                    b = one_HKL(h, k, l, TTheta_0,TTheta_0_weight, TTheta,TTheta_weight, chi, chi_, Omega, phi, sigma,
                                lam)
                    if math.isnan(b.psi) != True and math.isnan(b.PHI) != True:
                        print  "%i%i%i\tTheta: %.3f\tchi: %.0f\tomega: %.0f\tPsi: %.3f\tPhi: %.3f"% \
                                (h, k, l, TTheta_0 / 2, chi, Omega, r_t_d(b.psi), r_t_d(b.PHI))
                        a.append(b)

            HKL_object_list += a

        print HKL_object_list

        for i in HKL_object_list:
            # if i.chi != 90:
            self.__epsilon_list.append(i.d_epsilon)
            self.__epsilon_weight_list.append(i.d_epsilon_weight)
            h = i.h
            k = i.k
            l = i.l

            # print [h, k, l, j.PHI, j.psi, j.d_epsilon]
            self.__hkl_phi_psi_list.append([i.PHI, i.psi, i.h, i.k, i.l])
        return HKL_object_list


'''
def delta_epsilon(Theta, Theta0):
'''
# '''
# This funktion determines the strain. epsilon = (sin(Theta0)/sin(Theta))-1
# '''
'''
    Theta = deg_to_rad(Theta)
    Theta0= deg_to_rad(Theta0)
    e1=np.sin(Theta0)/np.sin(Theta)-1.
    e2=(Theta0-Theta)/np.tan(Theta)
    #print 'Epsilon: ',e1, e2, e2-e1
    return e1
'''


def sigma3(Force  # in N
           , diameter  # in mm
           ):
    '''
        This function determines the 33-component of the macro-straintensor due
        to the Force in this (and only in this) direction with respekt to
        the specimen reference frame.
        The component is:
            sigma33 = F/A
        where A is the crosssection of the Specimen
    '''
    A = (diameter * np.power(10., -3.)) ** 2 / 4 * np.pi
    return Force / A  # *np.power(10.,-9)


def deg_to_rad(deg):
    return np.pi * deg / 180.


def r_t_d(rad):
    return 180. * rad / np.pi


# just some useful coeficients
def Gama(h, k, l):
    g = ((h ** 2) * (k ** 2) + (h ** 2) * (l ** 2) + (k ** 2) * (l ** 2)) / (((h ** 2) + (k ** 2) + (l ** 2)) ** 2)
    return g


def BRAGG_Eq(laticeparam, wavelength, h, k, l):
    '''
        caululates the Theta angle of the refracting plane.
        !!!Laticeparam and wavelength must be of the same unit!!!
    '''
    return np.arcsin(wavelength * np.sqrt(h ** 2 + k ** 2 + l ** 2) / (2 * laticeparam))


'''
def psii(chi, chi_, omega, Theta, phi):
    chi = deg_to_rad(chi)
    chi_ = deg_to_rad(chi_)
    omega = deg_to_rad(omega)
    Theta = deg_to_rad(Theta)
    phi = deg_to_rad(phi)

    x = (np.cos(phi)*np.cos(chi)*np.cos(omega)-np.sin(phi)*np.sin(omega))*np.sin(chi_)*np.cos(Theta) \
        -(np.cos(phi)*np.cos(chi)*np.sin(omega)+np.sin(phi)*np.cos(omega))*np.sin(chi_)*np.sin(Theta) \
        +np.cos(phi)*np.sin(chi)*np.cos(chi_)

    y = (np.sin(phi)*np.cos(chi)*np.cos(omega)+np.cos(phi)*np.sin(omega))*np.sin(chi_)*np.cos(Theta) \
        +(-np.sin(phi)*np.cos(chi)*np.sin(omega)+np.cos(phi)*np.cos(omega))*np.sin(chi_)*np.sin(Theta) \
        +np.sin(phi)*np.sin(chi)*np.cos(chi_)

    z = -np.sin(chi)*np.cos(omega)*np.sin(chi_)*np.cos(Theta) \
        +np.sin(chi)*np.sin(omega)*np.sin(chi_)*np.sin(Theta) \
        +np.cos(chi)*np.cos(chi_)
    #print "xyz: ",x,y,z
    r = np.sqrt(x**2+y**2+z**2)
    #print 'r: ', r
    return np.arccos(z/r)


def alpha(chi, chi_, omega, Theta, phi,h,k,l,Force, diameter, nu):
    psi = psii(chi, chi_, omega, Theta, phi)
    Gamma_ = Gama(h,k,l)
    sigma = sigma3(Force, diameter)
    #print "psii", psi
    #assume sigma1=sigma2=0
    #res = (1-Gamma_)*sigma+(3*Gamma_-1)*sigma*(np.cos(psi)**2)
    #asume sigma1=sigma2!=0
    res = ((1-Gamma_)*(1-2*nu)+(3*Gamma_-1)*(np.cos(psi)**2-nu*np.sin(psi)**2))*sigma
    return res

def beta(chi, chi_, omega, Theta, phi,h,k,l,Force, diameter,nu):
    psi = psii(chi, chi_, omega, Theta, phi)
    Gamma_ = Gama(h,k,l)
    sigma = sigma3(Force, diameter)
    #assume sigma1=sigma2=0
    #res = Gamma_*sigma+(1-3*Gamma_)*sigma*(np.cos(psi)**2)
    #asume sigma1=sigma2!=0
    res = (Gamma_*(1-2*nu)+(1-3*Gamma_)*(np.cos(psi)**2-nu*np.sin(psi)**2))*sigma
    return res

def gamma(chi, chi_, omega, Theta, phi,h,k,l,Force, diameter,nu):
    psi = psii(chi, chi_, omega, Theta, phi)
    Gamma_ = Gama(h,k,l)
    sigma = sigma3(Force, diameter)
    #print "psii", psi
    #print "hkl, gamma: ", h,k,l,Gamma_
    #assume sigma1=sigma2=0
    #res = -2*Gamma_*sigma+6*Gamma_*sigma*(np.cos(psi)**2)
    #asume sigma1=sigma2!=0
    res = (-2*Gamma_*(1-2*nu)+6*Gamma_*(np.cos(psi)**2-nu*np.sin(psi)**2))*sigma
    return res

def s1_const(chi, chi_, omega, Theta, phi,h,k,l,Force, diameter,nu):
    #psi = psii(chi, chi_, omega, Theta, phi)
    #Gamma_ = Gama(h,k,l)
    sigma = sigma3(Force, diameter)
    #print "hkl, gamma: ", h,k,l,Gamma_
    #assume sigma1=sigma2=0
    #res = -2*Gamma_*sigma+6*Gamma_*sigma*(np.cos(psi)**2)
    #asume sigma1=sigma2!=0
    res = (1-2*nu)*sigma
    return res

def s2_2_const(chi, chi_, omega, Theta, phi,h,k,l,Force, diameter,nu):
    psi = psii(chi, chi_, omega, Theta, phi)
    #Gamma_ = Gama(h,k,l)
    sigma = sigma3(Force, diameter)
    #print "hkl, gamma: ", h,k,l,Gamma_
    #assume sigma1=sigma2=0
    #res = -2*Gamma_*sigma+6*Gamma_*sigma*(np.cos(psi)**2)
    #asume sigma1=sigma2!=0
    res = (np.cos(psi)**2-nu*np.sin(psi)**2)*sigma
    return res
'''
