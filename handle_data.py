"""
This contains the classstructure to treat the Data sets.

The class data is the base class. It has an
- epsilon list
- an hkl phi psi list
- a method to calc the macrostress
- a method to set the path

The other classes inherit form this class

"""
import re
import numpy as np
import Fittingtools
import Modells
import os
from glob import glob


class Data(object):
    def __init__(self, sample_diameter, odf_phase_1_file=None, odf_phase_2_file=None):
        self.data_list = []
        self.odf_phase_1 = self.get_odf(odf_phase_1_file)
        self.odf_phase_2 = self.get_odf(odf_phase_2_file)
        self.sample_diameter = sample_diameter

    @staticmethod
    def get_odf(odf_file):
        odf = Modells.ODF()
        odf_file = os.path.normpath(odf_file)
        path, filename = os.path.split(odf_file)
        path = "{0}\\".format(path)
        odf.read_data(path=path, filename=filename)
        return odf

    def calc_applied_stress(self, force):
        """
        :param force: applied force in kN
        :return: stress
        """
        area = (self.diameter * np.power(10., -3.)) ** 2 / 4 * np.pi
        return force * np.power(10., 3) / area  # *np.power(10.,-9)


class SPODIData(Data):
    def __init__(self, sample_diameter):
        super(SPODIData, self).__init__(sample_diameter)
        self.data_dic_raw = {}  # a dictionary containing the raw data, the key is the applied force, the values
        # are lists of all measured orientation angles with the measured two_theta, intens, error data
        self.data_dic_phases = {}  # dictionary containing dictionaries of all peaks of the different phases
        # data_dic_phases key=phase, val=dic containing all reflexes (key = force)

    def load_data(self, list_of_dirs):
        """
        loading the datafiles in the folders in the list "load the data"
        :param list_of_dirs:
        :return:
        """
        if type(list_of_dirs) is not list:
            list_of_dirs = [list_of_dirs]

        for i in list_of_dirs:
            data_files = i + "*.eth"
            data_files = glob(data_files)
            data_files.sort()
            help_list = []
            force = 0
            for j in data_files:
                # force in kN, omega and Chi in deg
                force, omega, chi = self.parse_filename(j)
                # read the data from the file
                two_theta, intens, error = self.read_data(j)
                help_list.append([force, omega, chi, two_theta, intens, error])
                # force, omega, Chi, two theta, intensity, error of intensity
            self.data_dic_raw[force] = help_list

    @staticmethod
    def parse_filename(filename):
        '''
        splits the filename and returns a list with the Angles and the Force
        Returnvalue:
        Parameters = [Force, Omega, Chi]
        :param filename:
        '''
        name_split = re.split('_', filename)
        force = name_split[1][0:-2]
        force = float(force)
        omega = 0.
        Chi = 0.
        for i in xrange(1, len(name_split)):
            if re.match('Ome', name_split[i]) is not None:
                omega = float(re.split('[a-z]+', name_split[i])[1])
            elif re.match('Chi', name_split[i]) is not None:
                Chi = float(re.split('[a-z]*', name_split[i])[1])
        Parameters = [force, omega, Chi]
        return Parameters

    @staticmethod
    def read_data(filename):
        '''
        This Function reads the data stord in filename and returns the 2Thetavalue
        and the intensity and the error in the list [TTheta, Intens, err]
        :param filename:
        '''
        data = open(filename, "r")
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

    def fit_the_peaks_for_on_diffraktin_pattern(self, data, peak_regions, plot=False):
        """
        fitting the indiwidual peaks of each difraction pattern
        :param peak_regions:
        :return:
        """
        # hkl_2_theta = [[h, k, l, 2Theta, 2Theta_error, phase], [], ...
        hkl_2_theta = np.zeros((len(peak_regions), 5))
        for j in xrange(len(peak_regions)):
            hkl_2_theta[j][0] = peak_regions[j][0]
            hkl_2_theta[j][1] = peak_regions[j][1]
            hkl_2_theta[j][2] = peak_regions[j][2]
            x = int(peak_regions[j][3])
            y = int(peak_regions[j][4])
            dax = data[0][x:y]
            day = data[1][x:y]
            day_err = data[2][x:y]
            gauss = Fittingtools.gauss_lin_fitting_2(dax, day, day_err, plot=plot)
            hkl_2_theta[j][3] = gauss[0]
            hkl_2_theta[j][4] = gauss[1]
        return hkl_2_theta

    def fit_all_data(self, peak_regions_phase, plot=False):
        # peak_region_phase is a list: [[phase_1, peak_region], [phase_2, peak_region]]

        for i in peak_regions_phase:  # loop over all phases
            phase, peak_regions = i
            # self.data_dic_phases[phase] = list()
            force_dic = {}
            for j in self.data_dic_raw:  # loop over all forces
                help_list = []
                for n in self.data_dic_raw[j]:  # loop over all orientations
                    force, omega, chi, two_theta, intens, error = n
                    data = [two_theta, intens, error]
                    hkl_2_theta = self.fit_the_peaks_for_on_diffraktin_pattern(data=data, peak_regions=peak_regions)
                    hkl_2_theta = hkl_2_theta.tolist()
                    for m in hkl_2_theta:
                        help_list.append([phase, force, omega, chi, m])
                force_dic[j] = help_list  # [[phase, force, omega, chi, h, k, l, 2theta, 2theta_error], ...]
            self.data_dic_phases[phase] = force_dic


class DataContainer(object):
    def __init__(self):
        self.phi_psi_hkl_list = []
        self.epsilon_list = []
        self.epsilon_weight_list = []
        self.stress = 0.

    def insert_data_point(self, phi, psi, h, k, l, epsilon, delta_epsilon, stress):
        self.phi_psi_hkl_list.append([phi, psi, h, k, l])
        self.epsilon_list.append(epsilon)
        self.epsilon_weight_list.append(delta_epsilon)
        self.stress = stress
