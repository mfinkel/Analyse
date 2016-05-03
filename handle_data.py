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
        self.fitted_data = DataContainer()  # contains the phi_psi_hkl and the strain_stress list for all phases and
        # applied forces
        self.odf_phase_1 = self.get_odf(odf_phase_1_file)
        self.odf_phase_2 = self.get_odf(odf_phase_2_file)
        self.sample_diameter = sample_diameter

    @staticmethod
    def get_odf(odf_file):
        """
        this function reads the odf from the odf_file and creates a ODF() object.
        :param odf_file:
        :return:
        """
        if odf_file is None:
            return None
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
        area = (self.sample_diameter * np.power(10., -3.)) ** 2 / 4 * np.pi
        stress = force * np.power(10., 3) / area
        stress_error = stress * (2 * 0.01 / self.sample_diameter + 0.05)
        return stress, stress_error


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
        """
        Fit all data and 
        :param peak_regions_phase: is a list: [[phase_1, peak_region], [phase_2, peak_region]]
        :param plot:
        :return:
        """
        for i in peak_regions_phase:  # loop over all phases
            phase, peak_regions = i
            # self.data_dic_phases[phase] = list()
            force_dic = {}
            for j in self.data_dic_raw:  # loop over all forces
                help_list = []
                for n in self.data_dic_raw[j]:  # loop over all orientations of the sample
                    force, omega, chi, two_theta, intens, error = n
                    data = [two_theta, intens, error]
                    hkl_2_theta = self.fit_the_peaks_for_on_diffraktin_pattern(data=data, peak_regions=peak_regions)
                    hkl_2_theta = hkl_2_theta.tolist()
                    for m in hkl_2_theta:
                        help_list.append([phase, force, omega, chi, m])
                force_dic[j] = help_list  # [[phase, force, omega, chi, h, k, l, 2theta, 2theta_error], ...]
            self.data_dic_phases[phase] = force_dic

        # calculate phi, psi, strain and stress and store it in
        self.calc_phi_psi_epsilon()

    def calc_phi_psi_epsilon(self):
        """
        calculate phi and psi angles of the scattering direction with respect to the specimen frame using the
        angles chi and omega, and calculating the strain
        :return: None
        """
        for i, force_dict in self.data_dic_phases:  # loop over all phases, i is the key of the dic
            key_list = sorted(force_dict.keys())  # key list of the force dict
            self.fitted_data.data_dict[i] = []
            # self.data_dict[i] = []

            force_stress_dict = {}
            for j, v in enumerate(key_list):
                # loop over all forces, v is the key of the force dic (it is equal to the force)
                # j is the position of the key, j+1 is the next force, j-1 is the previous force
                phi_psi_hkl = []
                strain_stress = []
                if v == 0:
                    pass
                else:
                    for n in xrange(len(force_dict[v])):  # loop over all data of each force
                        # values of the straind data:
                        phase, force, omega, chi, hkl_2_theta = force_dict[v][n]
                        h, k, l, two_theta, two_theta_error = hkl_2_theta

                        # vals of the unstraind data
                        phase_0, force_0, omega_0, chi_0, hkl_2_theta_0 = force_dict[key_list[0]][n]
                        h_0, k_0, l_0, two_theta_0, two_theta_error_0 = hkl_2_theta_0

                        phi = self.PHII(chi_of_scatteringvector=90, theta=two_theta / 2, theta_o=two_theta_0 / 2,
                                        chi=chi, omega=omega)

                        psi = self.psii(chi_of_scatteringvector=90, theta=two_theta / 2, theta_o=two_theta_0 / 2,
                                        chi=chi, omega=omega)

                        phi_psi_hkl.append([phi, psi, h, k, l])

                        strain, strain_error = self.delta_epsilon(two_theta=two_theta,
                                                                  two_theta_0=two_theta_0,
                                                                  two_theta_err=two_theta_error,
                                                                  two_theta_0_err=two_theta_error_0)

                        stress, stress_err = self.calc_applied_stress(force=force)
                        strain_stress.append([strain, strain_error, stress, stress_err])
                    force_stress_dict[v] = [phi_psi_hkl, strain_stress]
            self.fitted_data.data_dict[i] = force_stress_dict

    @staticmethod
    def delta_epsilon(two_theta, two_theta_0, two_theta_err, two_theta_0_err):
        """
            This method determines the strain. epsilon = (sin(Theta0)/sin(Theta))-1
        """
        Theta = deg_to_rad(two_theta)
        Theta0 = deg_to_rad(two_theta_0)
        Theta_weight = two_theta_err / 180. * np.pi / 2
        Theta_0_weight = two_theta_0_err * np.pi / 180. / 2

        e1 = np.sin(Theta0) / np.sin(Theta) - 1.

        weight = np.sqrt((np.cos(Theta0) / np.sin(Theta) * Theta_0_weight) ** 2 +
                         (np.sin(Theta0) / (np.sin(Theta) ** 2) * np.cos(Theta) * Theta_weight) ** 2)
        # e2=(Theta0-Theta)/np.tan(Theta)
        # print 'Epsilon: ',e1, "weight: ", weight, self.Theta_0_weight, Theta_0_weight
        # print "strain: ", e1, weight
        return e1, weight

    @staticmethod
    def transformation_L_from_I_to_P(chi, omega):
        """
           define transformation from lab.frame I to specimen frame P
           acording to Graesslin
        """
        chi = deg_to_rad(chi)
        omega = -deg_to_rad(omega)
        phi = 0

        # rot around z_L'
        O = np.array([[np.cos(phi), np.sin(phi), 0.],
                      [-np.sin(phi), np.cos(phi), 0.],
                      [0., 0., 1.]
                      ]
                     )
        # rotation around y_L' axis (lefthanded if chi<0)
        X = np.array([[np.cos(chi), 0., -np.sin(chi)],
                      [0., 1., 0.],
                      [np.sin(chi), 0., np.cos(chi)]
                      ]
                     )
        # rot around z_L
        W = np.array([[np.cos(omega), np.sin(omega), 0.],
                      [-np.sin(omega), np.cos(omega), 0.],
                      [0., 0., 1.]
                      ]
                     )
        res = W.dot(X.dot(O))
        # res = O.dot(X.dot(W))
        # L_1 = np.dot(res, np.array([[1], [0], [0]]))
        # L_2 = np.dot(res, np.array([[0], [1], [0]]))
        # L_3 = np.dot(res, np.array([[0], [0], [1]]))
        # titel = "chi: {}, omega: {}, phi: {}".format(r_t_d(chi), r_t_d(omega), r_t_d(phi))
        # cplot.plot_coordinatframe(L_1, L_2, L_3, Q=self.q(), titel=titel)
        # plt.show()
        # # print self.chi, self.Omega, W*X*O
        # print "PSI: ", r_t_d(np.arccos(np.dot(L_3.transpose(), self.q())))
        return res  # O.dot(X.dot(W))  # .transpose()

    def z_P_in_lab_frame(self, chi, omega):
        '''
            z of probsys in the lab frame
        '''
        # chi = -deg_to_rad(chi)
        # omega = deg_to_rad(Omega)
        # x = np.sin(chi) * np.cos(omega)
        # y = np.sin(chi) * np.sin(omega)
        # z = np.cos(chi)
        return np.dot(self.transformation_L_from_I_to_P(chi, omega), np.array([[0], [0], [1]])).transpose()

    @staticmethod
    def q(chi_of_scatteringvector, theta, theta_o):
        """Direction of the scattering vector in the Lab.sys."""
        chi_ = deg_to_rad(chi_of_scatteringvector)
        Theta = deg_to_rad(-(90. + (theta + theta_o) / 2.))
        q = np.array([[np.sin(chi_) * np.cos(Theta)],
                      [np.sin(chi_) * np.sin(Theta)],
                      [np.cos(chi_)]])
        # q=np.array([[-1/np.sqrt(2)], [-1/np.sqrt(2)], [0]])
        return q

    def LQ(self, chi_of_scatteringvector, theta, theta_o, chi, omega):
        '''q in specimen frame'''
        L = self.transformation_L_from_I_to_P(chi, omega).transpose()
        # Q = np.array(self.q())
        # print Q
        return L.dot(self.q(chi_of_scatteringvector, theta, theta_o))

    def psii(self, chi_of_scatteringvector, theta, theta_o, chi, omega):
        '''polar angle of q in the Specimen frame'''
        r = 1.
        # psi = np.arccos(self.LQ()[2] / r)
        psi = np.arccos(np.dot(self.z_P_in_lab_frame(chi, omega), self.q(chi_of_scatteringvector, theta, theta_o)))
        return float(psi)

    def PHII(self, chi_of_scatteringvector, theta, theta_o, chi, omega):
        '''Acimut angle of q in the Specimen frame'''
        lq = self.LQ(chi_of_scatteringvector, theta, theta_o, chi, omega)
        # print "lq", lq
        # lq = self.LQ()
        phi = 'nan'
        if (lq[0] > 0):
            phi = np.arctan((lq[1] / lq[0]))
        elif (lq[0] == 0):
            phi = np.sign(lq[1]) * np.pi / 2
        elif lq[0] < 0 and lq[1] >= 0:
            phi = np.arctan((lq[1] / lq[0])) + np.pi
        elif lq[0] < 0 and lq[1] < 0:
            phi = np.arctan((lq[1] / lq[0])) - np.pi
            # print "PHI: ",float(phi)

            # try an other way to compute this
        # phi = np.arccos(np.dot(self.x_I(), self.q_()))
        # print "phi: ", r_t_d(phi)
        return float(phi)


class DataContainer(object):
    def __init__(self):
        self.data_dict = {}  # dictionary of the data. Key is the phase. val is the list [hkl_phi_psi, epsilon_sigma]

    def get_data_phase_1_force(self, force):
        return self.data_dict[1][force]

    def get_data_phase_2_force(self, force):
        return self.data_dict[2][force]


def deg_to_rad(deg):
    return np.pi * deg / 180.


def r_t_d(rad):
    return 180. * rad / np.pi