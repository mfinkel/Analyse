from __future__ import print_function
import functools
import sys
import os
import lmfit as lm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import matplotlib as mpl
from matplotlib.backends import qt4_compat
import matplotlibwidget

use_pyside = qt4_compat.QT_API == qt4_compat.QT_API_PYSIDE

if use_pyside:
    from PySide.QtCore import *
    from PySide.QtGui import *
else:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *

from functools import partial
import handle_data
import Modells
import thread
from multiprocessing import Process, Queue  # processing ,freezeSupport
from time import sleep


class insert_const_widget(QWidget):
    def __init__(self, *args):
        QWidget.__init__(self)
        QWidget.resize(self, 500, 400)
        QWidget.setWindowTitle(self, "default")
        self.name = "default"
        sym = "default"
        self.layout = QGridLayout()

        self.params = lm.Parameters()
        self.add_params(self.params, sym)
        self.sym = sym
        self.create_text_output(self, self.sym + ":", 20, 20, 100, 30)
        self.array = self.create_arry()
        self.set_preprepare_array()
        # self.exitButton = QAction(QIcon('exit24.png'), 'Exit', self)
        # self.exitButton.setShortcut('Ctrl+Q')
        # self.exitButton.setStatusTip('Exit application')
        # self.exitButton.triggered.connect(self.exit)
        # self.exitButton.triggered.connect(self.close)

        self.confirm_button = QPushButton("confirm", self)
        self.confirm_button.move(140, 20)
        self.send_and_close_button = QPushButton("send and close", self)
        self.send_and_close_button.move(240, 20)
        self.send_and_close_button.clicked.connect(self.send_and_close)
        self.confirm_button.clicked.connect(self.set_params)

        # self.close()

    def show_of(self, name, sym):
        self.name = name
        self.sym = sym
        self.set_preprepare_array()
        self.show()

    @staticmethod
    def add_params(params, sym):
        if sym == "isotope":
            params.add('c_11', value=220 * np.power(10., 9), min=0 * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_12', value=126 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
        elif sym == "m-3m":
            params.add('c_11', value=240 * np.power(10., 9), min=10. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_12', value=120 * np.power(10., 9), min=10. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_44', value=120 * np.power(10., 9), min=10. * np.power(10., 9), max=600. * np.power(10., 9))
        elif sym == "hexagonal":
            params.add('c_11', value=217 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_12', value=120 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_13', value=120 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_33', value=126 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))
            params.add('c_44', value=126 * np.power(10., 9), min=0. * np.power(10., 9), max=600. * np.power(10., 9))

    def create_arry(self):
        x, y = 50, 30
        text_obj = []
        for j in xrange(7):
            l = []
            for i in xrange(7):
                if i == j == 0:
                    pass
                elif i == 0 or j == 0:
                    obj = self.create_text_output(self, "{}".format(i | j), 20 + i * (x + 10), 60 + j * (y + 10), x, y)
                    l.append(obj)
                else:
                    obj = self.create_text_output(self, "c_{1}{0}".format(i, j), 20 + i * (x + 10), 60 + j * (y + 10),
                                                  x, y)
                    l.append(obj)
            text_obj.append(l)
        return text_obj

    def set_params(self):
        if self.sym == "isotrop":
            try:
                self.params['c_11'].value = float(self.array[1][1].toPlainText()) * np.power(10., 9.)
                self.params['c_12'].value = float(self.array[1][2].toPlainText()) * np.power(10., 9.)
            except KeyError:
                self.params.add('c_11', value=float(self.array[1][1].toPlainText()) * np.power(10., 9.))
                self.params.add('c_12', value=float(self.array[1][2].toPlainText()) * np.power(10., 9.))

            c_11 = float(self.array[1][1].toPlainText())
            c_12 = float(self.array[1][2].toPlainText())
            c_44 = 0.5 * (c_11 - c_12)
            self.array[2][2].setText(str(c_11))
            self.array[3][3].setText(str(c_11))
            self.array[1][3].setText(str(c_12))
            self.array[2][3].setText(str(c_12))
            self.array[4][4].setText(str(c_44))
            self.array[5][5].setText(str(c_44))
            self.array[6][6].setText(str(c_44))

        elif self.sym == "m-3m":
            try:
                self.params['c_11'].value = float(self.array[1][1].toPlainText()) * np.power(10., 9.)
                self.params['c_12'].value = float(self.array[1][2].toPlainText()) * np.power(10., 9.)
                self.params['c_44'].value = float(self.array[4][4].toPlainText()) * np.power(10., 9.)
            except KeyError:
                self.params.add('c_11', value=float(self.array[1][1].toPlainText()) * np.power(10., 9.))
                self.params.add('c_12', value=float(self.array[1][2].toPlainText()) * np.power(10., 9.))
                self.params.add('c_44', value=float(self.array[4][4].toPlainText()) * np.power(10., 9.))

            c_11 = float(self.array[1][1].toPlainText())
            c_12 = float(self.array[1][2].toPlainText())
            c_44 = float(self.array[4][4].toPlainText())
            self.array[2][2].setText(str(c_11))
            self.array[3][3].setText(str(c_11))
            self.array[1][3].setText(str(c_12))
            self.array[2][3].setText(str(c_12))
            self.array[5][5].setText(str(c_44))
            self.array[6][6].setText(str(c_44))

        elif self.sym == "hexagonal":
            try:
                self.params['c_11'].value = float(self.array[1][1].toPlainText()) * np.power(10., 9.)
                self.params['c_12'].value = float(self.array[1][2].toPlainText()) * np.power(10., 9.)
                self.params['c_13'].value = float(self.array[1][3].toPlainText()) * np.power(10., 9.)
                self.params['c_33'].value = float(self.array[3][3].toPlainText()) * np.power(10., 9.)
                self.params['c_44'].value = float(self.array[4][4].toPlainText()) * np.power(10., 9.)
            except KeyError:
                self.params.add('c_11', value=float(self.array[1][1].toPlainText()) * np.power(10., 9.))
                self.params.add('c_12', value=float(self.array[1][2].toPlainText()) * np.power(10., 9.))
                self.params.add('c_13', value=float(self.array[1][3].toPlainText()) * np.power(10., 9.))
                self.params.add('c_33', value=float(self.array[3][3].toPlainText()) * np.power(10., 9.))
                self.params.add('c_44', value=float(self.array[4][4].toPlainText()) * np.power(10., 9.))

            c_11 = float(self.array[1][1].toPlainText())
            c_12 = float(self.array[1][2].toPlainText())
            c_13 = float(self.array[1][3].toPlainText())
            c_33 = float(self.array[3][3].toPlainText())
            c_44 = float(self.array[4][4].toPlainText())
            c_66 = 0.5 * (c_11 - c_12)
            self.array[2][2].setText(str(c_11))
            self.array[2][3].setText(str(c_13))
            self.array[5][5].setText(str(c_44))
            self.array[6][6].setText(str(c_66))
        self.exit()

    def set_preprepare_array(self):
        if self.sym == "isotrop":
            for i in xrange(1, 7):
                for j in xrange(1, 7):
                    if i == j == 1 or (i == 1 and j == 2):
                        self.array[i][j].setReadOnly(False)
                    else:
                        self.array[i][j].setTextBackgroundColor(QColor(255, 200, 0))
                        self.array[i][j].setText(str(0))

        elif self.sym == "m-3m":
            for i in xrange(1, 7):
                for j in xrange(1, 7):
                    if i == j == 1 or (i == 1 and j == 2) or i == j == 4:
                        self.array[i][j].setReadOnly(False)
                    else:
                        self.array[i][j].setTextBackgroundColor(QColor(255, 200, 0))
                        self.array[i][j].setText(str(0))

        elif self.sym == "hexagonal":
            for i in xrange(1, 7):
                for j in xrange(1, 7):
                    if i == j == 1 or (i == 1 and j == 2) or i == j == 4 or i == j == 3 or (i == 1 and j == 3):
                        self.array[i][j].setReadOnly(False)
                    else:
                        self.array[i][j].setTextBackgroundColor(QColor(255, 200, 0))
                        self.array[i][j].setText(str(0))
        self.update()

    @staticmethod
    def create_text_output(widgate, text, posx, posy, sizex, sizey):
        text_obj = QTextEdit(widgate)
        text_obj.setReadOnly(True)
        text_obj.setTextBackgroundColor(QColor(255, 255, 255))
        text_obj.resize(sizex, sizey)
        text_obj.setLineWidth(QTextEdit.NoWrap)
        text_obj.move(posx, posy)
        text_obj.setText(text)
        return text_obj

    def send_and_close(self):
        self.exit()
        self.close()

    def exit(self):
        self.emit(SIGNAL(self.name), self.params)


class gui(QMainWindow):
    def __init__(self, name, params_matrix, params_inclusion="hallo"):
        QMainWindow.__init__(self)
        QMainWindow.resize(self, 700, 700)
        QMainWindow.setWindowTitle(self, name)
        self.params_matrix = params_matrix
        self.params_inclusion = params_inclusion
        self.widget1 = QWidget()
        # self.matrix_sym = ""
        # self.inclusion_sym = ""
        # return filename
        # Create an PyQT4 application object.
        # a = QApplication(sys.argv)

        # The QWidget widget is the base class of all user interface objects in PyQt4.
        # self.widgate = QMainWindow()
        # self.widgate1 = QWidget()
        # self.string_class = tacke_string(self.textbox.text_matrix(), self.widgate1)
        # Set window size.
        # self.widgate.resize(320, 240)
        # self.s = ""

        # Set window title
        # self.widgate.setWindowTitle("Hello World!")

        # Create main menu
        self.mainMenu = self.menuBar()
        self.mainMenu.setNativeMenuBar(True)
        self.fileMenu = self.mainMenu.addMenu('&amp;File')
        # self.toolbar = NavigationToolbar(self.canvas, self) #


        # Add exit button
        self.exitButton = QAction(QIcon('exit24.png'), 'Exit', self)
        self.exitButton.setShortcut('Ctrl+Q')
        self.exitButton.setStatusTip('Exit application')
        self.exitButton.triggered.connect(self.close)
        self.fileMenu.addAction(self.exitButton)
        #
        # # Create textbox
        # self.textbox = QLineEdit(self.widgate)
        # self.textbox.move(20, 20)
        # self.textbox.resize(280, 40)

        # Create combobox
        self.combo_matrix = self.create_comb_box(widgate=self, posx=180, posy=20)
        self.text_matrix = self.create_text_output(widgate=self, text="insert symmetry of the matrix",
                                                   posx=20, posy=20,
                                                   sizex=140, sizey=40)
        self.combo_matrix.activated[str].connect(self.set_matrix_text)

        # self.text = self.create_text_output(self, self.combo_matrix.currentText(), 20, 100, 50, 30)
        # Create a button in the window
        self.insert_csym_matrix = QPushButton("insert sym", self)
        self.insert_csym_matrix.move(300, 20)
        self.matrix_sym = self.combo_matrix.currentText()
        self.connect_csym_matrix()

        self.odf_button = QPushButton("insert ODF", self)
        self.odf_button.move(20, 150)
        self.textbox_odf = self.create_text_output(self, "odf path", 140, 150, 300, 30)
        self.odf_button.clicked.connect(self.open_odf_file)

        # widget, name, params, sym

        # self.emit(self, SIGNAL, "hallo", params_matrix)
        # handle the inclusion
        if self.params_inclusion is not None:
            self.combo_inclusion = self.create_comb_box(widgate=self, posx=180, posy=60)
            self.text_inclusion = self.create_text_output(widgate=self, text="insert symmetry of the inclusion",
                                                          posx=20, posy=60,
                                                          sizex=140, sizey=40)
            self.insert_csym_inclusion = QPushButton("insert sym", self)
            self.insert_csym_inclusion.move(300, 60)
            self.inclusion_sym = self.combo_inclusion.currentText()
            self.connect_csym_inclusion()
            self.combo_inclusion.activated[str].connect(self.set_inclusion_text)

        self.widget1.connect(self, SIGNAL("Matrix"), self.set_matrix_params)

        # Send Data outside
        self.button = QPushButton("send Data", self)
        self.button.move(80, 100)
        self.text = self.create_text_output(self, "hallo", 200, 100, 100, 30)
        self.button.clicked.connect(self.set_text)
        self.button.clicked.connect(self.send_data)

    def open_odf_file(self):
        filename = QFileDialog.getOpenFileName(self, 'Open ODF File', '/')
        self.textbox_odf.setText(filename)
        self.emit(SIGNAL('ODF_matrix_path'), filename)

    def set_text(self):
        self.text.setText(str(self.params_matrix['c_11'].value))

    def set_matrix_text(self, text):
        self.matrix_sym = text
        self.connect_csym_matrix()
        print(self.matrix_sym)

    def set_matrix_params(self, params):
        print(params['c_11'].value)
        self.params_matrix = params

    def connect_csym_matrix(self):
        self.insert_csym_matrix.clicked.connect(functools.partial(self.insert_const,
                                                                  "Matrix", self.params_matrix,
                                                                  self.matrix_sym))

    def connect_csym_inclusion(self):
        self.insert_csym_inclusion.clicked.connect(functools.partial(self.insert_const,
                                                                     "Inclusion", self.params_inclusion,
                                                                     self.inclusion_sym))

    def set_inclusion_text(self):
        self.inclusion_sym = self.combo_inclusion.currentText()
        self.connect_csym_inclusion()

    @staticmethod
    def create_text_output(widgate, text, posx, posy, sizex, sizey):
        text_obj = QTextEdit(widgate)
        text_obj.setReadOnly(True)
        text_obj.resize(sizex, sizey)
        text_obj.setLineWidth(QTextEdit.NoWrap)
        text_obj.move(posx, posy)
        text_obj.setText(text)
        return text_obj

    def create_Button(self):
        pass

    @staticmethod
    def create_comb_box(widgate, posx, posy):
        combo = QComboBox(widgate)
        combo.addItem("isotrop")
        combo.addItem("cubic")
        combo.addItem("hexagonal")
        combo.addItem("tetragonal_1 (not implementet jet)")
        combo.addItem("tetragonal_2 (not implementet jet)")
        combo.move(posx, posy)
        return combo

    def insert_const(self, *args):
        # name = kwargs["name"]
        # sym = kwargs["sym"]
        self.widget1 = insert_const_widget(*args)
        self.update()

    def send_data(self):
        self.emit(SIGNAL("params_matrix"), self.params_matrix)
        if self.params_inclusion is not None:
            self.emit(SIGNAL("params"), self.params_matrix, self.params_inclusion)


class Main(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.create_main_enviroment()
        # self.add_buton()

    def create_main_frame(self):
        self.main_frame = QWidget()

        self.fig = Figure((5.0, 4.0), dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        self.canvas.setFocusPolicy(Qt.StrongFocus)
        self.canvas.setFocus()

        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

        self.canvas.mpl_connect('key_press_event', self.on_key_press)

        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)  # the matplotlib canvas
        vbox.addWidget(self.mpl_toolbar)
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def create_main_enviroment(self):
        """
        creates the main enviroment (just the layout1)
        :return:
        """
        # self.layout=QGridLayout()
        self.resize(900, 900)
        self.centralWidget = CentralWidget(self)  # CentralWidget(self)
        self.setCentralWidget(self.centralWidget)

        # self.toolbar = QToolBar(self)
        # self.addToolBar(self.toolbar)

        # self.setLayout(self.layout)
        self.setWindowTitle("Fitting elastic constants")

    def add_buton(self):
        # self.button1 = QPushButton("hallo", self.centralWidget)
        # self.centralWidget.layout.addWidget(self.button1)
        pass

    def get_data2(self):
        return np.arange(20).reshape([4, 5]).copy()

    def on_draw(self):
        self.fig.clear()
        self.axes = self.fig.add_subplot(111)
        # self.axes.plot(self.x, self.y, 'ro')
        self.axes.imshow(self.data, interpolation='nearest')
        # self.axes.plot([1,2,3])
        self.canvas.draw()

    def on_key_press(self, event):
        key_press_handler(event, self.canvas, self.mpl_toolbar)


class CentralWidget(QWidget):
    def __init__(self, parent):
        super(CentralWidget, self).__init__(parent)
        # object containing the data (scattering data end texture data)
        self.data_object = handle_data.Data(sample_diameter=None)
        self.name_of_phase_dic = {1: "alpha", 2: "betha"}
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")
        self.region = phase_region_class()
        self.checkBoxPlotFits = QCheckBox("Plot fit's ")
        self.material = QLineEdit("iron")
        self.modi = self.create_modi_comb_box()
        self.output_filename = QLineEdit("Result_" + str(self.material.text()) + "_" + str(self.modi.currentText()))

        try:
            self.phase_peak_region = self.region.load()
            self.material.setText(self.region.material)
            self.change_outputfile_name()
            for i in xrange(len(self.phase_peak_region)):
                phasenr, phase, dat = self.phase_peak_region[i]
                self.name_of_phase_dic[phasenr] = phase
            self.loaded_peak_region = True
        except IOError:
            self.loaded_peak_region = False
            self.phase_peak_region = []

        # Load the data and sample specifications
        self.open_data_folders = QPushButton("open data")

        self.choose_experiment_comb_box = self.create_choose_exp_combbox()
        self.number_of_phases_selecttion = self.create_n_o_phases_combbox()
        self.number_of_datasets_under_strain = self.create_number_of_datasets_combbox()
        self.diameter = QLineEdit("6")
        # self.automate = QLabel("select hkl manualy?")  # QCheckBox("select hkl manually")
        self.select_hkl_setting_manualy_coice = QCheckBox("select hkl manually")  # self.create_jn_combbox()
        self.load_data_button = QPushButton('load Data')
        self.load_data_button.setEnabled(False)

        self.open_data_folders.clicked.connect(self.set_data_path_func)

        # select the experimental setup (SPODI, POLDI, ....)
        # self.load_data_button.clicked.connect(self.read_scattering_data_SPODI_case)
        self.choose_experiment_comb_box.currentIndexChanged.connect(self.connect_read_scattering_data)

        # add the central plot to display the data
        self.central_plot = matplotlibwidget.Preview("Select_Data")
        self.connect(self.central_plot, SIGNAL("phase_peak_region"), self.set_hkl_setting_and_fit_the_peaks_SPODI)
        self.connect(self.central_plot, SIGNAL("MIN_MAX_COLOR"), self.color_peaks)

        # handel the fitting process
        self.Fit_phase = QLabel("Fit phases: ")
        self.fit_phase_combbox = self.create_n_o_phases_combbox()
        self.name_of_phase = QLineEdit(self.name_of_phase_dic[int(str(self.fit_phase_combbox.currentText()))])
        self.name_of_phase.returnPressed.connect(self.change_name_of_phase)
        self.fit_phase_combbox.currentIndexChanged.connect(self.change_name_of_phase_qlineedit)
        self.modi_text = QLabel("Theory")

        self.ODF_text = QLabel("Texture?")
        self.text_jn = self.create_jn_combbox()

        self.do_the_fit_button = QPushButton('fitting Data')
        self.do_the_fit_button.setEnabled(False)
        self.do_the_fit_button.clicked.connect(self.fit_the_data)

        self.do_the_fit_gh_button = QPushButton('fit gh')
        self.do_the_fit_gh_button.setEnabled(False)
        self.do_the_fit_gh_button.clicked.connect(self.do_the_fit_gh)

        self.do_the_fit_tt_button = QPushButton('fit TT')
        self.do_the_fit_tt_button.setEnabled(False)
        self.do_the_fit_tt_button.clicked.connect(self.do_the_fit_tensile_test)

        self.material.returnPressed.connect(self.change_outputfile_name)
        self.modi.currentIndexChanged.connect(self.change_outputfile_name)
        # self.connect(self, SIGNAL("data"), self.central_plot.add_data)

        self.insert_startvals_button = QPushButton("insert params")
        self.insert_startvals_button.setEnabled(False)
        self.insert_startvals_button.clicked.connect(self.set_const_vals)

        self.plot_polefig_button = QPushButton("plot pole figure")
        self.plot_polefig_button.setEnabled(False)
        self.pfmaxval = QLineEdit('None')
        self.pfminval = QLineEdit('None')
        self.plot_data_button = QPushButton("plot data")
        self.plot_data_button.setEnabled(False)
        self.miller_h = QLineEdit("h")
        self.miller_k = QLineEdit("k")
        self.miller_l = QLineEdit("l")
        self.with_fit_combbox = self.create_jn_combbox()
        self.plot_polefig_button.clicked.connect(self.plot_pole_figure)
        self.plot_data_button.clicked.connect(self.plot_data_fuc)
        self.connect(self, SIGNAL('polefigurevals'), self.show_pole_figure)
        self.connect(self, SIGNAL('result_of_fit'), self.show_result_of_fit)
        self.layout_handling()

    def plot_data_fuc(self):
        h = int(str(self.miller_h.text()))
        k = int(str(self.miller_k.text()))
        l = int(str(self.miller_l.text()))
        method = str(self.modi.currentText())
        with_fit_t = str(self.with_fit_combbox.currentText())
        if with_fit_t == "Yes":
            with_fit = True
        else:
            with_fit = False
        if str(self.text_jn.currentText()) == "Yes":
            texture = True
        else:
            texture = False

        phase = int(str(self.fit_phase_combbox.currentText()))
        self.fit_object.plot_data(h, k, l, phase=phase, with_fit=with_fit, method=method, texture=texture)

    def set_params_phase_1(self, params):
        self.fit_object.set_params_phase_1(params)

    def set_params_phase_2(self, params):
        self.fit_object.set_params_phase_2(params)

    def plot_pole_figure_thread(self):
        self.do_the_fit_button.setEnabled(False)
        self.do_the_fit_gh_button.setEnabled(False)
        self.do_the_fit_tt_button.setEnabled(False)
        h = int(str(self.miller_h.text()))
        k = int(str(self.miller_k.text()))
        l = int(str(self.miller_l.text()))
        if str(self.fit_phase_combbox.currentText()) == "1":
            theta, r, VAL = self.data_object.odf_phase_1.plot_polfigure(h, k, l)
        if str(self.fit_phase_combbox.currentText()) == "2":
            theta, r, VAL = self.data_object.odf_phase_2.plot_polfigure(h, k, l)
        self.emit(SIGNAL('polefigurevals'), (h, k, l, theta, r, VAL))

        # fig, axs = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
        # p1 = axs.contourf(theta, r, VAL, 100)
        #
        # cbar = plt.colorbar(p1, ax=axs)
        # axs.set_title("pole figure {}{}{}\n".format(h, k, l))
        #
        # plt.show()
        thread.exit()

    def show_pole_figure(self, *args):
        # print(args)
        h, k, l, theta, r, VAL = args[0]
        self.do_the_fit_button.setEnabled(True)
        self.do_the_fit_gh_button.setEnabled(True)
        self.do_the_fit_tt_button.setEnabled(True)


        # print(r)
        if self.pfmaxval.text() == 'None' or self.pfminval.text() == 'None':
            v = np.linspace(VAL.min(), VAL.max(), 100, endpoint=True)
        else:
            v = np.linspace(float(str(self.pfminval.text())), float(str(self.pfmaxval.text())), 100, endpoint=True)

        # set ticks for the colorbar
        ticks = []
        dmaxmin = max(v) - min(v)
        number = 10
        step = float("{:.1f}".format(dmaxmin / 10))
        i = 0
        ii = []
        while min(v) < 1 - i * step:
            ii.append(-i)
            i += 1
        lll = []
        for m in ii:
            lll.insert(0, m)
        i = 1
        while len(lll) <= number:
            lll.append(i)
            i += 1
        ticks = 1+np.array(lll)*step

        # set ticks for psi
        ticklables = np.rad2deg(r[0])

        print('min: {}, max: {}'.format(VAL.min(), VAL.max()))

        def mapr(r):
            """
            remap the radial axis
            :param r: radial value
            """
            return np.rad2deg(np.arctan(r)*2)

        # make the plot
        fig, axs = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))

        axs.grid(True)
        print("yticks: ", axs.get_yticks())
        ticklables = []
        for i in mapr(axs.get_yticks()):
            ticklables.append(float('{:.1f}'.format(i)))
        axs.yaxis.set_ticklabels(ticklables)
        # p1.ax.set_tichlabels(ticklables)
        # axs.set_yticks(ticklabels=ticklables)
        axs.set_title("pole figure {}{}{}\n".format(h, k, l))

        # draw the contourplot
        p1 = axs.contourf(theta, r, VAL, v)  # 100,,  vmin=0.7, vmax=1.8
        cbar = plt.colorbar(p1, ax=axs, ticks=ticks)  # , ticks=v)  # norm=mpl.colors.Normalize(vmin=0.7, vmax=1.8))
        # cbar.set_clim(0.7, 1.8)

        # axs.set_theta_zero_location("S")
        # axs.set_theta_offset(pi)

        plt.show()

    def plot_pole_figure(self):
        thread.start_new_thread(self.plot_pole_figure_thread, ())

    def change_name_of_phase_qlineedit(self):
        self.name_of_phase.setText(self.name_of_phase_dic[int(str(self.fit_phase_combbox.currentText()))])

    def change_name_of_phase(self):
        self.name_of_phase_dic[int(str(self.fit_phase_combbox.currentText()))] = str(self.name_of_phase.text())
        self.phase_peak_region[int(str(self.fit_phase_combbox.currentText())) - 1][1] = str(self.name_of_phase.text())

    def change_outputfile_name(self):
        text = "Result_" + str(self.material.text()) + "_" + str(self.modi.currentText())
        try:
            self.fit_object.material = str(self.material)
            self.fit_object_gh.material = str(self.material)
            self.fit_object_TT.material = str(self.material)
        except AttributeError:
            pass
        self.output_filename.setText(text)

    def connect_read_scattering_data(self):
        """
        select the setup (SPODI, POLDI, ...)
        :return:
        """
        if self.choose_experiment_comb_box.currentText() == "SPODI":
            self.load_data_button.clicked.connect(self.read_scattering_data_SPODI_case)
        elif self.choose_experiment_comb_box.currentText() == "POLDI":
            print("need to be implemented")

    def set_data_path_func(self):
        """
        open a widget to select the directory's of the input data
        use different methods vor the different experiments
        :return:
        """
        n_o_p = self.number_of_phases_selecttion.currentText()
        n_o_p = int(n_o_p)
        n_o_d = int(self.number_of_datasets_under_strain.currentText())
        print("phases", n_o_p)
        if self.choose_experiment_comb_box.currentText() == "SPODI":
            self.widget_set_data_path = LOAD_SPODI_DATA("set data path", number_of_phases=n_o_p,
                                                        number_of_straind_datasets=n_o_d)
            self.connect(self.widget_set_data_path, SIGNAL("data_dir_list"), self.receve_the_pathes_SPODI_case)
            try:
                self.phase_peak_region = self.region.load()
                self.material.setText(self.region.material)
                self.loaded_peak_region = True
            except IOError:
                self.loaded_peak_region = False
                self.phase_peak_region = []

        elif self.choose_experiment_comb_box.currentText() == "POLDI":
            self.widget_set_data_path = LOAD_STANDARD_DATA("set data path", number_of_phases=n_o_p)
            self.connect(self.widget_set_data_path, SIGNAL("data_dir_list"),
                         self.recive_the_pathes_of_standard_data_format)

    def recive_the_pathes_of_standard_data_format(self, *args):
        print(args)
        odf1, odf2, data_file, phase_key_dict, material = args[0]
        odf1 = str(odf1)
        odf2 = str(odf2)
        data_file = str(data_file)
        self.name_of_phase_dic[1] = phase_key_dict[1]
        material = str(material).strip()
        self.material.setText(material)
        self.change_outputfile_name()
        self.change_name_of_phase_qlineedit()
        try:
            self.name_of_phase_dic[2] = phase_key_dict[2]
        except KeyError:
            pass

        if odf2 == "None":
            odf2 = None
        print(odf1)
        print(odf2)
        self.path_of_odf_phase1 = odf1
        self.path_of_odf_phase2 = odf2
        self.path_of_standad_datafile = data_file
        self.load_data_button.setEnabled(True)
        print(self.name_of_phase_dic)
        self.load_data_button.clicked.connect(self.read_scattering_data_standard_case)

    def read_scattering_data_standard_case(self):
        print("#########################################\n",
              "odf1_path:", self.path_of_odf_phase1,
              "\nodf2_path:", self.path_of_odf_phase2,
              "\n############################################")
        self.data_object = handle_data.AllData(odf_phase_1_file=self.path_of_odf_phase1,
                                               odf_phase_2_file=self.path_of_odf_phase2)

        self.data_object.read_data(self.path_of_standad_datafile, phase_name_dict=self.name_of_phase_dic)

        self.do_the_fit_button.setEnabled(True)
        self.do_the_fit_gh_button.setEnabled(True)
        self.do_the_fit_tt_button.setEnabled(True)
        self.insert_startvals_button.setEnabled(True)
        self.plot_polefig_button.setEnabled(True)
        self.plot_data_button.setEnabled(True)

        self.fit_object = Modells.FitStrainWithTexture(data_object=self.data_object, material=self.material.text())
        self.fit_object_gh = Modells.FitGneupelHerold(data_object=self.data_object, material=self.material.text())
        self.fit_object_TT = Modells.TensileTest(data_object=self.data_object, material=self.material.text())
        self.fit_object.print_params()

    def receve_the_pathes_SPODI_case(self, *args):
        """
        receve the data path
        :param args:
        :return:
        """
        print(args, len(args[0]))
        odf1, odf2, Data = args[0]
        odf1 = str(odf1)
        odf2 = str(odf2)
        Data = str(Data)
        # for i in data_dir_list:
        #     i = str(i)
        #     print(i, type(i))
        print(odf1)
        if odf2 == "None":
            odf2 = None
        print(odf1, odf2, Data)
        self.path_of_odf_phase1 = odf1
        self.path_of_odf_phase2 = odf2
        self.path_of_unstraind_data = Data
        # self.path_of_data_under_strain = data_dir_list
        self.data_object = handle_data.SPODIData(sample_diameter=int(str(self.diameter.text())),
                                                 odf_phase_1_file=self.path_of_odf_phase1,
                                                 odf_phase_2_file=self.path_of_odf_phase2)
        self.data_object.load_data(self.path_of_unstraind_data)
        self.select_hkl_SPODI_Data()  # plot the data
        if self.loaded_peak_region:
            self.color_peakregion()
        self.load_data_button.setEnabled(True)
        self.load_data_button.clicked.connect(self.read_scattering_data_SPODI_case)

    def read_scattering_data_SPODI_case(self):
        # Bool = self.select_hkl_setting_manualy_coice.isChecked()
        # if self.select_hkl_setting_manualy_coice.isChecked()==True:  # .currentText() == "Yes":
        #     Bool = True



        # self.Data_Iron = methods.Data_old(str(self.odf_phase_1_path.text()), 6)
        # self.Data_Iron.read_scattering_SPODI_data(path_of_unstraind_data=str(self.path_of_unstraind_data.text()),
        #                                           path_of_straind_data=str(self.path_of_straind_data_1.text()))
        self.select_hkl_SPODI_Data()
        self.central_plot.roi_Button.setEnabled(False)
        if self.select_hkl_setting_manualy_coice.isChecked() or not self.loaded_peak_region:
            self.central_plot.roi_Button.setEnabled(True)
        else:
            self.data_object.fit_all_data(peak_regions_phase=self.phase_peak_region,
                                          plot=self.checkBoxPlotFits.isChecked(),
                                          material=str(self.material.text()))
            self.do_the_fit_button.setEnabled(True)
            self.do_the_fit_gh_button.setEnabled(True)
            self.do_the_fit_tt_button.setEnabled(True)
            self.insert_startvals_button.setEnabled(True)
            self.plot_polefig_button.setEnabled(True)
            self.plot_data_button.setEnabled(True)
            # self.Data_Iron.fit_all_peaks()
        self.fit_object = Modells.FitStrainWithTexture(data_object=self.data_object, material=self.material.text())
        self.fit_object_gh = Modells.FitGneupelHerold(data_object=self.data_object, material=self.material.text())
        self.fit_object_TT = Modells.TensileTest(data_object=self.data_object, material=self.material.text())
        # self.connect(self, SIGNAL('1'), self.fit_object.set_params_phase_1)
        # self.connect(self, SIGNAL('2'), self.fit_object.set_params_phase_2)
        self.fit_object.print_params()

    def select_hkl_SPODI_Data(self):
        x_data, y_data = self.data_object.get_sum_data()
        # print("data_x:", x_data)
        self.central_plot.add_xy_data(x_data, y_data)

    def color_peaks(self, MIN, MAX, color):
        x_data, y_data = self.data_object.get_sum_data()
        x_data = x_data[MIN:MAX]
        y_data = y_data[MIN:MAX]
        for i in xrange(len(x_data)):
            y_data[i] = max(y_data)
        y_data[0] = y_data[-1] = 0
        # print("data_x:", x_data)
        self.central_plot.color_xy_data(x_data, y_data, color=color)

    def color_peakregion(self):
        color = matplotlibwidget.ColorGenerator()
        for phasenr, phase, hkl_region in self.phase_peak_region:
            for h, k, l, MIN, MAX, double, peak in hkl_region:
                self.color_peaks(MIN, MAX, color.get_color())

    def set_hkl_setting_and_fit_the_peaks_SPODI(self, value):
        self.phase_peak_region = value
        print("------------------------------------")
        print("peak region: ", self.phase_peak_region)

        self.region.save(self.phase_peak_region)
        # self.region.load()
        # np.save(".\\phase_peak_region", np.array(self.phase_peak_region))
        print("saved peak region")
        print("------------------------------------")
        self.select_hkl_SPODI_Data()
        self.color_peakregion()
        self.data_object.fit_all_data(peak_regions_phase=self.phase_peak_region)
        # self.Data_Iron.set_hkl_setting(self.phase_peak_region)
        # self.Data_Iron.fit_all_peaks()

        self.do_the_fit_button.setEnabled(True)
        self.do_the_fit_gh_button.setEnabled(True)
        self.do_the_fit_tt_button.setEnabled(True)
        self.insert_startvals_button.setEnabled(True)
        self.plot_polefig_button.setEnabled(True)
        self.plot_data_button.setEnabled(True)
        # print("coming from other class:", "\n", self.phase_peak_region)

    def do_the_fit(self):
        self.plot_polefig_button.setEnabled(False)
        Bool = False
        if self.text_jn.currentText() == "Yes":
            Bool = True
        print(self.modi.currentText(), "\n",
              Bool)
        print("----------------------------\n",
              "Fit the data using model: " +
              self.modi.currentText() +
              "\n",
              "With texture: ", self.text_jn.currentText(), "\n",
              "Fitting phase: ", self.fit_phase_combbox, "\n",
              "----------------------------")

        # self.fit_object = Modells.FitStrainWithTexture(data_object=self.data_object)

        result = self.fit_object.do_the_fitting(filename=str(self.output_filename.text()),
                                                material=self.material.text(),
                                                method=str(self.modi.currentText()),
                                                phase=int(str(self.fit_phase_combbox.currentText())),
                                                phase_name=self.name_of_phase_dic[
                                                    int(str(self.fit_phase_combbox.currentText()))],
                                                texture=Bool)
        text = "Finnished calculation\nresults are stored under {}".format(result[1])
        self.plot_polefig_button.setEnabled(True)
        self.emit(SIGNAL('result_of_fit'), (text, result))
        thread.exit()

    def do_the_fit_gh(self):
        # self.plot_polefig_button.setEnabled(False)
        Bool = False
        if self.text_jn.currentText() == "Yes":
            Bool = True
        print(self.modi.currentText(), "\n",
              Bool)
        print("----------------------------\n",
              "Fit the data using model: " +
              self.modi.currentText() +
              "\n",
              "With texture: ", self.text_jn.currentText(), "\n",
              "Fitting phase: ", self.fit_phase_combbox, "\n",
              "----------------------------")

        # self.fit_object = Modells.FitStrainWithTexture(data_object=self.data_object)

        result = self.fit_object_gh.do_the_fitting_gneupel_herold(filename=str(self.output_filename.text()),
                                                                  material=self.material.text(),
                                                                  method=str(self.modi.currentText()),
                                                                  phase=int(str(self.fit_phase_combbox.currentText())),
                                                                  phase_name=self.name_of_phase_dic[
                                                                      int(str(self.fit_phase_combbox.currentText()))],
                                                                  texture=Bool)
        text = "Finnished calculation\nresults are stored under {}".format(result[1])
        plot_dic = result[2]
        self.show_result_of_fit(text, result, plot_dic)
        # self.plot_polefig_button.setEnabled(True)

    def cos2psi_plot(self, plots_dic):
        for figname, data in plots_dic.iteritems():
            xdata, ydata, yerr = data[0]
            Psi, val, s1, s1err, s2, s2err = data[1]
            plt.figure(figname)
            plt.errorbar(xdata, ydata, yerr=yerr, fmt='bo', label="Data")
            plots_dic[figname].append([xdata, ydata, yerr])

            plt.plot(Psi, val, 'r-',
                     label="s1 = {:.3g} $\pm$ {:.3g}\ns2 = {:.3g} $\pm$ {:.3g}".format(s1, s1err, s2, s2err))

            plt.xlabel('$\cos^2(\Psi)$')
            plt.ylabel('$\epsilon/\sigma$')
            plt.legend(loc='upper left')
            plt.xlim([0, 1])
            print("savefig, ", figname, ".svg")
            filename = ".\\sin2psi-plots\\" + figname + ".svg"
            if not os.path.exists(os.path.dirname(filename)):
                os.makedirs(os.path.dirname(filename))
            plt.savefig(".\\sin2psi-plots\\" + figname + ".svg", format="svg")
            plt.savefig(".\\sin2psi-plots\\" + figname + ".pdf", format="pdf")
            plt.savefig(".\\sin2psi-plots\\" + figname + ".png", format="png")
        plt.show()

    def show_result_of_fit(self, *args):
        try:
            try:
                text, result = args[0]
            except ValueError:
                print(len(args))
                text, result = args

        except ValueError:
            try:
                text, result, plots_dic = args[0]
            except ValueError:
                print(len(args))
                text, result, plots_dic = args
            print(text)
            self.cos2psi_plot(plots_dic=plots_dic)
        mbox = QMessageBox()
        mbox.standardButtons()
        mbox.setIcon(QMessageBox.Information)
        mbox.setText(text)
        mbox.setDetailedText(lm.fit_report(result[0].params))
        mbox.exec_()

    def fit_the_data(self):
        thread.start_new_thread(self.do_the_fit, ())

    def do_the_fit_tensile_test(self):
        # self.plot_polefig_button.setEnabled(False)
        Bool = False
        if self.text_jn.currentText() == "Yes":
            Bool = True
        print(self.modi.currentText(), "\n",
              Bool)
        print("----------------------------\n",
              "Fit the data using model: " +
              self.modi.currentText() +
              "\n",
              "With texture: ", self.text_jn.currentText(), "\n",
              "Fitting phase: ", self.fit_phase_combbox, "\n",
              "----------------------------")

        # self.fit_object = Modells.FitStrainWithTexture(data_object=self.data_object)

        result = self.fit_object_TT.do_the_fitting_gneupel_herold(filename=str(self.output_filename.text()),
                                                                  material=self.material.text(),
                                                                  method=str(self.modi.currentText()),
                                                                  phase=int(str(self.fit_phase_combbox.currentText())),
                                                                  phase_name=self.name_of_phase_dic[
                                                                      int(str(self.fit_phase_combbox.currentText()))],
                                                                  texture=Bool)
        text = "Finnished calculation\nresults are stored under {}".format(result[1])
        plot_dic = result[2]
        self.show_result_of_fit(text, result, plot_dic)
        # self.plot_polefig_button.setEnabled(True)

    def layout_handling(self):
        # Layout handling
        layout = QVBoxLayout()
        layout_ok_and_cancel_button = QHBoxLayout()
        layout_fitting = QHBoxLayout()
        layout_fitting_2 = QHBoxLayout()

        # layout.addWidget(self.toolbar)
        # layout.addWidget(self.canvas)
        layout.addWidget(self.central_plot)
        layout_ok_and_cancel_button.addStretch(1)
        layout_ok_and_cancel_button.addWidget(self.ok_button)
        layout_ok_and_cancel_button.addWidget(self.cancel_button)

        # Load the data
        HLine = []
        for i in xrange(10):
            HLine.append(self.HLine())

        layout_load_data_h1 = QHBoxLayout()
        layout_load_data_h2 = QHBoxLayout()
        layout_load_data_h3 = QHBoxLayout()
        layout_load_data_h1.addWidget(self.label("Load the Data: "))
        layout_load_data_h2.addWidget(self.label("chose expermient: "))
        layout_load_data_h2.addWidget(self.choose_experiment_comb_box)
        layout_load_data_h2.addWidget(self.label("# of phases: "))
        layout_load_data_h2.addWidget(self.number_of_phases_selecttion)
        # layout_load_data_h2.addWidget(self.label("# of dataset's under strain: "))
        # layout_load_data_h2.addWidget(self.number_of_datasets_under_strain)
        layout_load_data_h2.addWidget(self.open_data_folders)
        layout_load_data_h2.addWidget(self.label("set diameter in mm: "))
        layout_load_data_h2.addWidget(self.diameter)
        # layout_load_data_h3.addWidget(self.automate)
        layout_load_data_h2.addWidget(self.select_hkl_setting_manualy_coice)
        layout_load_data_h2.addWidget(self.checkBoxPlotFits)
        layout_load_data_h2.addWidget(self.load_data_button)
        layout.addLayout(layout_load_data_h1)
        layout.addLayout(layout_load_data_h2)
        # layout.addLayout(layout_load_data_h3)
        layout.addWidget(HLine[2])

        # handel the fitting process
        layout.addWidget(self.label("Fit the data: "))
        layout_fitting.addWidget(self.Fit_phase)
        layout_fitting.addWidget(self.fit_phase_combbox)
        layout_fitting.addWidget(self.label("phase name:"))
        layout_fitting.addWidget(self.name_of_phase)
        # layout_fitting.addWidget(self.insert_startvals_button)
        layout_fitting.addWidget(self.modi_text)
        layout_fitting.addWidget(self.modi)
        layout_fitting.addWidget(self.ODF_text)
        layout_fitting.addWidget(self.text_jn)
        layout_fitting.addWidget(self.label("material:"))
        layout_fitting.addWidget(self.material)
        layout_fitting.addWidget(self.label("output filename:"))
        layout_fitting.addWidget(self.output_filename)
        layout_fitting.addWidget(self.do_the_fit_button)

        layout_fitting_2.addWidget(self.label('h:'))
        layout_fitting_2.addWidget(self.miller_h)
        layout_fitting_2.addWidget(self.label('k:'))
        layout_fitting_2.addWidget(self.miller_k)
        layout_fitting_2.addWidget(self.label('l:'))
        layout_fitting_2.addWidget(self.miller_l)
        layout_fitting_2.addWidget(self.label('pfminval:'))
        layout_fitting_2.addWidget(self.pfminval)
        layout_fitting_2.addWidget(self.label('pfmsxval:'))
        layout_fitting_2.addWidget(self.pfmaxval)
        layout_fitting_2.addWidget(self.plot_polefig_button)
        layout_fitting_2.addWidget(self.label('with_fit:'))
        layout_fitting_2.addWidget(self.with_fit_combbox)
        layout_fitting_2.addWidget(self.plot_data_button)
        layout_fitting_2.addWidget(self.do_the_fit_gh_button)
        layout_fitting_2.addWidget(self.do_the_fit_tt_button)
        layout_fitting_2.addWidget(self.insert_startvals_button)

        layout.addLayout(layout_fitting)
        layout.addLayout(layout_fitting_2)
        layout.addLayout(layout_ok_and_cancel_button)

        self.setLayout(layout)

    def label(self, text):
        label = QLabel(text)
        return label

    @staticmethod
    def create_modi_comb_box():
        combo = QComboBox()
        combo.addItem("reus")
        combo.addItem("voigt")
        combo.addItem("hill")
        combo.addItem("eshelby")
        return combo

    @staticmethod
    def create_jn_combbox():
        combo = QComboBox()
        combo.addItem("Yes")
        combo.addItem("No")
        return combo

    @staticmethod
    def create_number_of_datasets_combbox():
        combo = QComboBox()
        for i in range(1, 11):
            combo.addItem("{0}".format(i))
        return combo

    @staticmethod
    def create_choose_exp_combbox():
        combo = QComboBox()
        combo.addItem("SPODI")
        combo.addItem("POLDI")
        return combo

    @staticmethod
    def create_n_o_phases_combbox():
        combo = QComboBox()
        combo.addItem("1")
        combo.addItem("2")
        return combo

        # def select_ODF_phase_1_func(self):
        #     filename = QFileDialog.getOpenFileName(self, 'Open ODF File', '/')  # (self, 'Open ODF File', '/')
        #     filename = os.path.normpath(str(filename))
        #     self.odf_phase_1_path.setText(filename)
        #     print(self.odf_phase_1_path.text())
        #
        # def select_unstraind_func(self):
        #     filename = QFileDialog.getExistingDirectory(self, 'Open unstraind data', '/')  # (self, 'Open ODF File', '/')
        #     filename = os.path.normpath(str(filename))
        #     self.path_of_unstraind_data.setText(filename + "\\")
        #     print(self.path_of_unstraind_data.text())
        #
        # def select_straind_func_1(self):
        #     filename = QFileDialog.getExistingDirectory(self, 'Open straind data', '/')  # (self, 'Open ODF File', '/')
        #     filename = os.path.normpath(str(filename))
        #     self.path_of_straind_data.setText(filename + "\\")
        #     print(self.path_of_straind_data.text())

    def HLine(self):
        horizontal_line = QFrame()
        horizontal_line.setFrameStyle(QFrame.HLine)
        horizontal_line.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Expanding)
        return horizontal_line

    def set_const_vals(self):

        if int(str(self.fit_phase_combbox.currentText())) == 1:
            sym = self.data_object.odf_phase_1.crystal_symmetry
        else:
            sym = self.data_object.odf_phase_2.crystal_symmetry
        print(str(self.fit_phase_combbox.currentText()))
        self.insert_params = insert_const_widget()

        self.insert_params.show_of(name=str(self.fit_phase_combbox.currentText()),
                                   sym=sym)
        # self.connect(self, SIGNAL('1'), self.fit_object.set_params_phase_1)
        # self.connect(self, SIGNAL('2'), self.fit_object.set_params_phase_2)
        self.connect(self.insert_params, SIGNAL("1"), self.set_params_phase_1)
        self.connect(self.insert_params, SIGNAL("2"), self.set_params_phase_2)
        self.fit_object.print_params()

        # def generate


class HandleCVals(object):
    def __init__(self, number, center):
        self.vals = np.zeros((number, number, number))  # index of [c44, c12, c11]
        self.c440, self.c120, self.c110 = center
        self.lookupc11 = self.generate_lookup(number, self.c110)
        self.lookupc12 = self.generate_lookup(number, self.c120)
        self.lookupc44 = self.generate_lookup(number, self.c440)

    def generate_lookup(self, number, center):
        diff = abs(center * (1 - 0.15) - center * (1 + 0.15)) / number
        return np.arange(center * (1 - 0.15), center * (1 + 0.15), diff)

    def get_val(self, ic44, ic12, ic11):
        return self.lookupc44[ic44], self.lookupc44[ic12], self.lookupc44[ic11], self.vals[ic44, ic12, ic11]

    def iterator(self):
        ic44, ic12, ic11 = 0, 0, 0
        while ic44 < len(self.lookupc44):
            while ic12 < len(self.lookupc12):
                while ic11 < len(self.lookupc11):
                    yield self.lookupc44[ic44], self.lookupc44[ic12], self.lookupc44[ic11], self.vals[ic44, ic12, ic11]


class LOAD_SPODI_DATA(QWidget):
    def __init__(self, name, number_of_phases=1, number_of_straind_datasets=1):
        super(LOAD_SPODI_DATA, self).__init__()
        self.setWindowTitle(name)

        self.ok_button = QPushButton("OK")
        # self.cancel_button = QPushButton("Cancel")
        self.hkl_setting = []

        self.odf_phase_1_button = QPushButton("select ODF of phase 1")
        self.odf_phase_2_button = QPushButton("select ODF of phase 2")
        self.unstraind_button = QPushButton("select unstraind")

        self.odf_phase_1_path = QLineEdit(
            "H:\Masterarbeit STRESS-SPEC\Daten\Daten-bearbeitet\Daten\ODF_Daten\Duplex gezogen\DUBNA_pol\FE_ODF_compleat.txt")
        self.odf_phase_2_path = QLineEdit("None")  # "AL_textur_complet.txt"
        self.path_of_unstraind_data = QLineEdit(
            "H:\\Masterarbeit STRESS-SPEC\\Daten\\Daten-bearbeitet\\Daten\\Duplex_Stahl_gezogen\\200N\\")

        self.straind_data = []
        self.straind_data_lable = []
        self.straind_data_button = []

        # create boxes for all datafiles of different load
        # for i in xrange(number_of_straind_datasets):
        #     self.straind_data.append(QLineEdit("H:\\Masterarbeit STRESS-SPEC\\Daten\\Daten-bearbeitet\\Daten\\Duplex_Stahl_gezogen\\5KN\\"))
        #     self.straind_data_lable.append(QLabel("select straind data %i:" % (i)))
        #     self.straind_data_button.append(QPushButton("select straind %i:" % (i)))

        self.select_odf_phase_1 = QLabel("select odf phase 1:   ")
        self.select_odf_phase_2 = QLabel("select odf phase 2:   ")
        self.select_unstraind = QLabel("select unstraind data:")

        self.odf_phase_1_button.clicked.connect(self.select_odf_phase_1_func)
        self.odf_phase_2_button.clicked.connect(self.select_odf_phase_2_func)
        self.unstraind_button.clicked.connect(self.select_unstraind_func)

        # for i in xrange(len(self.straind_data_button)):
        #     self.straind_data_button[i].clicked.connect(partial(self.select_straind_func, i))
        # self.straind_button_1.clicked.connect(self.select_straind_func_1)

        if number_of_phases == 1:
            self.odf_phase_2_button.setDisabled(True)
            self.odf_phase_2_path.setReadOnly(True)
        if number_of_phases == 2:
            self.odf_phase_2_button.setDisabled(False)
            self.odf_phase_2_path.setReadOnly(False)

        self.ok_button.clicked.connect(self.emit_and_quit)

        self.layout_handling()
        self.show()

    def layout_handling(self):
        # Layout handling
        self.resize(800, 800)
        layout = QVBoxLayout()
        layout1 = QHBoxLayout()
        layout_odf_phase_1_input = QHBoxLayout()
        layout_odf_phase_2_input = QHBoxLayout()
        # layout_straind_data = []
        # for i in xrange(len(self.straind_data)):
        #     layout_straind_data.append(QHBoxLayout())
        # layout_straind_data_1 = QHBoxLayout()
        # layout_straind_data_2 = QHBoxLayout()
        # layout_straind_data_3 = QHBoxLayout()
        # layout_straind_data_4 = QHBoxLayout()
        layout_unstraind_data = QHBoxLayout()

        # layout.addWidget(self.toolbar)
        # layout.addWidget(self.canvas)

        layout1.addStretch(1)
        layout1.addWidget(self.ok_button)
        # layout1.addWidget(self.cancel_button)

        # insert odf phase 1 path
        layout_odf_phase_1_input.addWidget(self.select_odf_phase_1)
        layout_odf_phase_1_input.addWidget(self.odf_phase_1_path)
        layout_odf_phase_1_input.addWidget(self.odf_phase_1_button)

        # insert odf phase 2 path
        layout_odf_phase_2_input.addWidget(self.select_odf_phase_2)
        layout_odf_phase_2_input.addWidget(self.odf_phase_2_path)
        layout_odf_phase_2_input.addWidget(self.odf_phase_2_button)

        # insert straind datasets
        # for i in xrange(len(layout_straind_data)):
        #     layout_straind_data[i].addWidget(self.straind_data_lable[i])
        #     layout_straind_data[i].addWidget(self.straind_data[i])
        #     layout_straind_data[i].addWidget(self.straind_data_button[i])

        # insert unstraind data
        layout_unstraind_data.addWidget(self.select_unstraind)
        layout_unstraind_data.addWidget(self.path_of_unstraind_data)
        layout_unstraind_data.addWidget(self.unstraind_button)

        layout.addLayout(layout_odf_phase_1_input)
        layout.addLayout(layout_odf_phase_2_input)
        layout.addLayout(layout_unstraind_data)
        # for i in layout_straind_data:
        #     layout.addLayout(i)
        # layout.addLayout(layout_straind_data_1)
        # layout.addLayout(layout_straind_data_2)
        # layout.addLayout(layout_straind_data_3)
        # layout.addLayout(layout_straind_data_4)

        layout.addLayout(layout1)

        self.setLayout(layout)

    def select_odf_phase_1_func(self):
        path = os.path.split(str(self.odf_phase_1_path.text()))[0]
        filename = QFileDialog.getOpenFileName(self, 'Open ODF of phase 1 File', path)  # (self, 'Open ODF File', '/')
        filename = os.path.normpath(str(filename))
        self.odf_phase_1_path.setText(filename)
        print(self.odf_phase_1_path.text())

    def select_odf_phase_2_func(self):
        path = os.path.split(str(self.odf_phase_2_path.text()))[0]
        filename = QFileDialog.getOpenFileName(self, 'Open ODF of phase 2 File', path)  # (self, 'Open ODF File', '/')
        filename = os.path.normpath(str(filename))
        self.odf_phase_2_path.setText(filename)
        print(self.odf_phase_2_path.text())

    def select_unstraind_func(self):
        path = os.path.split(str(self.path_of_unstraind_data.text()))[0]
        filename = QFileDialog.getExistingDirectory(self, 'Open data', path)  # (self, 'Open ODF File', '/')
        filename = os.path.normpath(str(filename))
        self.path_of_unstraind_data.setText(filename + "\\")
        print(self.path_of_unstraind_data.text())

    def select_straind_func(self, i):
        path = os.path.split(str(self.straind_data[i].text()))[0]
        filename = QFileDialog.getExistingDirectory(self, 'Open straind data', path)  # (self, 'Open ODF File', '/')
        filename = os.path.normpath(str(filename))
        self.straind_data[i].setText(filename + "\\")
        print(self.straind_data[i].text())

    def emit_and_quit(self):
        """
        emit the selected path's
        :return:
        """
        odf1 = self.odf_phase_1_path.text()
        odf2 = self.odf_phase_2_path.text()
        Data = self.path_of_unstraind_data.text()
        data_dir_list = []
        data_dir_list.append(self.path_of_unstraind_data.text())
        # for i in xrange(len(self.straind_data)):
        #     data_dir_list.append(self.straind_data[i].text())
        self.emit(SIGNAL("data_dir_list"), (odf1, odf2, Data))
        self.close()


class LOAD_STANDARD_DATA(QWidget):
    def __init__(self, name, number_of_phases=1):
        super(LOAD_STANDARD_DATA, self).__init__()
        self.setWindowTitle(name)

        self.ok_button = QPushButton("OK")
        # self.cancel_button = QPushButton("Cancel")
        self.phase_key_dict = {}

        self.odf_phase_1_button = QPushButton("select ODF of phase 1")
        self.odf_phase_2_button = QPushButton("select ODF of phase 2")
        self.data_file_button = QPushButton("select data file")

        self.odf_phase_1_path = QLineEdit(
            "H:\\Masterarbeit STRESS-SPEC\\Daten\\Daten-bearbeitet\\Stahl ST37\\ST37_textur_complet_recalc.txt")
        self.odf_phase_2_path = QLineEdit("None")  # "AL_textur_complet.txt"
        self.path_of_data_file = QLineEdit("H:\Masterarbeit STRESS-SPEC\Daten\create_input_file\DUPLEX_gezogen.dat")

        self.select_odf_phase_1 = QLabel("select odf phase 1:   ")
        self.select_odf_phase_2 = QLabel("select odf phase 2:   ")
        self.select_data_file = QLabel("select data file:")

        self.odf_phase_1_button.clicked.connect(self.select_odf_phase_1_func)
        self.odf_phase_2_button.clicked.connect(self.select_odf_phase_2_func)
        self.data_file_button.clicked.connect(self.select_data_file_func)

        self.odf_phase_1_path.setReadOnly(True)
        self.odf_phase_1_button.setEnabled(False)

        # if number_of_phases == 1:
        self.odf_phase_2_button.setDisabled(True)
        self.odf_phase_2_path.setReadOnly(True)
        # if number_of_phases == 2:
        #     self.odf_phase_2_button.setDisabled(False)
        #     self.odf_phase_2_path.setReadOnly(False)

        self.ok_button.clicked.connect(self.emit_and_quit)

        self.layout_handling()
        self.show()

    def layout_handling(self):
        # Layout handling
        self.resize(800, 800)
        layout = QVBoxLayout()
        layout1 = QHBoxLayout()
        layout_odf_phase_1_input = QHBoxLayout()
        layout_odf_phase_2_input = QHBoxLayout()
        layout_data_file_data = QHBoxLayout()

        layout1.addStretch(1)
        layout1.addWidget(self.ok_button)
        # layout1.addWidget(self.cancel_button)

        # insert odf phase 1 path
        layout_odf_phase_1_input.addWidget(self.select_odf_phase_1)
        layout_odf_phase_1_input.addWidget(self.odf_phase_1_path)
        layout_odf_phase_1_input.addWidget(self.odf_phase_1_button)

        # insert odf phase 2 path
        layout_odf_phase_2_input.addWidget(self.select_odf_phase_2)
        layout_odf_phase_2_input.addWidget(self.odf_phase_2_path)
        layout_odf_phase_2_input.addWidget(self.odf_phase_2_button)

        # insert data file
        layout_data_file_data.addWidget(self.select_data_file)
        layout_data_file_data.addWidget(self.path_of_data_file)
        layout_data_file_data.addWidget(self.data_file_button)

        layout.addLayout(layout_data_file_data)
        layout.addLayout(layout_odf_phase_1_input)
        layout.addLayout(layout_odf_phase_2_input)

        layout.addLayout(layout1)

        self.setLayout(layout)

    def select_odf_phase_1_func(self):
        path = os.path.split(str(self.odf_phase_1_path.text()))[0]
        filename = QFileDialog.getOpenFileName(self, 'Open ODF of phase 1 File', path)  # (self, 'Open ODF File', '/')
        filename = os.path.normpath(str(filename))
        self.odf_phase_1_path.setText(filename)
        print(self.odf_phase_1_path.text())

    def select_odf_phase_2_func(self):
        path = os.path.split(str(self.odf_phase_2_path.text()))[0]
        filename = QFileDialog.getOpenFileName(self, 'Open ODF of phase 2 File', path)  # (self, 'Open ODF File', '/')
        filename = os.path.normpath(str(filename))
        self.odf_phase_2_path.setText(filename)
        print(self.odf_phase_2_path.text())

    def select_data_file_func(self):
        path = os.path.split(str(self.path_of_data_file.text()))[0]

        filename = QFileDialog.getOpenFileName(self, 'Open data file ', path)  # '/')  # (self, 'Open ODF File', '/')
        filename = os.path.normpath(str(filename))
        self.path_of_data_file.setText(filename)
        data = handle_data.AllData()
        self.phase_keys, self.material = data.just_read_data(filename)
        if len(self.phase_keys) == 2:
            self.odf_phase_1_path.setReadOnly(False)
            self.odf_phase_1_button.setEnabled(True)
            self.select_odf_phase_1.setText('select odf phase 1 (= {}): '.format(self.phase_keys[0]))
            self.phase_key_dict[1] = self.phase_keys[0]
            self.odf_phase_2_path.setReadOnly(False)
            self.odf_phase_2_button.setEnabled(True)
            self.select_odf_phase_2.setText('select odf phase 2 (= {}): '.format(self.phase_keys[1]))
            self.phase_key_dict[2] = self.phase_keys[1]

        if len(self.phase_keys) == 1:
            self.odf_phase_1_path.setReadOnly(False)
            self.odf_phase_1_button.setEnabled(True)
            self.select_odf_phase_1.setText('"select odf phase 1 (= {}):   "'.format(self.phase_keys[0]))
            self.select_odf_phase_1.setText('select odf phase 1 (= {}): '.format(self.phase_keys[0]))
            self.phase_key_dict[1] = self.phase_keys[0]

        print(self.path_of_data_file.text())

    def emit_and_quit(self):
        """
        emit the selected path's
        :return:
        """
        odf1 = self.odf_phase_1_path.text()
        odf2 = self.odf_phase_2_path.text()
        data_file = self.path_of_data_file.text()

        self.emit(SIGNAL("data_dir_list"), (odf1, odf2, data_file, self.phase_key_dict, self.material))
        self.close()


class phase_region_class(object):
    def __init__(self):
        self.region = [[1, 'alpha', []],
                       [2, 'betha', []]]  # phase_region_list, [[phase, [[h, k, l, Tmin, Tmax, double, peak]],...], ...]
        # double=0 if single peak
        # double=1 else
        # if double = 1 peak in [1,2]
        # if peak = 1 use the first of the double peaks
        # else use the second

    def save(self, phase_region_list):
        f = open(".\\phase_region_list.dat", "w")
        for i in phase_region_list:
            f.write("phase: " + str(i[0]))
            f.write("\nhkl:    2Theta_min_pos:          2Theta_max_pos:            double:        peak:\n")
            # print("i: ", i)
            for j in i[1]:
                # print("hkl: ", j)
                # print(j[0], j[1], j[2], j[3], j[4], j[5], j[6])
                try:
                    string = "{} {} {} {} {} {} {}".format(int(j[0]), int(j[1]), int(j[2]), j[3], j[4],
                                                           int(float(j[5])), int(float(j[6])))
                    f.write(string + "\n")
                    # print("string: ", string)
                except IndexError:
                    pass
            f.write("########\n")
        f.close()

    def load(self):
        f = open(".\\phase_region_list.dat", "r")
        lines = f.readlines()
        phase = 0
        hkl_2_theta = []
        for i in lines:
            i.strip()
            split = i.split(" ")

            # print(split)
            if split[0] == "phase:":
                phase = int(split[1])
                # self.region[phase-1][0]
                self.region[phase - 1][1] = str(split[2])[0:-1]
                # test = (str(split[2][0:-1])=="BCC")

            elif split[0] == 'Material:':
                self.material = split[1].strip()

            elif split[0] == "hkl:":
                # print("hkl: ")
                pass

            elif split[0] == "########\n":
                self.region[phase - 1][2] = hkl_2_theta
                hkl_2_theta = []

            else:
                # print("else_case: ", split)
                h = int(split[0])
                k = int(split[1])
                l = int(split[2])
                theta_min = int(split[3])
                theta_max = int(split[4])
                double = int(float(split[5]))
                peak = int(float(split[6]))
                hkl_2_theta.append([h, k, l, theta_min, theta_max, double, peak])
        print("REGION: ", self.region)
        return self.region


def main(Thread_queue):
    Thread_queue = Thread_queue
    app = QApplication(sys.argv)
    form = Main()
    form.show()
    app.exec_()
