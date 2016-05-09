from __future__ import print_function
import sys
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from PyQt4.QtGui import *
import lmfit

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
from matplotlib.backends import qt4_compat
import matplotlibwidget

use_pyside = qt4_compat.QT_API == qt4_compat.QT_API_PYSIDE

if use_pyside:
    from PySide.QtCore import *
    from PySide.QtGui import *
else:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
import methods
from functools import partial
import handle_data


class insert_const_widget(QWidget):
    def __init__(self, name, params, sym, *args):
        QWidget.__init__(self)
        QWidget.resize(self, 500, 400)
        QWidget.setWindowTitle(self, name)
        self.name = name
        self.layout = QGridLayout()
        self.params = params
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
        self.confirm_button.clicked.connect(self.set_params)
        self.show()

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

        elif self.sym == "cubic":
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
        elif self.sym == "cubic":
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
        self.data_object = handle_data.Data(sample_diameter=None)
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")
        self.region = phase_region_class()
        try:
            self.phase_peak_region = self.region.load()
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
        self.automate = QLabel("select hkl manualy?")
        self.select_hkl_setting_manualy_coice = self.create_jn_combbox()
        self.load_data_button = QPushButton('load Data')
        self.load_data_button.setEnabled(False)

        self.open_data_folders.clicked.connect(self.set_data_path_func)

        # select the experimental setup (SPODI, POLDI, ....)
        self.load_data_button.clicked.connect(self.read_scattering_data_SPODI_case)
        self.choose_experiment_comb_box.currentIndexChanged.connect(self.connect_read_scattering_data)

        # add the central plot to display the data
        self.central_plot = matplotlibwidget.Preview("Select_Data")
        self.connect(self.central_plot, SIGNAL("phase_peak_region"), self.set_hkl_setting_and_fit_the_peaks_SPODI)

        # handel the fitting process
        self.Fit_phase = QLabel("Fit phases: ")
        self.fit_phase_combbox = self.create_n_o_phases_combbox()
        self.modi_text = QLabel("Theory")
        self.modi = self.create_modi_comb_box()
        self.ODF_text = QLabel("Texture?")
        self.text_jn = self.create_jn_combbox()

        self.do_the_fit_button = QPushButton('fitting Data')
        self.do_the_fit_button.setEnabled(False)
        self.do_the_fit_button.clicked.connect(self.fit_the_data)

        self.connect(self, SIGNAL("data"), self.central_plot.add_data)

        self.layout_handling()

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
        elif self.choose_experiment_comb_box.currentText() == "POLDI":
            print('is still to implement')

    def receve_the_pathes_SPODI_case(self, *args):
        """
        receve the data path
        :param args:
        :return:
        """
        print(args, len(args[0]))
        odf1, odf2, unstraind, data_dir_list = args[0]
        odf1 = str(odf1)
        odf2 = str(odf2)
        unstraind = str(unstraind)
        for i in data_dir_list:
            i = str(i)
            print(i, type(i))
        print(odf1)
        if odf2 == "None":
            odf2 = None
        print(odf1, odf2, unstraind, data_dir_list)
        self.path_of_odf_phase1 = odf1
        self.path_of_odf_phase2 = odf2
        self.path_of_unstraind_data = unstraind
        self.path_of_data_under_strain = data_dir_list
        self.load_data_button.setEnabled(True)

    def read_scattering_data_SPODI_case(self):
        Bool = False
        if self.select_hkl_setting_manualy_coice.currentText() == "Yes":
            Bool = True

        self.data_object = handle_data.SPODIData(sample_diameter=int(str(self.diameter.text())),
                                                 odf_phase_1_file=self.path_of_odf_phase1,
                                                 odf_phase_2_file=self.path_of_odf_phase2)
        self.data_object.load_data(self.path_of_data_under_strain)

        # self.Data_Iron = methods.Data_old(str(self.odf_phase_1_path.text()), 6)
        # self.Data_Iron.read_scattering_SPODI_data(path_of_unstraind_data=str(self.path_of_unstraind_data.text()),
        #                                           path_of_straind_data=str(self.path_of_straind_data_1.text()))

        if Bool or not self.loaded_peak_region:
            self.select_hkl_SPODI_Data()
        else:
            self.data_object.fit_all_data(peak_regions_phase=self.phase_peak_region, plot=False)
            self.do_the_fit_button.setEnabled(True)
            # self.Data_Iron.fit_all_peaks()

    def select_hkl_SPODI_Data(self):
        x_data, y_data = self.data_object.get_sum_data()
        print("data_x:", x_data)
        self.central_plot.add_xy_data(x_data, y_data)

    def set_hkl_setting_and_fit_the_peaks_SPODI(self, value):
        self.phase_peak_region = value
        print("------------------------------------")
        print("peak region: ", self.phase_peak_region)

        self.region.save(self.phase_peak_region)
        self.region.load()
        # np.save(".\\phase_peak_region", np.array(self.phase_peak_region))
        print("saved peak region")
        print("------------------------------------")
        self.data_object.fit_all_data(peak_regions_phase=self.phase_peak_region)
        # self.Data_Iron.set_hkl_setting(self.phase_peak_region)
        # self.Data_Iron.fit_all_peaks()

        self.do_the_fit_button.setEnabled(True)
        print("coming from other class:", "\n", self.phase_peak_region)

    def fit_the_data(self):
        Bool = False
        if self.text_jn.currentText() == "Yes":
            Bool = True
        print(self.modi.currentText(), "\n",
              Bool)
        try:
            self.Data_Iron.Fit_the_data_with_texture(filename="Result_iron_", method=str(self.modi.currentText()),
                                                     number_of_datapoints=None, texture=Bool)
        except AttributeError:
            self.read_scattering_data_SPODI_case()
            self.fit_the_data()

    def layout_handling(self):
        # Layout handling
        layout = QVBoxLayout()
        layout_ok_and_cancel_button = QHBoxLayout()
        layout_fitting = QHBoxLayout()

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
        layout_load_data_h2.addWidget(self.label("# of dataset's under strain: "))
        layout_load_data_h2.addWidget(self.number_of_datasets_under_strain)
        layout_load_data_h2.addWidget(self.open_data_folders)
        layout_load_data_h3.addWidget(self.label("set diameter in mm: "))
        layout_load_data_h3.addWidget(self.diameter)
        layout_load_data_h3.addWidget(self.automate)
        layout_load_data_h3.addWidget(self.select_hkl_setting_manualy_coice)
        layout_load_data_h3.addWidget(self.load_data_button)
        layout.addLayout(layout_load_data_h1)
        layout.addLayout(layout_load_data_h2)
        layout.addLayout(layout_load_data_h3)
        layout.addWidget(HLine[2])

        # handel the fitting process
        layout.addWidget(self.label("Fit the data: "))
        layout_fitting.addWidget(self.Fit_phase)
        layout_fitting.addWidget(self.fit_phase_combbox)
        layout_fitting.addWidget(self.modi_text)
        layout_fitting.addWidget(self.modi)
        layout_fitting.addWidget(self.ODF_text)
        layout_fitting.addWidget(self.text_jn)
        layout_fitting.addWidget(self.do_the_fit_button)

        layout.addLayout(layout_fitting)
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
            "..\\Daten-bearbeitet\\Stahl ST37\\" + "ST37_MTODF.txt")  # "AL_textur_complet.txt"
        self.odf_phase_2_path = QLineEdit("None")  # "AL_textur_complet.txt"
        self.path_of_unstraind_data = QLineEdit("..\\Daten-bearbeitet\\Stahl ST37\\" + "Euler-Scans ohne Last\\")

        self.straind_data = []
        self.straind_data_lable = []
        self.straind_data_button = []

        # create boxes for all datafiles of different load
        for i in xrange(number_of_straind_datasets):
            self.straind_data.append(QLineEdit("..\\Daten-bearbeitet\\Stahl ST37\\" + "Euler-Scans unter 5kN\\"))
            self.straind_data_lable.append(QLabel("select straind data %i:" % (i)))
            self.straind_data_button.append(QPushButton("select straind %i:" % (i)))

        self.select_odf_phase_1 = QLabel("select odf phase 1:   ")
        self.select_odf_phase_2 = QLabel("select odf phase 2:   ")
        self.select_unstraind = QLabel("select unstraind data:")

        self.odf_phase_1_button.clicked.connect(self.select_odf_phase_1_func)
        self.odf_phase_2_button.clicked.connect(self.select_odf_phase_2_func)
        self.unstraind_button.clicked.connect(self.select_unstraind_func)

        for i in xrange(len(self.straind_data_button)):
            self.straind_data_button[i].clicked.connect(partial(self.select_straind_func, i))
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
        layout_straind_data = []
        for i in xrange(len(self.straind_data)):
            layout_straind_data.append(QHBoxLayout())
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
        for i in xrange(len(layout_straind_data)):
            layout_straind_data[i].addWidget(self.straind_data_lable[i])
            layout_straind_data[i].addWidget(self.straind_data[i])
            layout_straind_data[i].addWidget(self.straind_data_button[i])

        # insert unstraind data
        layout_unstraind_data.addWidget(self.select_unstraind)
        layout_unstraind_data.addWidget(self.path_of_unstraind_data)
        layout_unstraind_data.addWidget(self.unstraind_button)

        layout.addLayout(layout_odf_phase_1_input)
        layout.addLayout(layout_odf_phase_2_input)
        layout.addLayout(layout_unstraind_data)
        for i in layout_straind_data:
            layout.addLayout(i)
        # layout.addLayout(layout_straind_data_1)
        # layout.addLayout(layout_straind_data_2)
        # layout.addLayout(layout_straind_data_3)
        # layout.addLayout(layout_straind_data_4)

        layout.addLayout(layout1)

        self.setLayout(layout)

    def select_odf_phase_1_func(self):
        filename = QFileDialog.getOpenFileName(self, 'Open ODF of phase 1 File', '/')  # (self, 'Open ODF File', '/')
        filename = os.path.normpath(str(filename))
        self.odf_phase_1_path.setText(filename)
        print(self.odf_phase_1_path.text())

    def select_odf_phase_2_func(self):
        filename = QFileDialog.getOpenFileName(self, 'Open ODF of phase 2 File', '/')  # (self, 'Open ODF File', '/')
        filename = os.path.normpath(str(filename))
        self.odf_phase_2_path.setText(filename)
        print(self.odf_phase_2_path.text())

    def select_unstraind_func(self):
        filename = QFileDialog.getExistingDirectory(self, 'Open unstraind data', '/')  # (self, 'Open ODF File', '/')
        filename = os.path.normpath(str(filename))
        self.path_of_unstraind_data.setText(filename + "\\")
        print(self.path_of_unstraind_data.text())

    def select_straind_func(self, i):
        filename = QFileDialog.getExistingDirectory(self, 'Open straind data', '/')  # (self, 'Open ODF File', '/')
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
        unstraind = self.path_of_unstraind_data.text()
        data_dir_list = []
        data_dir_list.append(self.path_of_unstraind_data.text())
        for i in xrange(len(self.straind_data)):
            data_dir_list.append(self.straind_data[i].text())
        self.emit(SIGNAL("data_dir_list"), (odf1, odf2, unstraind, data_dir_list))
        self.close()


class phase_region_class(object):
    def __init__(self):
        self.region = [[1, []], [2, []]]  # phase_region_list

    def save(self, phase_region_list):
        f = open(".\\phase_region_list.dat", "w")
        for i in phase_region_list:
            f.write("phase: " + str(i[0]))
            f.write("\nhkl:    2Theta_min_pos:          2Theta_max_pos:\n")
            print ("i: ", i)
            for j in i[1]:
                print ("hkl: ", j)
                print (j[0], j[1], j[2], j[3], j[4])
                try:
                    string = "{} {} {} {} {} ".format(int(j[0]), int(j[1]), int(j[2]), j[3], j[4])
                    f.write(string+"\n")
                    print ("string: ", string)
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

            print(split)
            if split[0] == "phase:":
                phase = int(split[1])

            elif split[0] == "hkl:":
                print("hkl: ")
                pass

            elif split[0] == "########\n":
                self.region[phase - 1][1] = hkl_2_theta
                hkl_2_theta = []


            else:
                print ("else_case: ", split)
                h = int(split[0])
                k = int(split[1])
                l = int(split[2])
                theta_min = int(split[3])
                theta_max = int(split[4])
                hkl_2_theta.append([h, k, l, theta_min, theta_max])
        print ("REGION: ", self.region)
        return self.region


def main():
    app = QApplication(sys.argv)
    form = Main()
    form.show()
    app.exec_()
