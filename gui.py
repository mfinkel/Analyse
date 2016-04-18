from __future__ import print_function
import sys
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from PyQt4.QtGui import *
import lmfit

import functools
import sys
import lmfit as lm

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backends import qt4_compat

use_pyside = qt4_compat.QT_API == qt4_compat.QT_API_PYSIDE

if use_pyside:
    from PySide.QtCore import *
    from PySide.QtGui import *
else:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *


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
        self.odf_button.move(20,150)
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

        #Send Data outside
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


class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        # self.x, self.y = self.get_data()
        self.data = self.get_data2()
        self.create_main_frame()
        self.on_draw()

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
        print('you pressed', event.key)
        # implement the default mpl key press events described at
        # http://matplotlib.org/users/navigation_toolbar.html#navigation-keyboard-shortcuts
        key_press_handler(event, self.canvas, self.mpl_toolbar)


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()
