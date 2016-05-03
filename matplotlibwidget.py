# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 10:18:04 2015

@author: tneuwirt
"""
from PyQt4.QtGui import QSizePolicy
from PyQt4.QtCore import QSize
from PyQt4 import QtCore, QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
# from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import os
from matplotlib import rcParams
from matplotlib.patches import Rectangle
import numpy as np

rcParams['font.size'] = 9

def get_index_from_arry(array, value):
    array= abs(np.array(array) - value)
    m = min(array)
    for i in xrange(len(array)):
        print abs(abs(m)-abs(array[i]))
        if abs(abs(m)-abs(array[i])) < np.power(10., -10):
            return i
    return -1


class MatplotlibWidget(FigureCanvas):
    """
    MatplotlibWidget inherits PyQt4.QtGui.QWidget
    and matplotlib.backend_bases.FigureCanvasBase

    Options: option_name (default_value)
    -------
    parent (None): parent widget
    title (''): figure title
    xlabel (''): X-axis label
    ylabel (''): Y-axis label
    xlim (None): X-axis limits ([min, max])
    ylim (None): Y-axis limits ([min, max])
    xscale ('linear'): X-axis scale
    yscale ('linear'): Y-axis scale
    width (4): width in inches
    height (3): height in inches
    dpi (100): resolution in dpi
    hold (False): if False, figure will be cleared each time plot is called

    Widget attributes:
    -----------------
    figure: instance of matplotlib.figure.Figure
    axes: figure axes

    Example:
    -------
    self.widget = MatplotlibWidget(self, yscale='log', hold=True)
    from numpy import linspace
    x = linspace(-10, 10)
    self.widget.axes.plot(x, x**2)
    self.wdiget.axes.plot(x, x**3)
    """

    def __init__(self, parent=None, title='', xlabel='', ylabel='',
                 xlim=None, ylim=None, xscale='linear', yscale='linear',
                 width=3, height=3, dpi=100, hold=False):
        self.figure = Figure(figsize=(width, height), dpi=dpi)

        self.axes = self.figure.add_subplot(111)
        self.axes.set_title(title)
        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)
        if xscale is not None:
            self.axes.set_xscale(xscale)
        if yscale is not None:
            self.axes.set_yscale(yscale)
        if xlim is not None:
            self.axes.set_xlim(*xlim)
        if ylim is not None:
            self.axes.set_ylim(*ylim)
        self.axes.hold(hold)

        self.canvas = FigureCanvas.__init__(self, self.figure)

        self.setParent(parent)

        # Canvas.setSizePolicy(self, QSizePolicy.Preferred, QSizePolicy.Preferred)
        policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        policy.setHeightForWidth(True)
        FigureCanvas.setSizePolicy(self, policy)
        FigureCanvas.updateGeometry(self)

    def heightForWidth(self, width):
        return width

    def widthForHeight(self, height):
        return height

    def sizeHint(self):
        w, h = self.get_width_height()
        return QSize(w, h)

    def minimumSizeHint(self):
        return QSize(10, 10)


class MyTable(QtGui.QTableWidget):
    def __init__(self, data, *args):
        QtGui.QTableWidget.__init__(self, *args)
        self.data = data
        self.setmydata()
        self.resizeColumnsToContents()
        self.resizeRowsToContents()

    def setmydata(self):

        horHeaders = []
        for n, key in enumerate(sorted(self.data.keys())):
            horHeaders.append(key)
            for m, item in enumerate(self.data[key]):
                newitem = QtGui.QTableWidgetItem(item)
                self.setItem(m, n, newitem)
        self.setHorizontalHeaderLabels(horHeaders)


class hkl_input_class(QtGui.QWidget):
    def __init__(self, start=0, *args):
        QtGui.QWidget.__init__(self, *args)
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.h = QtGui.QLabel("h")
        self.k = QtGui.QLabel("k")
        self.l = QtGui.QLabel("l")
        self.Theta_min = QtGui.QLabel("T_s")
        self.Theta_max = QtGui.QLabel("T_l")
        self.layout.addWidget(self.h, 0, 1)
        self.layout.addWidget(self.k, 0, 2)
        self.layout.addWidget(self.l, 0, 3)
        self.layout.addWidget(self.Theta_min, 0, 4)
        self.layout.addWidget(self.Theta_max, 0, 5)

        # t = QtGui.QLineEdit("some text")
        # t.setReadOnly(True)
        # t.changeEvent
        # t.setTextBackgroundColor(QtGui.QColor(255, 200, 0))

        self.lines = []
        if start != 0:
            start -= 1

        for i in xrange(start, start + 5):
            line = []
            line.append(QtGui.QLabel("#" + str(i + 1)))
            for j in xrange(5):
                line.append(QtGui.QLineEdit("-1"))
            self.lines.append(line)
        m = 1
        for i in self.lines:
            n = 0
            for j in i:
                try:
                    j.setReadOnly(True)
                except AttributeError:
                    pass
                j.setStyleSheet("color: rgb(255, 100, 0);")
                self.layout.addWidget(j, m, n)
                n += 1
            m += 1

    def get_last_1(self):
        for i in range(0, len(self.lines)):
            if self.lines[i][1].text() == "-1":
                for j in range(1, 4):
                    self.lines[i][j].setReadOnly(False)
                    self.lines[i][j].setStyleSheet("color: rgb(0, 0, 0);")
                for j in range(4, 6):
                    self.lines[i][j].setStyleSheet("color: rgb(0, 0, 0);")
                return i
            for j in range(0, 6):
                print self.lines[i][j].text()
                try:
                    self.lines[i][j].setReadOnly(True)
                except AttributeError:
                    pass
                self.lines[i][j].setStyleSheet("color: rgb(255, 100, 0);")
        return -1


class Preview(QtGui.QWidget):
    # roiSignal = QtCore.pyqtSignal()
    def __init__(self, name):
        QtGui.QWidget.__init__(self)
        QtGui.QWidget.resize(self, 700, 700)
        QtGui.QWidget.setWindowTitle(self, name)
        self.data_list = []
        self.ob_list = []
        self.dc_list = []
        self.filter_list = []
        self.load_filter_list = []
        self.img_filter_list = []
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.roi_Button = QtGui.QPushButton("ROI")
        self.roi_Button.setToolTip("Check to create a ROI in the image. Uncheck to fix the ROI")
        self.roi_Button.setCheckable(True)
        self.oscillation_Button = QtGui.QPushButton("Oscillation")
        self.oscillation_Button.setEnabled(False)
        self.oscillation_Button.setToolTip(
            "Plots the oscillation of the Data images and the OB images in the chosen ROI")

        self.toolbar.addWidget(self.roi_Button)
        self.toolbar.addWidget(self.oscillation_Button)

        # # Tabel with the hkl's in it
        data = {'col1': ['1', '2', '3'], 'col2': ['4', '5', '6'], 'col3': ['7', '8', '9']}
        self.table1_5 = hkl_input_class(start=0)  # MyTable(data=data)
        self.table6_10 = hkl_input_class(start=6)

        newitem = QtGui.QTableWidgetItem("hallo")
        # self.table.setItem(0, 0, newitem)
        self.add_hkl_button = QtGui.QPushButton("add hkl")
        self.add_hkl_button.clicked.connect(self.add_reflex)
        self.add_hkl_button.setEnabled(False)
        self.send_hkl_list_button = QtGui.QPushButton("finish")
        self.send_hkl_list_button.clicked.connect(self.create_hkl_list)

        self.roi_Button.clicked.connect(self.roi)
        self.oscillation_Button.clicked.connect(self.oscillation)

        self.ax = self.figure.add_subplot(111)

        # set the layout1
        layout = QtGui.QVBoxLayout()
        layout3_1 = QtGui.QVBoxLayout()
        layout3 = QtGui.QHBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        # layout.addWidget(self.table)

        # for i in self.table:
        layout3_1.addWidget(self.add_hkl_button)
        layout3_1.addWidget(self.send_hkl_list_button)
        # layout3.addLayout(layout3_1)
        layout3.addWidget(self.table1_5)
        layout3.addWidget(self.table6_10)
        layout3.addLayout(layout3_1)
        # layout3.addWidget(self.add_hkl_button)


        horizontal_line = QtGui.QFrame()
        horizontal_line.setFrameStyle(QtGui.QFrame.HLine)
        horizontal_line.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Expanding)

        layout.addLayout(layout3)
        layout.addWidget(horizontal_line)
        # layout.addLayout(layout2)
        self.setLayout(layout)

    def add_reflex(self):
        print("add")
        table2 = False
        buffer = self.table1_5.get_last_1()

        if buffer == -1:
            table2 = True
            buffer = self.table6_10.get_last_1()

        Theta_min = self.roi_list[2]
        Theta_max = self.roi_list[3]
        print Theta_min, Theta_max, buffer
        if not table2:
            self.table1_5.lines[buffer][4].setText(str(Theta_min))
            self.table1_5.lines[buffer][5].setText(str(Theta_max))

        else:
            self.table6_10.lines[buffer][4].setText(str(Theta_min))
            self.table6_10.lines[buffer][5].setText(str(Theta_max))
        self.dis_hkl_but()

    def dis_hkl_but(self):
        self.add_hkl_button.setEnabled(False)

    def create_hkl_list(self):
        hkl_setting = []
        table2 = False
        buffer = self.table1_5.get_last_1()
        if buffer == -1:
            table2 = True
            buffer = self.table6_10.get_last_1()

        if not table2:
            for i in xrange(0, buffer):
                l = []
                for j in xrange(1, 6):
                    l.append(float(self.table1_5.lines[i][j].text()))
                hkl_setting.append(l)
        else:
            for i in xrange(0, len(self.table1_5.lines)):
                l = []
                for j in xrange(1, 6):
                    l.append(float(self.table1_5.lines[i][j].text()))
                hkl_setting.append(l)
            if buffer==-1:
                buffer=len(self.table6_10.lines)

            for i in xrange(0, buffer):
                l = []
                for j in xrange(1, 6):
                    l.append(float(self.table6_10.lines[i][j].text()))
                hkl_setting.append(l)
        for i in xrange(len(hkl_setting)):
            min = int(get_index_from_arry(self.data_x, hkl_setting[i][3]))
            max = int(get_index_from_arry(self.data_x, hkl_setting[i][4]))
            hkl_setting[i][3] = min
            hkl_setting[i][4] = max

        self.emit(QtCore.SIGNAL("phase_peak_region"), hkl_setting)


    def roi(self):

        self.oscillation_Button.setEnabled(False)
        self.add_hkl_button.setEnabled(False)

        if self.roi_Button.isChecked() == True:
            ax = self.figure.add_subplot(111)
            plt.set_cmap('gray')
            self.ax.hold(False)
            self.ax.plot(self.data_x, self.data_y)
            # if self.typeCombo.currentText() == "Data Images":
            #     ax.imshow(self.data_list[1][self.data_list[0].index(str(self.imgCombo.currentText()))],
            #               vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())
            #
            # elif self.typeCombo.currentText() == "OB Images":
            #     ax.imshow(self.ob_list[1][self.ob_list[0].index(str(self.imgCombo.currentText()))],
            #               vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())
            #
            # elif self.typeCombo.currentText() == "DC Images":
            #     ax.imshow(self.dc_list[1][self.dc_list[0].index(str(self.imgCombo.currentText()))],
            #               vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())

            self.roi_create = ROI(ax)
            # self.typeCombo.setEnabled(False)
            # self.imgCombo.setEnabled(False)

        else:
            roi_x0, roi_x1, roi_y0, roi_y1 = self.roi_create._exit()
            print roi_x0, roi_x1
            len_data_l = len(self.data_list)
            len_ob_l = len(self.ob_list)
            roi_size = roi_x1 - roi_x0
            # if len_data_l is not 0 and len_ob_l is not 0 and roi_size is not 0 and self.roi_Button.isChecked() == False:
            #     self.oscillation_Button.setEnabled(True)
            # region_data = np.array(self.data_list[1][0][roi_y0:roi_y1, roi_x0:roi_x1])
            # print region_data
            # print np.median(region_data), type(np.median(region_data))
            self.roi_list = [roi_y0, roi_y1, roi_x0, roi_x1]
            # self.vmaxSpinBox.setValue(int(1.2 * int(np.max(region_data))))

            # self.roiSignal.emit(roi_list)
            self.emit(QtCore.SIGNAL('roi'), self.roi_list)
            self.add_hkl_button.setEnabled(True)
            # self.typeCombo.setEnabled(True)
            # self.imgCombo.setEnabled(True)

    def oscillation(self):

        len_data_l = len(self.data_list)
        len_ob_l = len(self.ob_list)
        file_number = 0
        oscillation_list = [[], [], []]
        if len_data_l is not 0 and len_ob_l is not 0:
            len_data = len(self.data_list[1])
            len_ob = len(self.ob_list[1])
            while (file_number < len_data):
                if (file_number < len_data):
                    region_data = np.array(self.data_list[1][file_number][self.roi_list[0]:self.roi_list[1],
                                           self.roi_list[2]:self.roi_list[3]])
                    region_ob = np.array(self.ob_list[1][file_number][self.roi_list[0]:self.roi_list[1],
                                         self.roi_list[2]:self.roi_list[3]])
                    av_counts_data = np.median(region_data)
                    av_counts_ob = np.median(region_ob)
                    oscillation_list[0].append(file_number)
                    oscillation_list[1].append(av_counts_data)
                    oscillation_list[2].append(av_counts_ob)
                    file_number += 1
            # print oscillation_list
            rcParams['toolbar'] = 'None'
            ax1 = plt.figure("Oscillation")
            ax2 = ax1.add_subplot(111)
            # ax2.hold(False)
            ax2.plot(oscillation_list[0], oscillation_list[1], 'r*-', label='Data oscillation')
            ax2.plot(oscillation_list[0], oscillation_list[2], 'b*-', label='OB oscillation')
            ax2.legend()
            # self.canvas.draw()
            ax1.show()

            ax1.canvas.manager.window.activateWindow()
            ax1.canvas.manager.window.raise_()

            # plt.close(ax1)
            # rcParams['toolbar']='None'
        else:
            self.dialog = QMessageBox(self)
            self.dialog.setStandardButtons(QMessageBox.Ok)
            self.dialog.setIcon(QMessageBox.Warning)
            self.dialog.setText("Pleas load both data files and OB files to plot the oscillation")
            self.dialog.exec_()

    def choose_type(self):

        if self.typeCombo.currentText() == "Data Images":
            self.imgCombo.clear()
            self.imgCombo.addItems(self.data_list[0])

        elif self.typeCombo.currentText() == "OB Images":
            self.imgCombo.clear()
            self.imgCombo.addItems(self.ob_list[0])
        elif self.typeCombo.currentText() == "DC Images":
            self.imgCombo.clear()
            self.imgCombo.addItems(self.dc_list[0])
        elif self.typeCombo.currentText() == "Filtered Images":
            self.imgCombo.clear()
            self.imgCombo.addItems(self.filter_list[0])
        # self.choose_img()
        len_data_l = len(self.data_list)
        len_ob_l = len(self.ob_list)

        # roi_size=roi_x1-roi_x0
        if len_data_l is not 0 and len_ob_l is not 0 and self.roi_Button.isChecked() == False:
            self.oscillation_Button.setEnabled(True)

    def choose_img(self):

        ax = self.figure.add_subplot(111)
        plt.set_cmap('gray')
        ax.hold(False)

        if self.typeCombo.currentText() == "Data Images":
            ax.imshow(self.data_list[1][self.data_list[0].index(str(self.imgCombo.currentText()))],
                      vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())

        elif self.typeCombo.currentText() == "OB Images":
            ax.imshow(self.ob_list[1][self.ob_list[0].index(str(self.imgCombo.currentText()))],
                      vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())

        elif self.typeCombo.currentText() == "DC Images":
            ax.imshow(self.dc_list[1][self.dc_list[0].index(str(self.imgCombo.currentText()))],
                      vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())

        elif self.typeCombo.currentText() == "Filtered Images":
            ax.imshow(self.filter_list[1][self.filter_list[0].index(str(self.imgCombo.currentText()))],
                      vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())
        self.rect = Rectangle((0, 0), 1, 1, facecolor='None', edgecolor='green')
        self.rect.set_width(self.roi_list[3] - self.roi_list[2])
        self.rect.set_height(self.roi_list[1] - self.roi_list[0])
        self.rect.set_xy((self.roi_list[2], self.roi_list[0]))
        self.rect.set_linestyle('solid')
        ax.add_patch(self.rect)

        # refresh canvas
        self.canvas.draw()

    def add_data(self, load_data_list, data_img_list, roi_list):
        self.roi_list = roi_list
        self.imgCombo.clear()
        load_data_list_temp = []
        for i in range(0, len(load_data_list)):
            temp_ind = load_data_list[i].rfind(str(os.path.sep))
            load_data_list_temp.append(load_data_list[i][temp_ind:])
        self.data_list = [load_data_list_temp, data_img_list]
        self.imgCombo.addItems(load_data_list_temp)
        self.typeCombo.removeItem(self.typeCombo.findText("Data Images"))
        self.typeCombo.addItem("Data Images")
        self.vminSpinBox.setValue(0.75 * data_img_list[0].min())
        self.vmaxSpinBox.setValue(1.25 * data_img_list[0].max())
        self.choose_type()

    def add_xy_data(self, data_x, data_y):
        self.roi_list = [0, 0, 0, 0]
        self.data_x = data_x
        self.data_y = data_y
        # ax = self.figure.add_subplot(111)
        self.ax.plot(self.data_x, self.data_y)

    def add_ob(self, load_ob_list, ob_img_list):
        self.imgCombo.clear()
        load_ob_list_temp = []
        for i in range(0, len(load_ob_list)):
            temp_ind = load_ob_list[i].rfind(str(os.path.sep))
            load_ob_list_temp.append(load_ob_list[i][temp_ind:])
        self.ob_list = [load_ob_list_temp, ob_img_list]
        self.imgCombo.addItems(load_ob_list_temp)
        self.typeCombo.removeItem(self.typeCombo.findText("OB Images"))
        self.typeCombo.addItem("OB Images")
        self.choose_type()

    def add_dc(self, load_dc_list, dc_img_list, dc_median):
        self.imgCombo.clear()
        load_dc_list_temp = []
        for i in range(0, len(load_dc_list)):
            temp_ind = load_dc_list[i].rfind(str(os.path.sep))
            load_dc_list_temp.append(load_dc_list[i][temp_ind:])
        if dc_median is not None:
            load_dc_list_temp.append("Median of DC")
            dc_img_list.append(dc_median)
        self.dc_list = [load_dc_list_temp, dc_img_list]
        self.imgCombo.addItems(load_dc_list_temp)
        self.typeCombo.removeItem(self.typeCombo.findText("DC Images"))
        self.typeCombo.addItem("DC Images")
        self.choose_type()

    def add_filtered(self, test_img, img):
        self.imgCombo.clear()
        temp_ind = test_img.rfind(str(os.path.sep))
        self.load_filter_list.append(test_img[temp_ind:])
        print self.load_filter_list
        self.img_filter_list.append(img)
        self.filter_list = [self.load_filter_list, self.img_filter_list]
        self.imgCombo.addItems(self.load_filter_list)
        self.typeCombo.removeItem(self.typeCombo.findText("Filtered Images"))
        self.typeCombo.addItem("Filtered Images")
        self.choose_type()


class ROI(object):
    def __init__(self, ax1, ax2=None):
        # global roi_x0
        # global roi_x1
        # global roi_y0
        # global roi_y1
        self.ax1 = ax1
        self.ax2 = ax2
        self.rect = Rectangle((0, 0), 1, 1, facecolor='None', edgecolor='green')
        # self.x0 = roi_x0
        # self.y0 = roi_y0
        # self.x1 = roi_x1
        # self.y1 = roi_y1
        # self.ax1=plt.gca()
        self.ax1.add_patch(self.rect)
        if self.ax2 != None:
            self.ax2.add_patch(self.rect)
        # self.is_pressed=True
        self.ax1_t1 = self.ax1.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax1_t2 = self.ax1.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.ax1_t3 = self.ax1.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        if self.ax2 != None:
            self.ax2_t1 = self.ax2.figure.canvas.mpl_connect('button_press_event', self.on_press)
            self.ax2_t2 = self.ax2.figure.canvas.mpl_connect('button_release_event', self.on_release)
            self.ax2_t3 = self.ax2.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.is_pressed = False

    def on_press(self, event):

        self.x0 = event.xdata
        self.y0 = event.ydata
        self.x1 = event.xdata
        self.y1 = event.ydata
        print 'press', event.xdata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.rect.set_linestyle('dashed')
        self.ax1.figure.canvas.draw()
        if self.ax2 != None:
            self.ax2.figure.canvas.draw()

        self.is_pressed = True

    def on_motion(self, event):
        # if self.on_press is True:
        # return
        if self.is_pressed == True:
            self.x1 = event.xdata
            self.y1 = event.ydata
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.rect.set_linestyle('dashed')
            self.ax1.figure.canvas.draw()
            if self.ax2 != None:
                self.ax2.figure.canvas.draw()

    def on_release(self, event):
        # global roi_x0
        # global roi_x1
        # global roi_y0
        # global roi_y1
        print 'release', event.xdata
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.rect.set_linestyle('solid')

        self.ax1.figure.canvas.draw()
        if self.ax2 != None:
            self.ax2.figure.canvas.draw()
        self.is_pressed = False

        # roi_x0 =int(min(self.x0,self.x1))
        # roi_x1 =int(max(self.x0,self.x1))
        # roi_y0 =int(min(self.y0,self.y1))
        # roi_y1 =int(max(self.y0,self.y1))
        # self.roi_list=[roi_y0,roi_y1,roi_x0,roi_x1]
        # self.emit(QtCore.SIGNAL('release'),self.roi_list)


        # print self.x0,self.x1,self.y0,self.y1
        # return [self.x0,self.x1,self.y0,self.y1]

    def _exit(self):
        self.ax1.figure.canvas.mpl_disconnect(self.ax1_t1)
        self.ax1.figure.canvas.mpl_disconnect(self.ax1_t2)
        self.ax1.figure.canvas.mpl_disconnect(self.ax1_t3)
        if self.ax2 != None:
            self.ax2.figure.canvas.mpl_disconnect(self.ax2_t1)
            self.ax2.figure.canvas.mpl_disconnect(self.ax2_t2)
            self.ax2.figure.canvas.mpl_disconnect(self.ax2_t3)
        roi_x0 = int(min(self.x0, self.x1))
        roi_x1 = int(max(self.x0, self.x1))
        roi_y0 = int(min(self.y0, self.y1))
        roi_y1 = int(max(self.y0, self.y1))
        return roi_x0, roi_x1, roi_y0, roi_y1


class Filter_Preview(QtGui.QWidget):
    # roiSignal = QtCore.pyqtSignal()
    def __init__(self, name):
        QtGui.QWidget.__init__(self)
        QtGui.QWidget.resize(self, 1000, 600)
        QtGui.QWidget.setWindowTitle(self, name)
        self.data_list = []
        self.ob_list = []
        self.dc_list = []
        self.filter_list = []
        self.load_filter_list = []
        self.img_filter_list = []
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)
        # self.roi_Button=QtGui.QPushButton("ROI")
        # self.roi_Button.setToolTip("Check to create a ROI in the image. Uncheck to fix the ROI")
        # self.roi_Button.setCheckable(True)
        # self.oscillation_Button=QtGui.QPushButton("Oscillation")
        # self.oscillation_Button.setEnabled(False)
        # self.oscillation_Button.setToolTip("Plots the oscillation of the Data images and the OB images in the chosen ROI")

        # self.toolbar.addWidget(self.roi_Button)
        # self.toolbar.addWidget(self.oscillation_Button)


        self.vminSpinBox = QtGui.QSpinBox()
        self.vminLabel = QtGui.QLabel("Min Grayvalue")
        self.vminSpinBox.setRange(0, 64000)
        self.vminSpinBox.setSingleStep(500)

        self.vmaxSpinBox = QtGui.QSpinBox()
        self.vmaxLabel = QtGui.QLabel("Max Grayvalue")
        self.vmaxSpinBox.setRange(0, 64000)
        self.vmaxSpinBox.setValue(32000)
        self.vmaxSpinBox.setSingleStep(500)

        # Just some button connected to `plot` method
        self.typeCombo = QtGui.QComboBox()

        self.imgCombo = QtGui.QComboBox()

        # self.roi_Button.clicked.connect(self.roi)
        # self.oscillation_Button.clicked.connect(self.oscillation)

        self.typeCombo.currentIndexChanged.connect(self.choose_type)
        self.imgCombo.currentIndexChanged.connect(self.choose_img)
        self.vminSpinBox.valueChanged.connect(self.choose_img)
        self.vmaxSpinBox.valueChanged.connect(self.choose_img)
        ax1 = self.figure.add_subplot(221)
        ax2 = self.figure.add_subplot(222)
        ax3 = self.figure.add_subplot(223)

        # set the layout1
        layout = QtGui.QVBoxLayout()
        layout2 = QtGui.QHBoxLayout()
        layout3 = QtGui.QHBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        layout2.addWidget(self.typeCombo)
        layout2.addWidget(self.imgCombo)
        layout3.addWidget(self.vminLabel)
        layout3.addWidget(self.vminSpinBox)
        layout3.addWidget(self.vmaxLabel)
        layout3.addWidget(self.vmaxSpinBox)
        layout.addLayout(layout3)
        layout.addLayout(layout2)
        self.setLayout(layout)

    def roi(self):

        self.oscillation_Button.setEnabled(False)

        if self.roi_Button.isChecked() == True:
            ax = self.figure.add_subplot(111)
            plt.set_cmap('gray')
            ax.hold(False)
            if self.typeCombo.currentText() == "Data Images":
                ax.imshow(self.data_list[1][self.data_list[0].index(str(self.imgCombo.currentText()))],
                          vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())

            elif self.typeCombo.currentText() == "OB Images":
                ax.imshow(self.ob_list[1][self.ob_list[0].index(str(self.imgCombo.currentText()))],
                          vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())

            elif self.typeCombo.currentText() == "DC Images":
                ax.imshow(self.dc_list[1][self.dc_list[0].index(str(self.imgCombo.currentText()))],
                          vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())

            self.roi_create = ROI(ax)
            self.typeCombo.setEnabled(False)
            self.imgCombo.setEnabled(False)

        else:
            roi_x0, roi_x1, roi_y0, roi_y1 = self.roi_create._exit()
            len_data_l = len(self.data_list)
            len_ob_l = len(self.ob_list)
            roi_size = roi_x1 - roi_x0
            if len_data_l is not 0 and len_ob_l is not 0 and roi_size is not 0 and self.roi_Button.isChecked() == False:
                self.oscillation_Button.setEnabled(True)
            region_data = np.array(self.data_list[1][0][roi_y0:roi_y1, roi_x0:roi_x1])
            # print region_data
            # print np.median(region_data), type(np.median(region_data))
            self.roi_list = [roi_y0, roi_y1, roi_x0, roi_x1]
            self.vmaxSpinBox.setValue(int(1.2 * int(np.max(region_data))))

            # self.roiSignal.emit(roi_list)
            self.emit(QtCore.SIGNAL('roi'), self.roi_list)
            self.typeCombo.setEnabled(True)
            self.imgCombo.setEnabled(True)

    def oscillation(self):

        len_data_l = len(self.data_list)
        len_ob_l = len(self.ob_list)
        file_number = 0
        oscillation_list = [[], [], []]
        if len_data_l is not 0 and len_ob_l is not 0:
            len_data = len(self.data_list[1])
            len_ob = len(self.ob_list[1])
            while (file_number < len_data):
                if (file_number < len_data):
                    region_data = np.array(self.data_list[1][file_number][self.roi_list[0]:self.roi_list[1],
                                           self.roi_list[2]:self.roi_list[3]])
                    region_ob = np.array(self.ob_list[1][file_number][self.roi_list[0]:self.roi_list[1],
                                         self.roi_list[2]:self.roi_list[3]])
                    av_counts_data = np.median(region_data)
                    av_counts_ob = np.median(region_ob)
                    oscillation_list[0].append(file_number)
                    oscillation_list[1].append(av_counts_data)
                    oscillation_list[2].append(av_counts_ob)
                    file_number += 1
            # print oscillation_list
            rcParams['toolbar'] = 'None'
            ax1 = plt.figure("Oscillation")
            ax2 = ax1.add_subplot(111)
            # ax2.hold(False)
            ax2.plot(oscillation_list[0], oscillation_list[1], 'r*-', label='Data oscillation')
            ax2.plot(oscillation_list[0], oscillation_list[2], 'b*-', label='OB oscillation')
            ax2.legend()
            # self.canvas.draw()
            ax1.show()

            ax1.canvas.manager.window.activateWindow()
            ax1.canvas.manager.window.raise_()
            # plt.close(ax1)
            # rcParams['toolbar']='None'
        else:
            self.dialog = QMessageBox(self)
            self.dialog.setStandardButtons(QMessageBox.Ok)
            self.dialog.setIcon(QMessageBox.Warning)
            self.dialog.setText("Pleas load both data files and OB files to plot the oscillation")
            self.dialog.exec_()

    def choose_type(self):

        if self.typeCombo.currentText() == "Data Images":
            self.imgCombo.clear()
            self.imgCombo.addItems(self.data_list[0])

        elif self.typeCombo.currentText() == "OB Images":
            self.imgCombo.clear()
            self.imgCombo.addItems(self.ob_list[0])
        elif self.typeCombo.currentText() == "DC Images":
            self.imgCombo.clear()
            self.imgCombo.addItems(self.dc_list[0])
        elif self.typeCombo.currentText() == "Filtered Images":
            self.imgCombo.clear()
            self.imgCombo.addItems(self.filter_list[0])
        # self.choose_img()
        len_data_l = len(self.data_list)
        len_ob_l = len(self.ob_list)

        # roi_size=roi_x1-roi_x0
        if len_data_l is not 0 and len_ob_l is not 0 and self.roi_Button.isChecked() == False:
            self.oscillation_Button.setEnabled(True)

    def choose_img(self):

        ax = self.figure.add_subplot(111)
        plt.set_cmap('gray')
        ax.hold(False)

        if self.typeCombo.currentText() == "Data Images":
            ax.imshow(self.data_list[1][self.data_list[0].index(str(self.imgCombo.currentText()))],
                      vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())

        elif self.typeCombo.currentText() == "OB Images":
            ax.imshow(self.ob_list[1][self.ob_list[0].index(str(self.imgCombo.currentText()))],
                      vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())

        elif self.typeCombo.currentText() == "DC Images":
            ax.imshow(self.dc_list[1][self.dc_list[0].index(str(self.imgCombo.currentText()))],
                      vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())

        elif self.typeCombo.currentText() == "Filtered Images":
            ax.imshow(self.filter_list[1][self.filter_list[0].index(str(self.imgCombo.currentText()))],
                      vmin=self.vminSpinBox.value(), vmax=self.vmaxSpinBox.value())
        self.rect = Rectangle((0, 0), 1, 1, facecolor='None', edgecolor='green')
        self.rect.set_width(self.roi_list[3] - self.roi_list[2])
        self.rect.set_height(self.roi_list[1] - self.roi_list[0])
        self.rect.set_xy((self.roi_list[2], self.roi_list[0]))
        self.rect.set_linestyle('solid')
        ax.add_patch(self.rect)

        # refresh canvas
        self.canvas.draw()

    def add_data(self, load_data_list, data_img_list, roi_list):
        self.roi_list = roi_list
        self.imgCombo.clear()
        load_data_list_temp = []
        for i in range(0, len(load_data_list)):
            temp_ind = load_data_list[i].rfind(str(os.path.sep))
            load_data_list_temp.append(load_data_list[i][temp_ind:])
        self.data_list = [load_data_list_temp, data_img_list]
        self.imgCombo.addItems(load_data_list_temp)
        self.typeCombo.removeItem(self.typeCombo.findText("Data Images"))
        self.typeCombo.addItem("Data Images")
        self.vminSpinBox.setValue(0.75 * data_img_list[0].min())
        self.vmaxSpinBox.setValue(1.25 * data_img_list[0].max())
        self.choose_type()

    def add_ob(self, load_ob_list, ob_img_list):
        self.imgCombo.clear()
        load_ob_list_temp = []
        for i in range(0, len(load_ob_list)):
            temp_ind = load_ob_list[i].rfind(str(os.path.sep))
            load_ob_list_temp.append(load_ob_list[i][temp_ind:])
        self.ob_list = [load_ob_list_temp, ob_img_list]
        self.imgCombo.addItems(load_ob_list_temp)
        self.typeCombo.removeItem(self.typeCombo.findText("OB Images"))
        self.typeCombo.addItem("OB Images")
        self.choose_type()

    def add_dc(self, load_dc_list, dc_img_list, dc_median):
        self.imgCombo.clear()
        load_dc_list_temp = []
        for i in range(0, len(load_dc_list)):
            temp_ind = load_dc_list[i].rfind(str(os.path.sep))
            load_dc_list_temp.append(load_dc_list[i][temp_ind:])
        if dc_median is not None:
            load_dc_list_temp.append("Median of DC")
            dc_img_list.append(dc_median)
        self.dc_list = [load_dc_list_temp, dc_img_list]
        self.imgCombo.addItems(load_dc_list_temp)
        self.typeCombo.removeItem(self.typeCombo.findText("DC Images"))
        self.typeCombo.addItem("DC Images")
        self.choose_type()

    def add_filtered(self, test_img, img):
        self.imgCombo.clear()
        temp_ind = test_img.rfind(str(os.path.sep))
        self.load_filter_list.append(test_img[temp_ind:])
        print self.load_filter_list
        self.img_filter_list.append(img)
        self.filter_list = [self.load_filter_list, self.img_filter_list]
        self.imgCombo.addItems(self.load_filter_list)
        self.typeCombo.removeItem(self.typeCombo.findText("Filtered Images"))
        self.typeCombo.addItem("Filtered Images")
        self.choose_type()
