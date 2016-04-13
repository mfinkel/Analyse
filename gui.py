import sys
from PyQt4.QtGui import *

class gui(QWidget):
    def __init__(self, name):
        QWidget.__init__(self)
        QWidget.resize(self, 700, 700)
        QWidget.setWindowTitle(self, name)
        # return filename
        # Create an PyQT4 application object.
        a = QApplication(sys.argv)

        # The QWidget widget is the base class of all user interface objects in PyQt4.
        self.widgate = QMainWindow()
        self.widgate1 = QWidget()
        #self.string_class = tacke_string(self.textbox.text(), self.widgate1)
        # Set window size.
        self.widgate.resize(320, 240)
        self.s = ""

        # Set window title
        self.widgate.setWindowTitle("Hello World!")

        # Create main menu
        self.mainMenu = self.widgate.menuBar()
        self.mainMenu.setNativeMenuBar(True)
        self.fileMenu = self.mainMenu.addMenu('&amp;File')
        self.toolbar = self.widgate.NavigationToolbar(True)


        # Add exit button
        self.exitButton = QAction(QIcon('exit24.png'), 'Exit', self.widgate)
        self.exitButton.setShortcut('Ctrl+Q')
        self.exitButton.setStatusTip('Exit application')
        self.exitButton.triggered.connect(self.widgate.close)
        self.fileMenu.addAction(self.exitButton)

        # Create textbox
        self.textbox = QLineEdit(self.widgate)
        self.textbox.move(20, 20)
        self.textbox.resize(280, 40)

        # Create combobox
        self.combo = QComboBox(self.widgate)
        self.combo.addItem("isotrop")
        self.combo.addItem("cubic")
        self.combo.addItem("hexagonal")
        self.combo.addItem("tetragonal_1 (not implementet jet)")
        self.combo.addItem("tetragonal_2 (not implementet jet)")
        self.sym_matrix.addAction(self.combo)
        self.combo.move(20, 80)

        # Create a button in the window
        button = QPushButton('Click me', self.widgate)
        button.move(300, 20)
        button1 = QPushButton('call goo', self.widgate)
        button1.move(300, 80)
        button2 = QPushButton('call goo', self.widgate)
        button2.move(300, 120)
        button3 = QPushButton('print s', self.widgate)
        button3.move(300, 120)
        self.widgate.resize(400, 300)

        def on_click():
            self.textbox.setText("Button clicked.")

        # connect the signals to the slots
        button.clicked.connect(self.open_file)
        button1.clicked.connect(self.call_class)
        button2.clicked.connect(self.print_)
        button3.clicked.connect(self.prints_)

        # print filename
        # Show window
        self.widgate.show()
        a.exec_()

    def open_file(self):
        filename = QFileDialog.getOpenFileName(self.widgate, 'Open File', '/')
        self.textbox.setText(filename)

    def call_class(self):
        self.string_class = self.tacke_string(self.textbox.text(), self.widgate1)

    def print_(self):
        print self.string_class.string
        self.s = self.string_class.get_string()
        print self.string_class.get_string()

    def prints_(self):
        print self.s
