# import our libraries
import sys
import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc
import PyQt5.QtGui as qg

class AdvancedOptions(qw.QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = "Edit Advanced Options"
        self.left = 150
        self.top = 150
        self.width = 400
        self.height = 300
        self.initEditAdvancedOptions()

    def initEditAdvancedOptions(self):

        # make the window
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # make central widget
        self.wid = qw.QWidget(self)
        self.setCentralWidget(self.wid)

        # create layout grid
        layout = qw.QGridLayout()

        # title at the top of the screen
        title = qw.QLabel("Advanced Options")
        title.setFont(qg.QFont("Times", 20, 2))

        # highly variable label

        # highly variable checkbox
        self.hv_label = qw.QCheckBox("Use highly variable")
        self.hv_label.setFont(qg.QFont("Times", 11, 2))
        self.hv_label.setChecked(True)

        # ndarray checkbox


        # rho label
        rho = qw.QLabel("Rho value")
        rho.setFont(qg.QFont("Times", 11))

        # rho input box
        self.rho_box = qw.QLineEdit(self)
        self.rho_box.setText("0.00001")
        self.rho_box.setFont(qg.QFont("Times", 11))

        # kappa label
        kappa = qw.QLabel("Kappa value")
        kappa.setFont(qg.QFont("Times", 11))

        # kappa input box
        self.kappa_box = qw.QLineEdit(self)
        self.kappa_box.setText("0.00001")
        self.kappa_box.setFont(qg.QFont("Times", 11))

        # delta label
        delta = qw.QLabel("Delta value")
        delta.setFont(qg.QFont("Times", 11))

        # delta input box
        self.delta_box = qw.QLineEdit(self)
        self.delta_box.setText("0.1")
        self.delta_box.setFont(qg.QFont("Times", 11))

        # use weights checkbox
        self.weights_label = qw.QCheckBox("Use weights")
        self.weights_label.setFont(qg.QFont("Times", 11, 2))
        self.weights_label.setChecked(True)

        # top genes facot label
        genes = qw.QLabel("Top genes factor")
        genes.setFont(qg.QFont("Times", 11))

        # top genes factor box
        self.genes_box = qw.QLineEdit(self)
        self.genes_box.setText("50")
        self.genes_box.setFont(qg.QFont("Times", 11))

        # ok button
        ok_button = qw.QPushButton("Ok")
        ok_button.clicked.connect(self.ok_press)
        ok_button.setFont(qg.QFont("Times", 11))




        # add widgets to the layout
        layout.addWidget(title, 0, 0, 1, 5, qc.Qt.AlignCenter)
        layout.addWidget(self.hv_label, 1, 0, 1, 3)
        layout.addWidget(rho, 2, 0, 1, 3)
        layout.addWidget(self.rho_box, 2, 3, 1, 2)
        layout.addWidget(delta, 3, 0, 1, 3)
        layout.addWidget(self.delta_box, 3, 3, 1, 2)
        layout.addWidget(kappa, 4, 0, 1,3)
        layout.addWidget(self.kappa_box, 4, 3, 1, 2)
        layout.addWidget(self.weights_label, 5, 0, 1, 2)
        layout.addWidget(genes, 6, 0, 1, 3)
        layout.addWidget(self.genes_box, 6, 3, 1, 2)

        layout.addWidget(ok_button,7,4,1,1,qc.Qt.AlignBottom)


        layout.setRowStretch(5, 100)
        layout.setHorizontalSpacing(80)


        layout.setColumnStretch(0, 80)
        layout.setColumnStretch(1, 80)
        layout.setColumnStretch(2, 80)
        layout.setColumnStretch(3, 80)
        layout.setColumnStretch(4, 80)

        self.wid.setLayout(layout)


    def ok_press(self):
        self.close()


if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = AdvancedOptions()
    sys.exit(app.exec_())
