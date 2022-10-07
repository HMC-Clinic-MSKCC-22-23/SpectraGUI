# import our libraries
import sys
import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc
import PyQt5.QtGui as qg

class InputPage(qw.QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = "Spectra input screen"
        self.left = 100
        self.top = 100
        self.width = 400
        self.height = 300
        self.initUI()
    
    def initUI(self):
        
        # make the window
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # make central widget
        self.wid = qw.QWidget(self)
        self.setCentralWidget(self.wid)

        # create layout grid
        layout = qw.QGridLayout()

        # title at the top of the screen
        title = qw.QLabel("SPECTRA")
        title.setFont(qg.QFont("Times", 20, 2))

        # anndata label
        annd = qw.QLabel("AnnData file path:")
        annd.setFont(qg.QFont("Times", 11))

        # anndata filepath box
        self.annd_box = qw.QLineEdit(self)
        self.annd_box.setFont(qg.QFont("Times", 11))

        # pathway annotations label
        pann = qw.QLabel("Pathway annotations:")
        pann.setFont(qg.QFont("Times", 11))

        # pathway annotations button
        pann_edit_button = qw.QPushButton("Edit")
        pann_edit_button.clicked.connect(self.pathway_edit_press)
        pann_edit_button.setFont(qg.QFont("Times", 11))

        # lambda label
        lam = qw.QLabel("Lambda value:")
        lam.setFont(qg.QFont("Times", 11))

        # lambda input box
        self.lam_box = qw.QLineEdit(self)
        self.lam_box.setText("0.0001")
        self.lam_box.setFont(qg.QFont("Times", 11))

        # advanced options label
        adv = qw.QLabel("Advanced options:")
        adv.setFont(qg.QFont("Times", 11))

        # advanced options button
        adv_edit_button = qw.QPushButton("Edit")
        adv_edit_button.clicked.connect(self.adv_edit_press)
        adv_edit_button.setFont(qg.QFont("Times", 11))

        # run button
        run_button = qw.QPushButton("Run")
        run_button.clicked.connect(self.run_button_press)
        run_button.setFont(qg.QFont("Times", 14))

        # add widgets to the layout
        layout.addWidget(title, 0, 0, 1, 5, qc.Qt.AlignCenter)
        layout.addWidget(annd, 1, 0, 1, 3)
        layout.addWidget(self.annd_box, 1, 3, 1, 2)
        layout.addWidget(pann, 2, 0, 1, 3)
        layout.addWidget(pann_edit_button, 2, 4)
        layout.addWidget(lam, 3, 0, 1, 3)
        layout.addWidget(self.lam_box, 3, 4)
        layout.addWidget(adv, 4, 0, 1, 3)
        layout.addWidget(adv_edit_button, 4, 4)
        layout.addWidget(run_button, 6, 0, 1, 5, qc.Qt.AlignCenter)

        layout.setRowStretch(5, 100)
        layout.setHorizontalSpacing(80)


        layout.setColumnStretch(0, 80)
        layout.setColumnStretch(1, 80)
        layout.setColumnStretch(2, 80)
        layout.setColumnStretch(3, 80)
        layout.setColumnStretch(4, 80)


        self.wid.setLayout(layout)

        self.show()


    def adv_edit_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Open up advanced options screen!")
        newWindow.exec()

    def pathway_edit_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Open up pathway annotations screen!")
        newWindow.exec()

    def run_button_press(self):
        newWindow = qw.QMessageBox(self.wid)
        if self.annd_box.text() == "":
            newWindow.setText("Please input a filepath for the Annotated Data")
        else: 
            newWindow.setText("Running in progress, with AnnData from \'" + self.annd_box.text() + "\' and a lambda value of " + self.lam_box.text())
        newWindow.exec()
    
if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = InputPage()
    sys.exit(app.exec_())


