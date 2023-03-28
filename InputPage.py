# import our libraries
import sys
import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc
import PyQt5.QtGui as qg
import PathwayAnnotations as pa
import AdvancedOptions as ao

from os.path import exists

class InputPage(qw.QMainWindow):

    def __init__(self, screen_width, screen_height):
        super().__init__()
        self.title = "Spectra input screen"
        self.left = 100
        self.top = 100
        self.width = int(screen_width // 4)
        self.height = int(screen_height // 4)
        self.path_ann = pa.PathwayAnnotations(self.width, self.height)
        self.adv_op = ao.AdvancedOptions()
        self.initUI()

    def closeEvent(self, a0: qg.QCloseEvent) -> None:
        self.path_ann.close()
        self.adv_op.close()
        return super().closeEvent(a0)
    
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
        ctype = qw.QLabel("AnnData key for cell type:")
        ctype.setFont(qg.QFont("Times", 11))

        # lambda input box
        self.ctype_box = qw.QLineEdit(self)
        self.ctype_box.setText("cell_type")
        self.ctype_box.setFont(qg.QFont("Times", 11))

        # lambda label
        lam = qw.QLabel("Lambda value:")
        lam.setFont(qg.QFont("Times", 11))

        # lambda input box
        self.lam_box = qw.QLineEdit(self)
        self.lam_box.setText("0.0001")
        self.lam_box.setFont(qg.QFont("Times", 11))

        # epochs label
        epochs = qw.QLabel("Number of epochs:")
        epochs.setFont(qg.QFont("Times", 11))

        # epochs input box
        self.epochs_box = qw.QLineEdit(self)
        self.epochs_box.setText("10000")
        self.epochs_box.setFont(qg.QFont("Times", 11))

        # advanced options label
        adv = qw.QLabel("Advanced options:")
        adv.setFont(qg.QFont("Times", 11))

        # advanced options button
        adv_edit_button = qw.QPushButton("Edit")
        adv_edit_button.clicked.connect(self.adv_edit_press)
        adv_edit_button.setFont(qg.QFont("Times", 11))

        # run button
        self.run_button = qw.QPushButton("Run")
        self.run_button.setFont(qg.QFont("Times", 14))

        # add widgets to the layout
        layout.addWidget(title, 0, 0, 1, 5, qc.Qt.AlignCenter)
        layout.addWidget(annd, 1, 0, 1, 3)
        layout.addWidget(self.annd_box, 1, 2, 1, 3)
        layout.addWidget(ctype, 2, 0, 1, 3)
        layout.addWidget(self.ctype_box, 2, 3, 1, 2)
        layout.addWidget(pann, 3, 0, 1, 3)
        layout.addWidget(pann_edit_button, 3, 4)
        layout.addWidget(lam, 4, 0, 1, 3)
        layout.addWidget(self.lam_box, 4, 4)
        layout.addWidget(epochs, 5, 0, 1, 3)
        layout.addWidget(self.epochs_box, 5, 4)
        layout.addWidget(adv, 6, 0, 1, 3)
        layout.addWidget(adv_edit_button, 6, 4)
        layout.addWidget(self.run_button, 8, 0, 1, 5, qc.Qt.AlignCenter)

        layout.setRowStretch(7, self.height // 3)
        layout.setHorizontalSpacing(self.height // 4)


        layout.setColumnStretch(0, self.width // 3)
        layout.setColumnStretch(1, self.width // 8)
        layout.setColumnStretch(2, self.width // 8)
        layout.setColumnStretch(3, self.width // 5)
        layout.setColumnStretch(4, self.width // 5)


        self.wid.setLayout(layout)

        self.show()


    def adv_edit_press(self):
        self.adv_op.show()

    def pathway_edit_press(self):
        self.path_ann.show()

    def run_button_press(self):
        newWindow = qw.QMessageBox(self.wid)
        if self.annd_box.text() == "":
            newWindow.setText("Please input a filepath for the Annotated Data")
            newWindow.exec()
            return False
        elif not exists(self.annd_box.text()):
            newWindow.setText("The inputted filepath for the Annotated Data is invalid")
            newWindow.exec()
            return False
        else: 
            return True
        
    
if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = InputPage(app.primaryScreen().size().width(), app.primaryScreen().size().height())
    sys.exit(app.exec_())


