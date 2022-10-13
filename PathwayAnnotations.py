# import our libraries
import sys
import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc
import PyQt5.QtGui as qg

class PathwayAnnotations(qw.QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Edit Pathway Annotations"
        self.left = 100
        self.top = 100
        self.width = 800
        self.height = 900
        self.initEditAnnotations()


    def initEditAnnotations(self):

        # make the window
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # make central widget
        self.wid = qw.QWidget(self)
        self.setCentralWidget(self.wid)

        # create layout grid
        layout = qw.QGridLayout()

        # title at the top of the screen
        title = qw.QLabel("Pathway Annotations")
        title.setFont(qg.QFont("Times", 15, 2))

        # gene-set annotations section
        gene_set_box = qw.QWidget()
        gene_set_layout = qw.QGridLayout()

        gene_set_label = qw.QLabel("Gene-Set Annotations")
        gene_set_label.setFont(qg.QFont("Times",12,2))

        gene_set_temp = qw.QLabel("EDIT STUFF HERE")
        gene_set_temp.setFrameShape(qw.QFrame.Box)

        delete_gene_set_button = qw.QPushButton("Delete")
        delete_gene_set_button.clicked.connect(self.delete_gene_set_press)
        delete_gene_set_button.setFont(qg.QFont("Times", 11))

        edit_gene_set_button = qw.QPushButton("Edit")
        edit_gene_set_button.clicked.connect(self.edit_gene_set_press)
        edit_gene_set_button.setFont(qg.QFont("Times", 11))

        add_gene_set_button = qw.QPushButton("Add")
        add_gene_set_button.clicked.connect(self.add_gene_set_press)
        add_gene_set_button.setFont(qg.QFont("Times", 11))

        # genes section
        gene_box = qw.QWidget()
        gene_box_layout = qw.QGridLayout()

        gene_box_label = qw.QLabel("Genes")
        gene_box_label.setFont(qg.QFont("Times",12,2))

        gene_box_temp = qw.QLabel("EDIT STUFF HERE")
        gene_box_temp.setFrameShape(qw.QFrame.Box)

        delete_gene_box_button = qw.QPushButton("Delete")
        delete_gene_box_button.clicked.connect(self.delete_gene_press)
        delete_gene_box_button.setFont(qg.QFont("Times", 11))

        edit_gene_box_button = qw.QPushButton("Edit")
        edit_gene_box_button.clicked.connect(self.edit_gene_press)
        edit_gene_box_button.setFont(qg.QFont("Times", 11))

        add_gene_box_button = qw.QPushButton("Add")
        add_gene_box_button.clicked.connect(self.add_gene_press)
        add_gene_box_button.setFont(qg.QFont("Times", 11))


        # save button
        save_button = qw.QPushButton("Save")
        save_button.clicked.connect(self.save_edit_press)
        save_button.setFont(qg.QFont("Times", 11))


        # ok button
        upload_button = qw.QPushButton("Upload")
        upload_button.clicked.connect(self.upload_edit_press)
        upload_button.setFont(qg.QFont("Times", 11))

        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8,qc.Qt.AlignCenter)
        

        # gene set annotations layout
        gene_set_layout.addWidget(gene_set_label,0,0,1,1,qc.Qt.AlignBottom)
        gene_set_layout.addWidget(gene_set_temp,1,0,3,8)
        gene_set_layout.addWidget(delete_gene_set_button,4,7,1,1,qc.Qt.AlignBottom)
        gene_set_layout.addWidget(edit_gene_set_button,4,6,1,1,qc.Qt.AlignBottom)
        gene_set_layout.addWidget(add_gene_set_button,4,5,1,1,qc.Qt.AlignBottom)
        gene_set_box.setLayout(gene_set_layout)

        # genes layout
        gene_box_layout.addWidget(gene_box_label,0,0,1,1,qc.Qt.AlignBottom)
        gene_box_layout.addWidget(gene_box_temp,1,0,3,8)
        gene_box_layout.addWidget(delete_gene_box_button,4,7,1,1,qc.Qt.AlignBottom)
        gene_box_layout.addWidget(edit_gene_box_button,4,6,1,1,qc.Qt.AlignBottom)
        gene_box_layout.addWidget(add_gene_box_button,4,5,1,1,qc.Qt.AlignBottom)
        gene_box.setLayout(gene_box_layout)


        layout.addWidget(gene_set_box,1,0,3,8)
        layout.addWidget(gene_box,4,0,3,8)

        
        layout.addWidget(save_button, 8, 7, 1, 1, qc.Qt.AlignBottom)
        layout.addWidget(upload_button, 8, 6, 1, 1, qc.Qt.AlignBottom)

        layout.setRowStretch(1, 20)
        layout.setRowStretch(4, 20)
        # layout.setVerticalSpacing(50)
    
        self.wid.setLayout(layout)

        # self.show()

    def save_edit_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Closes and saves pathway annotations edit screen")
        newWindow.exec()

    def delete_gene_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("deletes gene annotation")
        newWindow.exec()

    def delete_gene_set_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("deletes gene-set annotation")
        newWindow.exec()

    def upload_edit_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("new window for uploading annotations")
        newWindow.exec()

    def edit_gene_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Edits gene annotation")
        newWindow.exec()
    
    def edit_gene_set_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Edits gene-set annotation")
        newWindow.exec()

    def add_gene_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Adds gene annotation")
        newWindow.exec()

    def add_gene_set_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Adds gene-set annotation")
        newWindow.exec()


if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = PathwayAnnotations()
    sys.exit(app.exec_())