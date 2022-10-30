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
        self.genes_dict = {} # should this get passed as a prop??
        self.edit_gene_set = editGeneSetAnnotationWindow()
        self.edit_gene = editGeneAnnotationWindow()
        self.new_gene = newGeneAnnotationWindow()
        self.new_gene_set = newGeneSetAnnotationWindow()
        self.initEditAnnotations()

    # probably needs a close event

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

        gene_set_table = qw.QTableWidget()
        gene_set_table.setColumnCount(2)
        gene_set_table.setHorizontalHeaderLabels(["Gene-Sets", "Pathway Names"])
        # do we need to update table on init? yes if it is a prop

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

        gene_table = qw.QTableWidget()
        gene_table.setColumnCount(2)
        gene_table.setHorizontalHeaderLabels(["Pathway Names", "Genes"])

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

        # cancel button
        cancel_button = qw.QPushButton("Cancel")
        cancel_button.clicked.connect(self.cancel_press)
        cancel_button.setFont(qg.QFont("Times", 11))

        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8,qc.Qt.AlignCenter)
        

        # gene set annotations layout
        gene_set_layout.addWidget(gene_set_label,0,0,1,1,qc.Qt.AlignBottom)
        gene_set_layout.addWidget(gene_set_table,1,0,3,8)
        gene_set_layout.addWidget(delete_gene_set_button,4,7,1,1,qc.Qt.AlignBottom)
        gene_set_layout.addWidget(edit_gene_set_button,4,6,1,1,qc.Qt.AlignBottom)
        gene_set_layout.addWidget(add_gene_set_button,4,5,1,1,qc.Qt.AlignBottom)
        gene_set_box.setLayout(gene_set_layout)

        # genes layout
        gene_box_layout.addWidget(gene_box_label,0,0,1,1,qc.Qt.AlignBottom)
        gene_box_layout.addWidget(gene_table,1,0,3,8)
        gene_box_layout.addWidget(delete_gene_box_button,4,7,1,1,qc.Qt.AlignBottom)
        gene_box_layout.addWidget(edit_gene_box_button,4,6,1,1,qc.Qt.AlignBottom)
        gene_box_layout.addWidget(add_gene_box_button,4,5,1,1,qc.Qt.AlignBottom)
        gene_box.setLayout(gene_box_layout)


        layout.addWidget(gene_set_box,1,0,3,8)
        layout.addWidget(gene_box,4,0,3,8)

        layout.addWidget(cancel_button, 8, 7, 1, 1, qc.Qt.AlignBottom)
        layout.addWidget(save_button, 8, 6, 1, 1, qc.Qt.AlignBottom)
        layout.addWidget(upload_button, 8, 5, 1, 1, qc.Qt.AlignBottom)

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
        self.edit_gene.show()
    
    def edit_gene_set_press(self):
        self.edit_gene_set.show()

    def add_gene_press(self):
        self.new_gene.show()

    def add_gene_set_press(self):
        self.new_gene_set.show()
    
    def cancel_press(self):
        self.close()

    def updateTable(self, dict, table):
        table.setRowCount(len(dict.keys()))

        for row in range(0,len(dict.keys())):
            keyVal = qw.QTableWidgetItem(str(dict.keys()[row]))
            valueVal = qw.QTableWidgetItem(str(dict.get(dict.keys()[row])))
            table.setItem(row, 0, keyVal)
            table.setItem(row, 1, valueVal)

        

class editGeneSetAnnotationWindow(qw.QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Edit Annotation"
        self.left = 100
        self.top = 100
        self.width = 900
        self.height = 200
        self.editGeneSetAnnotation()


    def editGeneSetAnnotation(self):
        # make the window
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # make central widget
        self.wid = qw.QWidget(self)
        self.setCentralWidget(self.wid)

        # create layout grid
        layout = qw.QGridLayout()

        # title at the top of the screen
        title = qw.QLabel("Edit Annotation")
        title.setFont(qg.QFont("Times", 12, 2))

        # gene set label
        gs_name = qw.QLabel("Gene-Set Name:")
        gs_name.setFont(qg.QFont("Times", 11))

        # gene set box
        gs_box = qw.QLineEdit(self)
        gs_box.setFont(qg.QFont("Times", 11))

        # pathway name label
        pathway_name = qw.QLabel("Pathway Annotation Name:")
        pathway_name.setFont(qg.QFont("Times", 11))

        # pathway name box
        pathway_box = qw.QLineEdit(self)
        pathway_box.setFont(qg.QFont("Times", 11))

        # ok button
        upload_button = qw.QPushButton("OK")
        upload_button.clicked.connect(self.ok_press)
        upload_button.setFont(qg.QFont("Times", 11))

        # cancel button
        cancel_button = qw.QPushButton("Cancel")
        cancel_button.clicked.connect(self.cancel_press)
        cancel_button.setFont(qg.QFont("Times", 11))


        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8,qc.Qt.AlignCenter)
        layout.addWidget(gs_name, 1, 0, 1, 1)
        layout.addWidget(gs_box, 1, 1, 1, 2)
        layout.addWidget(pathway_name,2, 0, 1, 1)
        layout.addWidget(pathway_box, 2, 1, 1, 2)

        layout.addWidget(cancel_button, 3, 2, 1, 1, qc.Qt.AlignRight)
        layout.addWidget(upload_button, 3, 1, 1, 1, qc.Qt.AlignRight)


        layout.setColumnStretch(1, 80)
        self.wid.setLayout(layout)


    def ok_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Closes and updates table")
        newWindow.exec()

    def cancel_press(self):
        self.close()

class newGeneSetAnnotationWindow(qw.QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Add New Annotation"
        self.left = 100
        self.top = 100
        self.width = 900
        self.height = 200
        self.newGeneSetAnnotation()


    def newGeneSetAnnotation(self):
        # make the window
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # make central widget
        self.wid = qw.QWidget(self)
        self.setCentralWidget(self.wid)

        # create layout grid
        layout = qw.QGridLayout()

        # title at the top of the screen
        title = qw.QLabel("New Annotation")
        title.setFont(qg.QFont("Times", 12, 2))

        # gene set label
        gs_name = qw.QLabel("Gene-Set Name:")
        gs_name.setFont(qg.QFont("Times", 11))

        # gene set box
        gs_box = qw.QLineEdit(self)
        gs_box.setFont(qg.QFont("Times", 11))

        # pathway name label
        pathway_name = qw.QLabel("Pathway Annotation Name:")
        pathway_name.setFont(qg.QFont("Times", 11))

        # pathway name box
        pathway_box = qw.QLineEdit(self)
        pathway_box.setFont(qg.QFont("Times", 11))

        # ok button
        upload_button = qw.QPushButton("OK")
        upload_button.clicked.connect(self.ok_press)
        upload_button.setFont(qg.QFont("Times", 11))

        # cancel button
        cancel_button = qw.QPushButton("Cancel")
        cancel_button.clicked.connect(self.cancel_press)
        cancel_button.setFont(qg.QFont("Times", 11))


        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8,qc.Qt.AlignCenter)
        layout.addWidget(gs_name, 1, 0, 1, 1)
        layout.addWidget(gs_box, 1, 1, 1, 2)
        layout.addWidget(pathway_name,2, 0, 1, 1)
        layout.addWidget(pathway_box, 2, 1, 1, 2)

        layout.addWidget(cancel_button, 3, 2, 1, 1, qc.Qt.AlignRight)
        layout.addWidget(upload_button, 3, 1, 1, 1, qc.Qt.AlignRight)


        layout.setColumnStretch(1, 80)
        self.wid.setLayout(layout)


    def ok_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Closes and updates table")
        newWindow.exec()

    def cancel_press(self):
        self.close()


class newGeneAnnotationWindow(qw.QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Add New Annotation"
        self.left = 100
        self.top = 100
        self.width = 900
        self.height = 200
        self.newAnnotation()


    def newAnnotation(self):
        #  make the window
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # make central widget
        self.wid = qw.QWidget(self)
        self.setCentralWidget(self.wid)

        # create layout grid
        layout = qw.QGridLayout()

        # title at the top of the screen
        title = qw.QLabel("New Annotation")
        title.setFont(qg.QFont("Times", 12, 2))

        # gene label
        gene_name = qw.QLabel("Gene Name:")
        gene_name.setFont(qg.QFont("Times", 11))

        # gene set box
        gene_box = qw.QLineEdit(self)
        gene_box.setFont(qg.QFont("Times", 11))

        # pathway name label
        pathway_name = qw.QLabel("Pathway Annotation Name:")
        pathway_name.setFont(qg.QFont("Times", 11))

        # pathway name box
        pathway_box = qw.QLineEdit(self)
        pathway_box.setFont(qg.QFont("Times", 11))

        # ok button
        upload_button = qw.QPushButton("OK")
        upload_button.clicked.connect(self.ok_press)
        upload_button.setFont(qg.QFont("Times", 11))

        # cancel button
        cancel_button = qw.QPushButton("Cancel")
        cancel_button.clicked.connect(self.cancel_press)
        cancel_button.setFont(qg.QFont("Times", 11))


        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8,qc.Qt.AlignCenter)
        layout.addWidget(pathway_name,1, 0, 1, 1)
        layout.addWidget(pathway_box, 1, 1, 1, 2)
        layout.addWidget(gene_name, 2, 0, 1, 1)
        layout.addWidget(gene_box, 2, 1, 1, 2)

        layout.addWidget(cancel_button, 3, 2, 1, 1, qc.Qt.AlignRight)
        layout.addWidget(upload_button, 3, 1, 1, 1, qc.Qt.AlignRight)
    

        layout.setColumnStretch(1, 80)
        self.wid.setLayout(layout)


        self.wid.setLayout(layout)

    def ok_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Closes and updates table")
        newWindow.exec()

    def cancel_press(self):
        self.close()


class editGeneAnnotationWindow(qw.QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Edit Annotation"
        self.left = 100
        self.top = 100
        self.width = 900
        self.height = 200
        self.newAnnotation()


    def newAnnotation(self):
        #  make the window
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # make central widget
        self.wid = qw.QWidget(self)
        self.setCentralWidget(self.wid)

        # create layout grid
        layout = qw.QGridLayout()

        # title at the top of the screen
        title = qw.QLabel("Edit Annotation")
        title.setFont(qg.QFont("Times", 12, 2))

        # gene label
        gene_name = qw.QLabel("Gene Name:")
        gene_name.setFont(qg.QFont("Times", 11))

        # gene set box
        gene_box = qw.QLineEdit(self)
        gene_box.setFont(qg.QFont("Times", 11))

        # pathway name label
        pathway_name = qw.QLabel("Pathway Annotation Name:")
        pathway_name.setFont(qg.QFont("Times", 11))

        # pathway name box
        pathway_box = qw.QLineEdit(self)
        pathway_box.setFont(qg.QFont("Times", 11))

        # ok button
        upload_button = qw.QPushButton("OK")
        upload_button.clicked.connect(self.ok_press)
        upload_button.setFont(qg.QFont("Times", 11))

        # cancel button
        cancel_button = qw.QPushButton("Cancel")
        cancel_button.clicked.connect(self.cancel_press)
        cancel_button.setFont(qg.QFont("Times", 11))


        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8,qc.Qt.AlignCenter)
        layout.addWidget(pathway_name,1, 0, 1, 1)
        layout.addWidget(pathway_box, 1, 1, 1, 2)
        layout.addWidget(gene_name, 2, 0, 1, 1)
        layout.addWidget(gene_box, 2, 1, 1, 2)

        layout.addWidget(cancel_button, 3, 2, 1, 1, qc.Qt.AlignRight)
        layout.addWidget(upload_button, 3, 1, 1, 1, qc.Qt.AlignRight)
    

        layout.setColumnStretch(1, 80)
        self.wid.setLayout(layout)


        self.wid.setLayout(layout)

    def ok_press(self):
        newWindow = qw.QMessageBox(self.wid)
        newWindow.setText("Closes and updates table")
        newWindow.exec()

    def cancel_press(self):
        self.close()


if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = PathwayAnnotations()
    sys.exit(app.exec_())