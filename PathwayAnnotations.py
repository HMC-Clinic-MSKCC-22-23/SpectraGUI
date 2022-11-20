# import our libraries
import sys
import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc
import PyQt5.QtGui as qg
import ast, os, csv, json


class PathwayAnnotations(qw.QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Edit Pathway Annotations"
        self.left = 150
        self.top = 150
        self.width = 800
        self.height = 900
        self.genes_dict = {}
        self.edit_cell_type = editCellTypeAnnotationWindow()
        self.edit_gene = editGeneAnnotationWindow()
        self.new_gene = newGeneAnnotationWindow()
        self.new_cell_type = newCellTypeAnnotationWindow()
        self.current_cell_type = None
        self.current_pathway = None
        self.initEditAnnotations()

    def closeEvent(self, a0: qg.QCloseEvent) -> None:
        self.edit_cell_type.close()
        self.edit_gene.close()
        self.new_gene.close()
        self.new_cell_type.close()
        return super().closeEvent(a0)

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

        gene_set_label = qw.QLabel("Cell-Type Annotations")
        gene_set_label.setFont(qg.QFont("Times", 12, 2))

        self.gene_set_table = qw.QTableWidget()
        self.gene_set_table.setColumnCount(2)
        self.gene_set_table.setHorizontalHeaderLabels(["Cell-Types", "Pathway Names"])
        header = self.gene_set_table.horizontalHeader()
        header.setSectionResizeMode(0, qw.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, qw.QHeaderView.Stretch)

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
        gene_box_label.setFont(qg.QFont("Times", 12, 2))

        self.gene_table = qw.QTableWidget()
        self.gene_table.setColumnCount(2)
        self.gene_table.setHorizontalHeaderLabels(["Pathway Names", "Genes"])
        header = self.gene_table.horizontalHeader()
        header.setSectionResizeMode(0, qw.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, qw.QHeaderView.Stretch)

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
        upload_button.clicked.connect(self.csv_press)
        upload_button.setFont(qg.QFont("Times", 11))

        # cancel button
        cancel_button = qw.QPushButton("Cancel")
        cancel_button.clicked.connect(self.cancel_press)
        cancel_button.setFont(qg.QFont("Times", 11))

        # setting up the table connections
        self.cell_type_table.itemDoubleClicked.connect(self.propogate_gene_table)
        self.cell_type_table.itemClicked.connect(self.update_current_cell_type)
        self.gene_table.itemClicked.connect(self.update_current_pathway)

        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8, qc.Qt.AlignCenter)

        # gene set annotations layout
        cell_type_layout.addWidget(cell_type_label, 0, 0, 1, 1, qc.Qt.AlignBottom)
        cell_type_layout.addWidget(self.cell_type_table, 1, 0, 3, 8)
        cell_type_layout.addWidget(delete_cell_type_button, 4, 7, 1, 1, qc.Qt.AlignBottom)
        cell_type_layout.addWidget(edit_cell_type_button, 4, 6, 1, 1, qc.Qt.AlignBottom)
        cell_type_layout.addWidget(add_cell_type_button, 4, 5, 1, 1, qc.Qt.AlignBottom)
        cell_type_box.setLayout(cell_type_layout)

        # genes layout
        gene_box_layout.addWidget(gene_box_label, 0, 0, 1, 1, qc.Qt.AlignBottom)
        gene_box_layout.addWidget(self.gene_table, 1, 0, 3, 8)
        gene_box_layout.addWidget(delete_gene_box_button, 4, 7, 1, 1, qc.Qt.AlignBottom)
        gene_box_layout.addWidget(edit_gene_box_button, 4, 6, 1, 1, qc.Qt.AlignBottom)
        gene_box_layout.addWidget(add_gene_box_button, 4, 5, 1, 1, qc.Qt.AlignBottom)
        gene_box.setLayout(gene_box_layout)

        layout.addWidget(gene_set_box, 1, 0, 3, 8)
        layout.addWidget(gene_box, 4, 0, 3, 8)

        layout.addWidget(cancel_button, 8, 7, 1, 1, qc.Qt.AlignBottom)
        layout.addWidget(save_button, 8, 6, 1, 1, qc.Qt.AlignBottom)
        layout.addWidget(upload_button, 8, 5, 1, 1, qc.Qt.AlignBottom)

        layout.setRowStretch(1, 20)
        layout.setRowStretch(4, 20)
        # layout.setVerticalSpacing(50)

        self.wid.setLayout(layout)

        # self.show()

    def csv_press(self,fileName):
        try:
            fileName, _ = qw.QFileDialog.getOpenFileName(self, "Open CSV",
                                                           (qc.QDir.homePath()), "CSV (*.csv *.tsv)")
            ff = open(fileName, 'r')
            mytext = ff.read()
        #            print(mytext)
            ff.close()
        # self.genes.insertPlainText(fileName)
            if fileName:
                f = open(fileName, 'r')
                with f:
                    self.fname = os.path.splitext(str(fileName))[0].split("/")[-1]
                    if mytext.count(';') <= mytext.count(','):  # tab?
                        reader = csv.reader(f, delimiter=',')
                        for row in reader:
                            items = [field for field in row]
                                if items[0] == "Cell_Type":
                                    continue
                        # self.genes.appendRow(items)

                            i = 2
                            if items[0] not in self.genes_dict:
                                self.genes_dict[items[0]] = {}
                            if items[1] not in self.genes_dict[items[0]]:
                                self.genes_dict[items[0]][items[1]] = []
                            while len(items) > i:
                                if items[i] != "":
                                    self.genes_dict[items[0]][items[1]] += [items[i]]
                                i += 1

                    # self.tableView.resizeColumnsToContents()
                    else:
                        reader = csv.reader(f, delimiter='\t')
                        for row in reader:
                            items = [field for field in row]
                                if items[0] == "Cell_Type":
                                    continue
                        # self.genes.appendRow(items)

                            i = 2
                            if items[0] not in self.genes_dict:
                                self.genes_dict[items[0]] = {}
                            if items[1] not in self.genes_dict[items[0]]:
                                self.genes_dict[items[0]][items[1]] = []
                            while len(items) > i:
                                if items[i] != "":
                                    self.genes_dict[items[0]][items[1]] += [items[i]]
                                i += 1
            self.updateTable(self.genes_dict, self.cell_type_table)
        except:
            newWindow = qw.QMessageBox(self.wid)
            newWindow.setText("Error in uploading file - try again")
            newWindow.exec()

    def save_edit_press(self):
        self.hide()

    def delete_gene_press(self):
        current_dict = self.genes_dict.get(self.current_cell_type)
        reply = qw.QMessageBox.question(self.wid, "Delete",
                                        "Are you sure you want to delete the pathway " + self.current_pathway + "?",
                                        qw.QMessageBox.Yes | qw.QMessageBox.No, qw.QMessageBox.No)
        if reply == qw.QMessageBox.Yes:
            current_dict.pop(self.current_pathway)
            self.updateTable(current_dict, self.gene_table)
            self.updateTable(self.genes_dict, self.cell_type_table)

    def delete_gene_set_press(self):
        reply = qw.QMessageBox.question(self.wid, "Delete",
                                        "Are you sure you want to delete the cell_type " + self.current_cell_type + "?",
                                        qw.QMessageBox.Yes | qw.QMessageBox.No, qw.QMessageBox.No)
        if reply == qw.QMessageBox.Yes:
            self.genes_dict.pop(self.current_cell_type)
            self.updateTable(self.genes_dict, self.cell_type_table)


    def edit_gene_press(self):
        self.edit_gene.show()
        self.edit_gene.ok_button.clicked.connect(self.save_gene_edit)

    def save_gene_edit(self):
        try:
            new_gene = {self.edit_gene.pathway_box.text(): self.edit_gene.gene_box.text()}
            current_dict = self.genes_dict.get(self.current_cell_type)
            current_dict = current_dict | new_gene
            self.genes_dict.update({self.current_cell_type: current_dict})
            self.updateTable(current_dict, self.gene_table)
            self.updateTable(self.genes_dict, self.cell_type_table)
            self.edit_gene.hide()
        except:
            newWindow = qw.QMessageBox(self.wid)
            newWindow.setText("Error saving gene annotation - try again")
            newWindow.exec()

    def edit_gene_set_press(self):
        if self.current_cell_type is not None:
            self.edit_cell_type.gs_box.setText(self.current_cell_type)
            self.edit_cell_type.pathway_box.setText(str(self.genes_dict.get(self.current_cell_type)))
        self.edit_cell_type.show()
        self.edit_cell_type.ok_button.clicked.connect(self.save_cell_type_edit)

    def save_gene_set_edit(self):
        try:
            new_pathway = {self.edit_cell_type.gs_box.text(): ast.literal_eval(self.edit_cell_type.pathway_box.text())}
            self.current_cell_type = self.edit_cell_type.gs_box.text()
            self.genes_dict.update(new_pathway)
            self.updateTable(new_pathway.get(self.current_cell_type), self.gene_table)
            self.updateTable(self.genes_dict, self.cell_type_table)
            self.edit_cell_type.hide()
        except:
            newWindow = qw.QMessageBox(self.wid)
            newWindow.setText("Error saving cell type - try again")
            newWindow.exec()

    def add_gene_press(self):
        self.new_gene.gene_box.setText("")
        self.new_gene_set.pathway_box.setText("")
        self.new_gene.show()
        self.new_gene.ok_button.clicked.connect(self.save_new_gene)

    def save_new_gene(self):
        try:
            new_gene = {self.new_gene.pathway_box.text(): self.new_gene.gene_box.text()}
            current_dict = self.genes_dict.get(self.current_cell_type)
            current_dict = current_dict | new_gene
            self.genes_dict.update({self.current_cell_type: current_dict})
            self.updateTable(current_dict, self.gene_table)
            self.updateTable(self.genes_dict, self.cell_type_table)
            self.new_gene.hide()
        except:
            newWindow = qw.QMessageBox(self.wid)
            newWindow.setText("Error adding new gene annotation - try again")
            newWindow.exec()

    def add_gene_set_press(self):
        self.new_cell_type.gs_box.setText("")
        self.new_cell_type.pathway_box.setText("")
        self.new_cell_type.show()
        self.new_cell_type.ok_button.clicked.connect(self.save_new_cell_type)

    def save_new_gene_set(self):
        try:
            new_gs = {self.new_cell_type.gs_box.text(): ast.literal_eval(self.new_cell_type.pathway_box.text())}
            self.genes_dict = self.genes_dict | new_gs
            self.updateTable(self.genes_dict, self.cell_type_table)
            self.new_gene_set.hide()
        except:
            newWindow = qw.QMessageBox(self.wid)
            newWindow.setText("Error adding new cell type - try again")
            newWindow.exec()

    def cancel_press(self):
        self.close()

    def updateTable(self, dict, table):
        table.setRowCount(len(dict.keys()))

        for row in range(0, len(dict.keys())):
            keyVal = qw.QTableWidgetItem(str(list(dict.keys())[row]))
            valueVal = qw.QTableWidgetItem(str(dict.get(list(dict.keys())[row])))
            table.setItem(row, 0, keyVal)
            table.setItem(row, 1, valueVal)

    def propogate_gene_table(self):
        self.updateTable(self.genes_dict.get(self.cell_type_table.item(self.cell_type_table.currentRow(), 0).text()),
                         self.gene_table)
        self.current_cell_type = self.cell_type_table.item(self.cell_type_table.currentRow(), 0).text()

    def update_current_cell_type(self):
        self.current_cell_type = self.cell_type_table.item(self.cell_type_table.currentRow(), 0).text()

    def update_current_pathway(self):
        self.current_pathway = self.gene_table.item(self.gene_table.currentRow(), 0).text()


class editCellTypeAnnotationWindow(qw.QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Edit Annotation"
        self.left = 100
        self.top = 100
        self.width = 900
        self.height = 200
        self.editCellTypeAnnotation()

    def editCellTypeAnnotation(self):
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
        ct_name = qw.QLabel("Cell-Type Name:")
        ct_name.setFont(qg.QFont("Times", 11))

        # gene set box
        self.ct_box = qw.QLineEdit(self)
        self.ct_box.setFont(qg.QFont("Times", 11))

        # pathway name label
        pathway_name = qw.QLabel("Pathway Annotation Name:")
        pathway_name.setFont(qg.QFont("Times", 11))

        # pathway name box
        self.pathway_box = qw.QLineEdit(self)
        self.pathway_box.setFont(qg.QFont("Times", 11))

        # ok button
        self.ok_button = qw.QPushButton("OK")
        self.ok_button.setFont(qg.QFont("Times", 11))

        # cancel button
        cancel_button = qw.QPushButton("Cancel")
        cancel_button.clicked.connect(self.cancel_press)
        cancel_button.setFont(qg.QFont("Times", 11))

        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8, qc.Qt.AlignCenter)
        layout.addWidget(gs_name, 1, 0, 1, 1)
        layout.addWidget(self.gs_box, 1, 1, 1, 2)
        layout.addWidget(pathway_name, 2, 0, 1, 1)
        layout.addWidget(self.pathway_box, 2, 1, 1, 2)

        layout.addWidget(cancel_button, 3, 2, 1, 1, qc.Qt.AlignRight)
        layout.addWidget(self.ok_button, 3, 1, 1, 1, qc.Qt.AlignRight)

        layout.setColumnStretch(1, 80)
        self.wid.setLayout(layout)

    def cancel_press(self):
        self.close()


class newCellTypeAnnotationWindow(qw.QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Add New Annotation"
        self.left = 100
        self.top = 100
        self.width = 900
        self.height = 200
        self.newCellTypeAnnotation()

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
        ct_name = qw.QLabel("Cell-Type Name:")
        ct_name.setFont(qg.QFont("Times", 11))

        # gene set box
        self.ct_box = qw.QLineEdit(self)
        self.ct_box.setFont(qg.QFont("Times", 11))

        # pathway name label
        pathway_name = qw.QLabel("Pathway Annotation Name:")
        pathway_name.setFont(qg.QFont("Times", 11))

        # pathway name box
        self.pathway_box = qw.QLineEdit(self)
        self.pathway_box.setFont(qg.QFont("Times", 11))

        # ok button
        self.ok_button = qw.QPushButton("OK")
        self.ok_button.setFont(qg.QFont("Times", 11))

        # cancel button
        cancel_button = qw.QPushButton("Cancel")
        cancel_button.clicked.connect(self.cancel_press)
        cancel_button.setFont(qg.QFont("Times", 11))

        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8, qc.Qt.AlignCenter)
        layout.addWidget(gs_name, 1, 0, 1, 1)
        layout.addWidget(self.gs_box, 1, 1, 1, 2)
        layout.addWidget(pathway_name, 2, 0, 1, 1)
        layout.addWidget(self.pathway_box, 2, 1, 1, 2)

        layout.addWidget(cancel_button, 3, 2, 1, 1, qc.Qt.AlignRight)
        layout.addWidget(self.ok_button, 3, 1, 1, 1, qc.Qt.AlignRight)

        layout.setColumnStretch(1, 80)
        self.wid.setLayout(layout)

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

        # pathway name label
        pathway_name = qw.QLabel("Pathway Annotation Name:")
        pathway_name.setFont(qg.QFont("Times", 11))

        # pathway name box
        self.pathway_box = qw.QLineEdit(self)
        self.pathway_box.setFont(qg.QFont("Times", 11))

        # gene label
        gene_name = qw.QLabel("Gene Name:")
        gene_name.setFont(qg.QFont("Times", 11))

        # gene set box
        self.gene_box = qw.QLineEdit(self)
        self.gene_box.setFont(qg.QFont("Times", 11))

        # ok button
        self.ok_button = qw.QPushButton("OK")
        self.ok_button.setFont(qg.QFont("Times", 11))

        # cancel button
        cancel_button = qw.QPushButton("Cancel")
        cancel_button.clicked.connect(self.cancel_press)
        cancel_button.setFont(qg.QFont("Times", 11))

        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8, qc.Qt.AlignCenter)
        layout.addWidget(pathway_name, 1, 0, 1, 1)
        layout.addWidget(self.pathway_box, 1, 1, 1, 2)
        layout.addWidget(gene_name, 2, 0, 1, 1)
        layout.addWidget(self.gene_box, 2, 1, 1, 2)

        layout.addWidget(cancel_button, 3, 2, 1, 1, qc.Qt.AlignRight)
        layout.addWidget(self.ok_button, 3, 1, 1, 1, qc.Qt.AlignRight)

        layout.setColumnStretch(1, 80)
        self.wid.setLayout(layout)

        self.wid.setLayout(layout)

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
        self.gene_box = qw.QLineEdit(self)
        self.gene_box.setFont(qg.QFont("Times", 11))

        # pathway name label
        pathway_name = qw.QLabel("Pathway Annotation Name:")
        pathway_name.setFont(qg.QFont("Times", 11))

        # pathway name box
        self.pathway_box = qw.QLineEdit(self)
        self.pathway_box.setFont(qg.QFont("Times", 11))

        # ok button
        self.ok_button = qw.QPushButton("OK")
        self.ok_button.setFont(qg.QFont("Times", 11))

        # cancel button
        cancel_button = qw.QPushButton("Cancel")
        cancel_button.clicked.connect(self.cancel_press)
        cancel_button.setFont(qg.QFont("Times", 11))

        # add widgets to layout(s)
        layout.addWidget(title, 0, 0, 1, 8, qc.Qt.AlignCenter)
        layout.addWidget(pathway_name, 1, 0, 1, 1)
        layout.addWidget(self.pathway_box, 1, 1, 1, 2)
        layout.addWidget(gene_name, 2, 0, 1, 1)
        layout.addWidget(self.gene_box, 2, 1, 1, 2)

        layout.addWidget(cancel_button, 3, 2, 1, 1, qc.Qt.AlignRight)
        layout.addWidget(self.ok_button, 3, 1, 1, 1, qc.Qt.AlignRight)

        layout.setColumnStretch(1, 80)
        self.wid.setLayout(layout)

        self.wid.setLayout(layout)

    def cancel_press(self):
        self.close()


if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = PathwayAnnotations()
    sys.exit(app.exec_())
