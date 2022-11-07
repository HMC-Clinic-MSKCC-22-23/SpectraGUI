# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'OutputUI.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.

import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
matplotlib.use("Qt5Agg")
from spectra import spectra as spc

from PyQt5 import QtCore, QtGui, QtWidgets
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

class OutputPage2(object):

    # def __init__(self, anndata = None, gene_dict = None, cell_type_key = None, lambda_val = None, highly_var = None, rho_val = None, delta_val = None, kappa_val = None, use_weights = None, top_genes = None):
    def __init__(self, screen_width, screen_height, anndata = None, model = None):
        sc.set_figure_params(facecolor="F0F0F0")

        # self.gene_dictionary = gene_dict
        # self.cell_type_key = cell_type_key
        # self.lam = lambda_val
        # self.highly_var = highly_var
        # self.rho = rho_val
        # self.delta = delta_val
        # self.kappa = kappa_val
        # self.use_weights = use_weights
        # self.top_genes = top_genes

        self.width = int(screen_width // 2.2)
        self.height = int(screen_height // 2.2)
        

        self.curr_factor = None

        self.anndata = anndata

        self.MainWindow = QtWidgets.QWidget()
        self.setupUi()
    
    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:
        return super().closeEvent(a0)

    def test_umap(self):
        if self.curr_factor is not None:
            self.curr_factor += 1
        else:
            self.curr_factor = 0

        if self.anndata:
            self.recolor_umap()

    def draw_umap(self):

        self.ax = self.canvas.figure.subplots()
        self.ax.grid(False)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_title("UMAP")
        # if we have loaded anndata, draw colorless umap
        if self.anndata:
            self.umap_plot = self.ax.scatter(self.anndata.obsm["X_umap"][:,0], self.anndata.obsm["X_umap"][:,1], color = "grey", s = 1, marker = ",")
        else: # just draw a placeholder
            self.umap_plot = self.ax.plot([2,5,4,3.5,4,5,2])
    
    def recolor_umap(self):

        if self.umap_plot:
            self.umap_plot.remove()
            
        self.umap_plot = self.ax.scatter(self.anndata.obsm["X_umap"][:,0], self.anndata.obsm["X_umap"][:,1], c = self.anndata.obsm["SPECTRA_cell_scores"][:,self.curr_factor], s = 1.1, marker = ",", cmap = "inferno")
        self.canvas.draw()
        self.canvas.print_figure("UMAP_plot.png")

    def setupUi(self):

        self.MainWindow.setObjectName("MainWindow")
        self.MainWindow.resize(self.width, self.height)


        self.umap = QtWidgets.QFrame(self.MainWindow)
        self.umap.setGeometry(QtCore.QRect(0, 0, self.width // 2, int(self.height // 1.5)))
        self.umap.setObjectName("umapFrame")

        self.umap_box = QtWidgets.QVBoxLayout()
        self.canvas = FigureCanvasQTAgg(plt.Figure())
        self.umap_box.addWidget(self.canvas)

        self.umap.setLayout(self.umap_box)

        self.draw_umap()


        self.heatMap = QtWidgets.QFrame(self.MainWindow)
        self.heatMap.setGeometry(QtCore.QRect(350, 30, 281, 231))
        self.heatMap.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.heatMap.setFrameShadow(QtWidgets.QFrame.Raised)
        self.heatMap.setObjectName("heatMapFrame")


        self.groupBox = QtWidgets.QGroupBox(self.MainWindow)
        # self.groupBox.setGeometry(QtCore.QRect(390, 320, 261, 131))
        self.groupBox.setGeometry(QtCore.QRect(390, 320, 261, 136))
        self.groupBox.setObjectName("groupBox")

        self.pushButton = QtWidgets.QPushButton(self.groupBox)
        self.pushButton.setGeometry(QtCore.QRect(140, 100, 113, 32))
        self.pushButton.setObjectName("pushButton")

        self.checkBox_4 = QtWidgets.QCheckBox(self.groupBox)
        self.checkBox_4.setGeometry(QtCore.QRect(130, 70, 131, 21))
        self.checkBox_4.setObjectName("checkBox_4")

        self.checkBox_3 = QtWidgets.QCheckBox(self.groupBox)
        self.checkBox_3.setGeometry(QtCore.QRect(10, 70, 141, 21))
        self.checkBox_3.setObjectName("checkBox_3")

        self.checkBox_2 = QtWidgets.QCheckBox(self.groupBox)
        self.checkBox_2.setGeometry(QtCore.QRect(130, 40, 87, 20))
        self.checkBox_2.setObjectName("checkBox_2")
        
        self.checkBox = QtWidgets.QCheckBox(self.groupBox)
        self.checkBox.setGeometry(QtCore.QRect(10, 40, 111, 21))
        self.checkBox.setObjectName("checkBox")


        # rerun button
        self.reRunButton = QtWidgets.QPushButton(self.MainWindow)
        self.reRunButton.setGeometry(QtCore.QRect(280, 420, 100, 32))
        self.reRunButton.setObjectName("reRunButton")
        


        # gene gene button
        self.geneGeneButton = QtWidgets.QPushButton(self.MainWindow)
        self.geneGeneButton.setGeometry(QtCore.QRect(20, 420, 113, 32))
        self.geneGeneButton.setObjectName("geneGeneButton")
        # self.geneGeneButton.clicked.connect(self.genePopUp(self))

        # self.geneGenePopUp = QtWidgets.QMessageBox(self.wid)

        # def genePopUp(self):
            # geneGenePopUp.show()
        #self.geneGeneButton.adjustSize()
        


        # self.scrollArea = QtWidgets.QScrollArea(MainWindow)
        # self.scrollArea.setGeometry(QtCore.QRect(40, 330, 241, 91))
        # self.scrollArea.setWidgetResizable(True)
        # self.scrollArea.setObjectName("scrollArea")

        self.dropdown = QtWidgets.QComboBox(self.MainWindow)
        # self.dropdown.setGeometry(QtCore.QRect(30, 300, 239, 89))
        self.dropdown.setGeometry(QtCore.QRect(37, 330, 239, 29))

        def pathwayNames(d):
            names = []
            for pathway in d.values():
                names.extend([*pathway])
            return names
        
        # for example 
        dict = {'sonia': {'brown': 'tan', 'braids': 'tired'}, 'elijah': {'black': 'blue', 'blonde': 'yawn'}, 'lucas': {'nike': 'finland', 'curly': 'sniffle'}, 'raffa': {'skate': 'fast', 'locs': 'long'}, 'amani': {'sleep': 'happy', 'work': 'sad'}}
        

        self.dropdown.addItems(pathwayNames(dict))
        self.dropdown.setObjectName("pathwayDropdown")
        # self.dropdown.setWidget(self.scrollAreaWidgetContents)


        # self.scrollAreaWidgetContents = QtWidgets.QWidget()
        # self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 239, 89))
        # self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        # self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        # self.label = QtWidgets.QLabel(self.MainWindow)
        # self.label.setGeometry(QtCore.QRect(140, 270, 41, 16))
        # self.label.setObjectName("label")


        self.label_2 = QtWidgets.QLabel(self.MainWindow)
        self.label_2.setGeometry(QtCore.QRect(460, 270, 100, 16))
        self.label_2.setObjectName("label_2")


        self.label_3 = QtWidgets.QLabel(self.MainWindow)
        self.label_3.setGeometry(QtCore.QRect(40, 310, 121, 16))
        self.label_3.setObjectName("label_3")


        self.pushButton_2 = QtWidgets.QPushButton(self.MainWindow)
        self.pushButton_2.setGeometry(QtCore.QRect(142, 420, 131, 32))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.pressed.connect(self.test_umap)


        self.groupBox.raise_()
        self.umap.raise_()
        self.heatMap.raise_()
        # self.scrollArea.raise_()
        # self.label.raise_()
        self.label_2.raise_()
        self.label_3.raise_()
        self.pushButton_2.raise_()
        self.geneGeneButton.raise_()
        self.reRunButton.raise_()
        self.umap.raise_()

        self.retranslateUi()
        QtCore.QMetaObject.connectSlotsByName(self.MainWindow)

    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.MainWindow.setWindowTitle(_translate("MainWindow", "SPECTRA"))
        self.groupBox.setTitle(_translate("MainWindow", "Save Options"))
        self.pushButton.setText(_translate("MainWindow", "Save "))
        self.checkBox_4.setText(_translate("MainWindow", "SPECTRA Model"))
        self.checkBox_3.setText(_translate("MainWindow", "Umap"))
        self.checkBox_2.setText(_translate("MainWindow", "Heat map"))
        self.checkBox.setText(_translate("MainWindow", "AnnData"))
        # self.label.setText(_translate("MainWindow", "Umap "))
        self.label_2.setText(_translate("MainWindow", "Heat map"))
        self.label_3.setText(_translate("MainWindow", "Umap Coloration"))
        self.pushButton_2.setText(_translate("MainWindow", "Recolor Umap"))
        self.geneGeneButton.setText(_translate("MainWindow", "Gene-gene graph"))
        self.reRunButton.setText(_translate("MainWindow", "Run Again"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    
    ui = OutputPage2(app.primaryScreen().size().width(), app.primaryScreen().size().height())
    
    ui.MainWindow.show()
    sys.exit(app.exec_())
