import scanpy as sc
import pandas as pd
import seaborn as sb
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

import os
from spectra import spectra as spc


from PyQt5 import QtCore, QtGui, QtWidgets, QtWebEngineWidgets
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
sb.set(font_scale = 0.75)

class OutputPage2(object):

    def __init__(self, screen_width, screen_height, anndata = None, model = None, gene_sets = None, cell_type_key = None):
        sc.set_figure_params(facecolor="F0F0F0")
        self.title = "SPECTRA Output"

        self.width = int(screen_width // 1.2)
        self.height = int(screen_height // 1.5)
        
        self.gene_sets = gene_sets

        self.curr_factor = None

        self.anndata = anndata

        self.photo_window = None

        self.cell_type_key = cell_type_key
        self.model = model
    
        self.MainWindow = QtWidgets.QWidget()
        self.setupUi()
    
    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:
        return super().closeEvent(a0)

    def draw_heatmap(self, y_axis):

        if self.anndata:
            df = pd.DataFrame(self.anndata.obsm["SPECTRA_cell_scores"])
            df["y_axis"] = self.anndata.obs[y_axis].values
            df = df.groupby("y_axis").mean()

            df_plot = df.iloc[:, self.factor_list]

            if len(self.factor_list) <= 20:
                x_ticks = True
            else:
                x_ticks = 0

            g = sb.clustermap(df_plot, row_cluster = False, xticklabels = x_ticks, col_cluster = True, dendrogram_ratio = (0, 0), cbar_pos = None, cmap = "viridis", standard_scale = 0, linewidth = 0)

            g.ax_heatmap.set_xlabel("Factors")
            g.ax_heatmap.set_ylabel(y_axis)

            g.figure.colorbar(g.ax_heatmap.collections[0], ax = g.ax_heatmap, location = 'top', fraction = 0.05, pad = 0.05)


            fig = g.figure
            
            return fig

        else:
            return plt.Figure()

    def redraw_heatmap(self):
        new_data = self.heatmap_dropdown.currentText()
        self.heatmap_canvas.figure.clear()
        self.plots.removeWidget(self.heatmap_canvas)
        self.heatmap_canvas = FigureCanvasQTAgg(self.draw_heatmap(new_data))
        self.plots.addWidget(self.heatmap_canvas)

    def recolor_umap(self):
        
        self.curr_factor = self.dropdown.currentIndex()

        if self.anndata:
            self.redraw_umap()

    def draw_umap(self):

        self.point_size = ((self.width * self.height) - 700000) / 100000

        self.ax = self.umap_canvas.figure.subplots()
        self.ax.grid(False)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_title("UMAP")
        # if we have loaded anndata, draw colorless umap
        if self.anndata:
            self.umap_plot = self.ax.scatter(self.anndata.obsm["X_umap"][:,0], self.anndata.obsm["X_umap"][:,1], color = "grey", s = self.point_size)
        else: # just draw a placeholder
            self.umap_plot = self.ax.plot([2,5,4,3.5,4,5,2])
    
    def redraw_umap(self):

        self.umap_canvas.figure.clear()
        self.ax = self.umap_canvas.figure.subplots()
        self.ax.grid(False)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_title("UMAP: Factor " + str(self.curr_factor))

        vmin = min(self.anndata.obsm["SPECTRA_cell_scores"][:,self.curr_factor])
        vmax = max(self.anndata.obsm["SPECTRA_cell_scores"][:,self.curr_factor])
        vdomain = vmax - vmin

        try:
            vmax = vmin + (vdomain * (float(self.vmax_box.text()) / 100))
        except:
            vmax = vmin + vdomain

        try:
            vmin += vdomain * (float(self.vmin_box.text()) / 100) 
        except:
            vmin += 0
            
            
        self.umap_plot = self.ax.scatter(self.anndata.obsm["X_umap"][:,0], self.anndata.obsm["X_umap"][:,1], vmin = vmin, vmax = vmax, c = self.anndata.obsm["SPECTRA_cell_scores"][:,self.curr_factor], s = self.point_size, cmap = "viridis")
        self.umap_canvas.figure.colorbar(self.umap_plot, ax = self.ax, pad = 0.01, format = "%4.3f")
        self.umap_canvas.draw()

    def genePopUp(self):
        if not self.model:

            newWindow = QtWidgets.QMessageBox(self.MainWindow)
            newWindow.setText("No Model  loaded")
            newWindow.exec()
        
        else:
            big_graph = self.model.return_graph()

            genes_set = []
            for i in self.anndata.uns["SPECTRA_markers"][self.curr_factor]:
                genes_set.append(str(i))
            out = spc.graph_network(self.anndata, big_graph, genes_set)
            out.show("ggGraph.html")

            self.photo_window = QtWidgets.QDialog()
            gg = ggg_window()
            gg.openWindow(self.photo_window)
            self.photo_window.show()

    def change_factors(self):

        new_list = []

        for i, checkBox in enumerate(self.factor_popup.box_list):
            if checkBox.isChecked():
                new_list.append(i)

        self.factor_list = new_list
        self.redraw_heatmap()
    
    def factorPopUp(self):
        if self.anndata:

            self.factor_popup = factor_window(self.colorByFactor(), self.factor_list, self.change_factors)
            self.factor_popup.show()

        else:
            newWindow = QtWidgets.QMessageBox(self.MainWindow)
            newWindow.setText("AnnData not properly formatted")
            newWindow.exec()
            
    def colorByFactor(self):
            if self.anndata:
                factorList = []

                for i in range(len(self.anndata.uns["SPECTRA_markers"])):
                    factorNames = ', '.join(self.anndata.uns["SPECTRA_markers"][i][:5])
                    factorString = "Factor " + str(i) + ": " + factorNames
                    factorList.append(factorString)

                return factorList
            else:
                return ["Factor 0", "Factor 1", "Factor 2", "Factor 3"]

    def setupUi(self):

        self.MainWindow.setObjectName("MainWindow")
        self.MainWindow.setWindowTitle(self.title)
        self.MainWindow.resize(self.width, self.height)

        self.main_layout = QtWidgets.QGridLayout()

        self.plots = QtWidgets.QHBoxLayout()

        self.factor_list = list(range(len(self.anndata.uns["SPECTRA_markers"])))

        self.umap_canvas = FigureCanvasQTAgg(plt.Figure())
        self.draw_umap()

        self.heatmap_canvas = FigureCanvasQTAgg(self.draw_heatmap(self.cell_type_key))

        self.plots.addWidget(self.umap_canvas)
        self.plots.addWidget(self.heatmap_canvas)

        self.main_layout.addLayout(self.plots, 0, 0, 1, 3)

        # frame to hold the output options
        self.output_options_frame = QtWidgets.QFrame()

        self.output_options = QtWidgets.QGridLayout()

        self.dropdown_label = QtWidgets.QLabel("UMAP coloration")
        self.dropdown_label.setFont(QtGui.QFont("Times", 11))

        self.output_options.addWidget(self.dropdown_label, 0, 0)

        self.dropdown = QtWidgets.QComboBox(self.MainWindow)



        self.dropdown.addItems(self.colorByFactor())

        self.output_options.addWidget(self.dropdown, 1, 0, 1, 2)


        self.vmin_max = QtWidgets.QHBoxLayout()

        self.vmin_label = QtWidgets.QLabel("v-min:")
        self.vmin_label.setFont(QtGui.QFont("Times", 11))
        self.vmin_label.setMaximumWidth(self.width // 30)

        self.vmin_box = QtWidgets.QLineEdit()
        self.vmin_box.setText("0")
        self.vmin_box.setFont(QtGui.QFont("Times", 11))
        self.vmin_box.setMaximumWidth(self.width // 40)

        self.vmin_max.addWidget(self.vmin_label)
        self.vmin_max.addWidget(self.vmin_box)

        self.vmax_label = QtWidgets.QLabel("v-max:")
        self.vmax_label.setFont(QtGui.QFont("Times", 11))
        self.vmax_label.setMaximumWidth(self.width // 30)

        self.vmax_box = QtWidgets.QLineEdit()
        self.vmax_box.setText("100")
        self.vmax_box.setFont(QtGui.QFont("Times", 11))
        self.vmax_box.setMaximumWidth(self.width // 40)

        self.vmin_max.addWidget(self.vmax_label)
        self.vmin_max.addWidget(self.vmax_box)

        self.vmin_max.setSpacing(20)
        
        self.output_options.addLayout(self.vmin_max, 2, 0)

        self.recolor_button = QtWidgets.QPushButton(self.MainWindow)
        self.recolor_button.setText("Recolor UMAP")
        self.recolor_button.pressed.connect(self.recolor_umap)

        self.output_options.addWidget(self.recolor_button, 3, 0)
    
        self.geneGeneButton = QtWidgets.QPushButton(self.MainWindow)
        self.geneGeneButton.setText("Gene-Gene Graph")
        self.geneGeneButton.clicked.connect(self.genePopUp)

        self.output_options.addWidget(self.geneGeneButton, 3, 1)

        self.reRunButton = QtWidgets.QPushButton(self.MainWindow)
        self.reRunButton.setText("Run again")

        self.output_options.addWidget(self.reRunButton, 3, 2)

        self.output_options.setHorizontalSpacing(int( (self.height * 5 / 12) / 4))
 

        self.output_options_frame.setLayout(self.output_options)
        self.output_options_frame.setMaximumWidth(int(self.width / 2.2))

        self.main_layout.addWidget(self.output_options_frame, 1, 0)


        self.heatmap_options_frame = QtWidgets.QFrame()

        self.heatmap_options = QtWidgets.QGridLayout()

        heatmap_options_label = QtWidgets.QLabel("Heatmap Y-axis")
        heatmap_options_label.setFont(QtGui.QFont("Times", 11))

        self.heatmap_options.addWidget(heatmap_options_label, 0, 0)

        self.heatmap_dropdown = QtWidgets.QComboBox(self.MainWindow)


        def getObsNames():
            if self.anndata:
                obs_list = list(self.anndata.obs.select_dtypes("category").columns.values)
                obs_list.remove(self.cell_type_key)
                obs_list.insert(0, self.cell_type_key)

                return obs_list

            else:
                return ["obs1", "obs2", "obs3"]

        self.heatmap_dropdown.addItems(getObsNames())

        self.heatmap_dropdown.currentTextChanged.connect(self.redraw_heatmap)

        self.heatmap_options.addWidget(self.heatmap_dropdown, 1, 0, 1, 2)

        heatmap_factors_button = QtWidgets.QPushButton(self.MainWindow)
        heatmap_factors_button.setText("Add / Remove Factors")

        self.factor_popup = None

        heatmap_factors_button.clicked.connect(self.factorPopUp)

        self.heatmap_options.addWidget(heatmap_factors_button)

        self.heatmap_options.setHorizontalSpacing(int( (self.height * 5 / 12) / 4))
   
        self.heatmap_options.setRowMinimumHeight(1, int(self.height / 10))


        self.heatmap_options_frame.setLayout(self.heatmap_options)

        self.main_layout.addWidget(self.heatmap_options_frame, 1, 1)


        self.save_options_frame = QtWidgets.QFrame()
       

        self.save_options = QtWidgets.QGridLayout()

        save_options_label = QtWidgets.QLabel("Save options")
        save_options_label.setFont(QtGui.QFont("Times", 11))

        self.save_options.addWidget(save_options_label, 0, 0, 1, 2)

        # connected checkboxes to save button
        # save button is connected to function that saves what is checked.

        self.checkBox_adata = QtWidgets.QCheckBox()
        self.checkBox_adata.setText("Updated AnnData")

        self.save_options.addWidget(self.checkBox_adata, 1, 0)

        self.checkBox_model = QtWidgets.QCheckBox()
        self.checkBox_model.setText("SPECTRA model")

        self.save_options.addWidget(self.checkBox_model, 1, 2)

        self.checkBox_umap = QtWidgets.QCheckBox()
        self.checkBox_umap.setText("UMAP figure")

        self.save_options.addWidget(self.checkBox_umap, 2, 0)

        self.checkBox_heatmap = QtWidgets.QCheckBox()
        self.checkBox_heatmap.setText("Heatmap figure")
        
        self.save_options.addWidget(self.checkBox_heatmap, 2, 2)


        def saveData():
            if self.checkBox_adata.isChecked():
                self.anndata.write_h5ad("new_annData.h5ad")

            if self.checkBox_model.isChecked():
                print("model save")


            if self.checkBox_umap.isChecked():
                self.umap_canvas.print_figure("UMAP_figure.png")

            if self.checkBox_heatmap.isChecked():
                self.heatmap_canvas.print_figure("Heatmap_figure.png")

        self.save_button = QtWidgets.QPushButton()
        self.save_button.setText("Save")
        self.save_button.clicked.connect(saveData)

        self.save_options.addWidget(self.save_button, 3, 2)

        self.save_options_frame.setLayout(self.save_options)

        self.main_layout.addWidget(self.save_options_frame, 1, 2)


        self.main_layout.setRowMinimumHeight(0, int(self.height * 2.2 / 3))

        self.main_layout.setColumnMinimumWidth(0, self.width // 100)
        self.main_layout.setColumnStretch(0, 0)
        self.main_layout.setColumnMinimumWidth(1, int(self.width / 3))


        self.MainWindow.setLayout(self.main_layout)


class factor_window(QtWidgets.QMainWindow):

    def __init__(self, factorList, checkedFactors, checkFunction):
        super().__init__()
        self.title = "Select Heatmap Factors"
        self.left = 150
        self.top = 150
        self.width = 600
        self.height = 300
        self.box_list = []
        self.factor_list = factorList
        self.checkedFactors = checkedFactors
        self.checkFunction = checkFunction
        self.initFactorWindow()
    
    def initFactorWindow(self):

        self.scroll = QtWidgets.QScrollArea()             # Scroll Area which contains the widgets, set as the centralWidget
        self.widget = QtWidgets.QWidget()                 # Widget that contains the collection of Vertical Box
        self.vbox = QtWidgets.QVBoxLayout()               # The Vertical Box that contains the Horizontal Boxes of  labels and buttons

        for i, label in enumerate(self.factor_list):
            box = QtWidgets.QCheckBox()
            box.setText(label)
            
            if i in self.checkedFactors:
                box.setChecked(True)
            
            box.stateChanged.connect(self.checkFunction)

            self.box_list.append(box)
            self.vbox.addWidget(box)

        self.widget.setLayout(self.vbox)

        #Scroll Area Properties
        self.scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.widget)

        self.setCentralWidget(self.scroll)

        self.setGeometry(150, 150, self.width, self.height)
        self.setWindowTitle('Scroll Area Demonstration')


class ggg_window(object):
    def openWindow(self, Dialog):
        Dialog.setObjectName("Gene-Gene Graph")
        Dialog.resize(500, 500)

        self.layout = QtWidgets.QVBoxLayout(Dialog)
        self.layout.setObjectName("layout")

        self.wid = QtWidgets.QWidget(Dialog)
        self.wid.setObjectName("widget")

        self.webEngineView = QtWebEngineWidgets.QWebEngineView(self.wid)
        self.webEngineView.load(QtCore.QUrl().fromLocalFile(os.path.abspath("ggGraph.html")))

        self.layout.addWidget(self.webEngineView)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)
    
    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Gene-Gene Graph", "Gene-Gene Graph"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    
    ui = OutputPage2(app.primaryScreen().size().width(), app.primaryScreen().size().height())
    
    ui.MainWindow.show()
    sys.exit(app.exec_())
