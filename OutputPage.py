
import os
from turtle import color
#from spectra import spectra as spc
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from PyQt5 import QtCore, QtGui, QtWidgets, QtWebEngineWidgets
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
sb.set(font_scale = 0.75)

class OutputPage(object):

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

        self.point_size = ((self.width * self.height) - 700000) / 100000
    
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
        try:
            self.point_size = float(self.point_size_box.text())
        except:
            newWindow = QtWidgets.QMessageBox()
            newWindow.setText("You must input a valid point size.")
            newWindow.exec()
            return

        if self.anndata:
            self.redraw_umap()

    def recolor_umap_by_gene(self):
        
        try:
            self.point_size = float(self.point_size_box.text())
        except:
            newWindow = QtWidgets.QMessageBox()
            newWindow.setText("You must input a valid point size.")
            newWindow.exec()
            return

        try:
            gene_text = self.gene_color_popup.input.text()
            gene_list = ( "".join(gene_text.split()) ).split(",")
            plot_var = np.sum(self.anndata.X[:, [self.anndata.var_names.get_loc(gene) for gene in gene_list]],axis=1)
            
            plot_var = plot_var.flatten().tolist()[0]
        except:
            newWindow = QtWidgets.QMessageBox()
            newWindow.setText("Genes not found in AnnData. Please use the exact name, and separate names with a comma.")
            newWindow.exec()
            return



        self.umap_canvas.figure.clear()
        self.ax = self.umap_canvas.figure.subplots()
        self.ax.grid(False)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_title(f"UMAP: {len(gene_list)} Genes")

        vmin = min(plot_var)
        vmax = max(plot_var)
        vdomain = vmax - vmin

        try:
            vmax = vmin + (vdomain * (float(self.vmax_box.text()) / 100))
        except:
            vmax = vmin + vdomain

        try:
            vmin += vdomain * (float(self.vmin_box.text()) / 100) 
        except:
            vmin += 0
            
            
        self.umap_plot = self.ax.scatter(self.anndata.obsm["X_umap"][:,0], self.anndata.obsm["X_umap"][:,1], vmin = vmin, vmax = vmax, c = plot_var, s = self.point_size, cmap = "viridis")
        self.umap_canvas.figure.colorbar(self.umap_plot, ax = self.ax, pad = 0.01, fraction = 0.05, format = "%4.3f")
        self.umap_canvas.draw()

    def draw_umap(self):

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
        self.umap_canvas.figure.colorbar(self.umap_plot, ax = self.ax, pad = 0.01, fraction = 0.05, format = "%4.3f")
        self.umap_canvas.draw()

    def genePopUp(self):
        if not self.model:

            newWindow = QtWidgets.QMessageBox(self.MainWindow)
            newWindow.setText("No Model  loaded")
            newWindow.exec()
        
        else:
            big_graph = self.model.return_graph()

            if not self.curr_factor:
                self.curr_factor = self.dropdown.currentIndex()

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

    def geneColorPopUp(self):
        self.gene_color_popup = umap_gene_window(self.recolor_umap_by_gene)
        self.gene_color_popup.show()

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

        if self.anndata:
            self.factor_list = list(range(len(self.anndata.uns["SPECTRA_markers"])))

        self.umap_canvas = FigureCanvasQTAgg(plt.Figure())
        self.draw_umap()

        self.umap_canvas.setContentsMargins(0,0,0,0)

        self.heatmap_canvas = FigureCanvasQTAgg(self.draw_heatmap(self.cell_type_key))

        self.heatmap_canvas.setContentsMargins(0,0,0,0)

        self.plots.addWidget(self.umap_canvas)
        self.plots.addWidget(self.heatmap_canvas)

        self.main_layout.addLayout(self.plots, 0, 0, 1, 5)

        # frame to hold the output options
        self.output_options_frame = QtWidgets.QFrame()

        self.output_options = QtWidgets.QGridLayout()

        self.dropdown_label = QtWidgets.QLabel("UMAP coloration")
        self.dropdown_label.setFont(QtGui.QFont("Times", 11))

        self.output_options.addWidget(self.dropdown_label, 0, 0)

        self.dropdown = QtWidgets.QComboBox(self.MainWindow)



        self.dropdown.addItems(self.colorByFactor())

        self.output_options.addWidget(self.dropdown, 1, 0, 1, 2)


        self.plot_options = QtWidgets.QHBoxLayout()

        self.vmin_label = QtWidgets.QLabel("v-min:")
        self.vmin_label.setFont(QtGui.QFont("Times", 11))
        self.vmin_label.setMaximumWidth(self.width // 30)

        self.vmin_box = QtWidgets.QLineEdit()
        self.vmin_box.setText("0")
        self.vmin_box.setFont(QtGui.QFont("Times", 11))
        self.vmin_box.setMaximumWidth(self.width // 40)
        

        self.plot_options.addWidget(self.vmin_label, 12)
        self.plot_options.addWidget(self.vmin_box, 8)

        self.vmax_label = QtWidgets.QLabel("v-max:")
        self.vmax_label.setFont(QtGui.QFont("Times", 11))
        self.vmax_label.setMaximumWidth(self.width // 28)

        self.vmax_box = QtWidgets.QLineEdit()
        self.vmax_box.setText("100")
        self.vmax_box.setFont(QtGui.QFont("Times", 11))
        self.vmax_box.setMaximumWidth(self.width // 40)

        self.plot_options.addWidget(self.vmax_label, 12)
        self.plot_options.addWidget(self.vmax_box, 8)

        self.point_size_label = QtWidgets.QLabel("point size:")
        self.point_size_label.setFont(QtGui.QFont("Times", 11))
        self.point_size_label.setMaximumWidth(self.width // 19)

        self.point_size_box = QtWidgets.QLineEdit()
        self.point_size_box.setText(str(int(self.point_size)))
        self.point_size_box.setFont(QtGui.QFont("Times", 11))
        self.point_size_box.setMaximumWidth(self.width // 45)

        self.plot_options.addWidget(self.point_size_label, 52)
        self.plot_options.addWidget(self.point_size_box, 8)

        self.plot_options.setSpacing(10)

        self.output_options.addLayout(self.plot_options, 2, 0, 1, 2)


        self.recolor_button = QtWidgets.QPushButton(self.MainWindow)
        self.recolor_button.setText("Redraw UMAP")
        self.recolor_button.pressed.connect(self.recolor_umap)

        self.output_options.addWidget(self.recolor_button, 3, 0)
    
        self.geneGeneButton = QtWidgets.QPushButton(self.MainWindow)
        self.geneGeneButton.setText("Gene-Gene Graph")
        self.geneGeneButton.clicked.connect(self.genePopUp)

        self.output_options.addWidget(self.geneGeneButton, 3, 1)

        gene_color_button = QtWidgets.QPushButton(self.MainWindow)
        gene_color_button.setText("Color UMAP by Gene")

        self.gene_color_popup = None

        gene_color_button.clicked.connect(self.geneColorPopUp)

        self.output_options.addWidget(gene_color_button, 3, 2)

        self.output_options.setHorizontalSpacing(15)

        # self.reRunButton = QtWidgets.QPushButton(self.MainWindow)
        # self.reRunButton.setText("Run again")
        # self.output_options.addWidget(self.reRunButton, 3, 2)

        # self.output_options.setHorizontalSpacing(int( (self.height * 5 / 12) / 4))
 
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

        self.heatmap_options.setRowMinimumHeight(2, self.vmin_box.sizeHint().height() + 9)

        heatmap_factors_button = QtWidgets.QPushButton(self.MainWindow)
        heatmap_factors_button.setText("Add / Remove Factors")

        self.factor_popup = None

        heatmap_factors_button.clicked.connect(self.factorPopUp)

        self.heatmap_options.addWidget(heatmap_factors_button, 3, 0)

        # self.heatmap_options.setRowStretch(0, 1)
        # self.heatmap_options.setRowStretch(1, 0)
        # self.heatmap_options.setRowStretch(2, 0)

        self.heatmap_options_frame.setLayout(self.heatmap_options)

        self.main_layout.addWidget(self.heatmap_options_frame, 1, 2)


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

        self.save_options.addWidget(self.checkBox_model, 2, 0)

        self.checkBox_cellscore = QtWidgets.QCheckBox()
        self.checkBox_cellscore.setText("Cell Score Matrix")

        self.save_options.addWidget(self.checkBox_cellscore, 1, 1)

        self.checkBox_genescore = QtWidgets.QCheckBox()
        self.checkBox_genescore.setText("Gene Score Matrix")

        self.save_options.addWidget(self.checkBox_genescore, 2, 1)

        self.checkBox_umap = QtWidgets.QCheckBox()
        self.checkBox_umap.setText("UMAP figure")

        self.save_options.addWidget(self.checkBox_umap, 1, 2)

        self.checkBox_heatmap = QtWidgets.QCheckBox()
        self.checkBox_heatmap.setText("Heatmap figure")
        
        self.save_options.addWidget(self.checkBox_heatmap, 2, 2)


        def saveData():
            if self.checkBox_adata.isChecked():
                self.anndata.write_h5ad("new_adata.h5ad")

            if self.checkBox_model.isChecked():
                self.model.save("SPECTRA_model.pt")

            if self.checkBox_cellscore.isChecked():
                cell_scores = pd.DataFrame(self.anndata.obsm["SPECTRA_cell_scores"])
                cell_scores.columns = [f"Factor_{x}" for x in range(len(cell_scores.columns))]
                cell_scores.index = self.anndata.obs_names
                cell_scores.to_csv("cell_scores.csv", header=True, index = True)

            if self.checkBox_genescore.isChecked():
                gene_scores = pd.DataFrame(self.anndata.uns["SPECTRA_factors"]).T
                gene_scores.columns = [f"Factor_{x}" for x in range(len(gene_scores.columns))]
                gene_scores.index = self.anndata.var_names[:len(gene_scores.index)]
                gene_scores.to_csv("gene_scores.csv", header=True, index = True)

            if self.checkBox_umap.isChecked():
                self.umap_canvas.print_figure("UMAP_figure.png")

            if self.checkBox_heatmap.isChecked():
                self.heatmap_canvas.print_figure("Heatmap_figure.png")

        self.save_button = QtWidgets.QPushButton()
        self.save_button.setText("Save")
        self.save_button.clicked.connect(saveData)

        self.save_options.addWidget(self.save_button, 3, 2)

        self.reRunButton = QtWidgets.QPushButton(self.MainWindow)
        self.reRunButton.setText("Run again")
        self.save_options.addWidget(self.reRunButton, 3, 0)

        self.save_options.setHorizontalSpacing(20)
        self.save_options.setColumnStretch(0, 0)
        self.save_options.setColumnStretch(1, 0)
        self.save_options.setColumnStretch(2, 0)

        self.save_options_frame.setLayout(self.save_options)

        self.main_layout.addWidget(self.save_options_frame, 1, 4)


        self.main_layout.setRowMinimumHeight(0, int(self.height * 2.2 / 3))

        self.main_layout.setColumnMinimumWidth(0, self.width // 100)
        self.main_layout.setColumnStretch(0, 0)
        self.main_layout.setColumnStretch(2, 0)
        self.main_layout.setColumnStretch(4, 0)
        self.main_layout.setColumnMinimumWidth(2, int(self.width / 3))


        self.MainWindow.setLayout(self.main_layout)


class umap_gene_window(QtWidgets.QMainWindow):

    def __init__(self, colorFunction):
        super().__init__()
        self.title = "Select UMAP Genes"
        self.left = 150
        self.top = 150
        self.width = 400
        self.height = 120
        self.colorFunction = colorFunction
        self.initGeneWindow()

    def initGeneWindow(self):

        basic_layout = QtWidgets.QVBoxLayout()

        gene_label = QtWidgets.QLabel("Genes to color by:", self)
        gene_label.move(20, 10)

        self.input = QtWidgets.QLineEdit(self)
        self.input.setPlaceholderText("GENE1, GENE2, GENE3, etc.")
        self.input.resize(350, 30)
        self.input.move(20, 40)
        
        color_button = QtWidgets.QPushButton(self)
        color_button.setText("Recolor UMAP")
        color_button.pressed.connect(self.colorFunction)
        color_button.resize(120, 30)
        color_button.move(20, 80)

        self.setGeometry(150, 150, self.width, self.height)
        self.setWindowTitle(self.title)



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
        self.setWindowTitle(self.title)


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
    
    ui = OutputPage(app.primaryScreen().size().width(), app.primaryScreen().size().height())
    
    ui.MainWindow.show()
    sys.exit(app.exec_())
