import scanpy as sc
import pandas as pd
import seaborn as sb
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
# matplotlib.use("Qt5Agg")
from spectra import spectra as spc

from html2image import Html2Image

from PyQt5 import QtCore, QtGui, QtWidgets
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
sb.set(font_scale = 0.75)

class OutputPage2(object):

    def __init__(self, screen_width, screen_height, anndata = None, model = None, gene_sets = None, cell_type_key = None):
        sc.set_figure_params(facecolor="F0F0F0")

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

    def draw_heatmap(self):

        if self.anndata:
            df = pd.DataFrame(self.anndata.obsm["SPECTRA_cell_scores"])
            df["cell_type"] = self.anndata.obs[self.cell_type_key].values
            df = df.groupby("cell_type").mean()

            g = sb.clustermap(df, row_cluster = False, xticklabels = 0, col_cluster = True, dendrogram_ratio = (0, 0), cbar_pos = None, standard_scale = 0, linewidth = 0)

            g.ax_heatmap.set_xlabel("Factors")
            g.ax_heatmap.set_ylabel("Cell type")

            g.figure.colorbar(g.ax_heatmap.collections[0], ax = g.ax_heatmap, location = 'top', fraction = 0.05, pad = 0.05)


            fig = g.figure
            
            return fig

        else:
            return plt.Figure()

    def recolor_umap(self):
        
        self.curr_factor = self.dropdown.currentIndex()

        if self.anndata:
            self.redraw_umap()

    def draw_umap(self):

        self.point_size = ((self.width * self.height) - 800000) / 100000

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
        self.ax.set_title("UMAP")
            
        self.umap_plot = self.ax.scatter(self.anndata.obsm["X_umap"][:,0], self.anndata.obsm["X_umap"][:,1], c = self.anndata.obsm["SPECTRA_cell_scores"][:,self.curr_factor], s = self.point_size, alpha = 0.9, cmap = "inferno")
        self.umap_canvas.figure.colorbar(self.umap_plot, ax = self.ax, pad = 0.01, format = "%4.2e")
        self.umap_canvas.draw()

    def genePopUp(self):
        if self.model == None:

            newWindow = QtWidgets.QMessageBox(self.MainWindow)
            newWindow.setText("No Model  loaded")
            newWindow.exec()
        
        else:
            big_graph = self.model.return_graph()

            # this works on kate's computer, with the local try/except for graph_network
            # needs more work
            # genes_set = self.anndata.var_names[self.anndata.var["spectra_vocab"]][0]
            # genes_set = []
            # out = spc.graph_network(self.anndata, big_graph, genes_set)
            # out.show("ggGraph.html")
            
            # hti = Html2Image()
            # hti.screenshot(html_file='ggGraph.html', save_as='ggGraph.png')
            
            self.photo_window = ggg_window()
            self.photo_window.show()

    def setupUi(self):

        self.MainWindow.setObjectName("MainWindow")
        self.MainWindow.resize(self.width, self.height)

        self.main_layout = QtWidgets.QGridLayout()

        self.plots = QtWidgets.QHBoxLayout()


        self.umap_canvas = FigureCanvasQTAgg(plt.Figure())
        self.draw_umap()

        self.heatmap_canvas = FigureCanvasQTAgg(self.draw_heatmap())

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


        def colorByFactor():
            if self.anndata:
                factorList = []

                for i in range(len(self.anndata.uns["SPECTRA_markers"])):
                    factorNames = ', '.join(self.anndata.uns["SPECTRA_markers"][i][:5])
                    factorString = "Factor " + str(i) + ": " + factorNames
                    factorList.append(factorString)

                return factorList
            else:
                return ["Factor 0", "Factor 1", "Factor 2", "Factor 3"]



        self.dropdown.addItems(colorByFactor())

        self.output_options.addWidget(self.dropdown, 1, 0, 1, 2)


       


        self.vmin_max = QtWidgets.QHBoxLayout()

        self.vmin_label = QtWidgets.QLabel("v-min:")
        self.vmin_label.setFont(QtGui.QFont("Times", 11))
        self.vmin_label.setMaximumWidth(self.width // 35)

        self.vmin_box = QtWidgets.QLineEdit()
        self.vmin_box.setText("0")
        self.vmin_box.setFont(QtGui.QFont("Times", 11))
        self.vmin_box.setMaximumWidth(self.width // 40)

        self.vmin_max.addWidget(self.vmin_label)
        self.vmin_max.addWidget(self.vmin_box)

        self.vmax_label = QtWidgets.QLabel("v-max:")
        self.vmax_label.setFont(QtGui.QFont("Times", 11))
        self.vmax_label.setMaximumWidth(self.width // 35)

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
        self.output_options.setRowMinimumHeight(0, int( self.height * 1 / 3 / 3))
        self.output_options.setRowMinimumHeight(2, int( self.height * 1 / 3 / 5))

        self.output_options_frame.setLayout(self.output_options)
        self.output_options_frame.setMaximumWidth(int(self.width / 2.5))

        self.main_layout.addWidget(self.output_options_frame, 1, 0)


        self.save_options_frame = QtWidgets.QFrame()
        # self.save_options_frame.setFrameStyle(QtWidgets.QFrame.StyledPanel)
        # self.save_options_frame.setLineWidth(2)

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
            if self.checkBox_adata.isChecked() == True:
                self.anndata.write_h5ad("new_annData.h5ad")

            if self.checkBox_model.isChecked() == True:
                print("model save")
                # self.model.save("SPECTRA_model")

            if self.checkBox_umap.isChecked() == True:
                self.umap_canvas.print_figure("UMAP_figure.png")

            if self.checkBox_heatmap.isChecked() == True:
                self.heatmap_canvas.print_figure("Heatmap_figure.png")

        self.save_button = QtWidgets.QPushButton()
        self.save_button.setText("Save")
        self.save_button.clicked.connect(saveData)

        self.save_options.addWidget(self.save_button, 3, 2)

        # self.save_options.setHorizontalSpacing(int( (self.height * 5 / 12) / 5))

        self.save_options_frame.setLayout(self.save_options)

        self.main_layout.addWidget(self.save_options_frame, 1, 2)


        self.main_layout.setRowMinimumHeight(0, int(self.height * 2 / 3))

        self.main_layout.setColumnMinimumWidth(0, self.width // 100)
        self.main_layout.setColumnStretch(0, 0)
        self.main_layout.setColumnMinimumWidth(1, int(self.width / 3))


        self.MainWindow.setLayout(self.main_layout)


class ggg_window(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Gene-Gene Graph"
        self.left = 100
        self.top = 100
        self.width = 400
        self.height = 400
        self.newgggWindow()

    def newgggWindow(self):
        # make the window
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # make central widget
        self.wid = QtWidgets.QWidget(self)
        self.setCentralWidget(self.wid)

        
        photo = QtWidgets.QLabel(self.wid)
        photo.setText("")
        try:
            photo.setPixmap(QtGui.QPixmap("ggGraph.png"))
        except:
            photo.setText("Image not Found")
        photo.setScaledContents(True)
        photo.setObjectName("photo")




if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    
    ui = OutputPage2(app.primaryScreen().size().width(), app.primaryScreen().size().height())
    
    ui.MainWindow.show()
    sys.exit(app.exec_())
