# import our libraries
import sys
import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc
import PyQt5.QtGui as qg
import spectra as spc
import scanpy as sc
import matplotlib
import matplotlib.pyplot
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
matplotlib.use("Qt5Agg")

class OutputPage(qw.QMainWindow):

    def __init__(self, anndata_path = None, gene_dict = None, lambda_val = 0.0001, highly_var = False, rho_val = 0, delta_val = 0, kappa_val = 0, use_weights = False, top_genes = 1):

        sc.set_figure_params(facecolor="F0F0F0")

        if anndata_path != None:
            anndata = sc.read_h5ad(anndata_path)
            self.umap_figure = sc.pl.umap(anndata, title = "UMAP", na_color = "gray", frameon = False, use_raw = False, show = False, return_fig= True)
            self.adata = anndata
        else:
            self.umap_figure = matplotlib.pyplot.figure()

     


        self.gene_dictionary = gene_dict
        self.lam = lambda_val
        self.highly_var = highly_var
        self.rho = rho_val
        self.delta = delta_val
        self.kappa = kappa_val
        self.use_weights = use_weights
        self.top_genes = top_genes

        super().__init__()
        self.title = "Spectra output screen"
        self.left = 100
        self.top = 100
        self.width = 1200
        self.height = 800
        self.initOutputUI()
    
    def initOutputUI(self):
        # make the window
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # make central widget
        self.wid = qw.QWidget(self)
        self.setCentralWidget(self.wid)
        layout = qw.QGridLayout()

        # make save widget (the one on the bottom)
        saveWidget = qw.QWidget(self)
        saveLayout = qw.QGridLayout()

        # top layout widget
        outputWidgets = qw.QWidget(self)
        topWidgetsLayout = qw.QGridLayout()

        # make pathway annotations widget (top left)
        pathwayWidget = qw.QWidget(self)
        pathwayOptionsLayout = qw.QGridLayout()

        # make umap/tsne widget (top middle)
        # maybe this doesn't need to be a widget
        # will be label box for now

        # make gene-gene graph widget (top right)
        # same as above, will be label box for now


        # save options
        save_this = qw.QLabel("Save this?")
        save_this.setFrameShape(qw.QFrame.Box)
        saveLayout.addWidget(save_this,0,0,1,1,qc.Qt.AlignCenter)

        save_that = qw.QLabel("Save that?")
        save_that.setFrameShape(qw.QFrame.Box)
        saveLayout.addWidget(save_that,1,0,1,1,qc.Qt.AlignCenter)



        # pathway options
        # this is also going to be a box for rn
        pathway_options_temp = qw.QLabel("NEEDS TO BE MAYBE ONE OF MANY CHECKBOXES?")
        pathway_options_temp.setFrameShape(qw.QFrame.Box)
        pathwayOptionsLayout.addWidget(pathway_options_temp,0,0,1,1,qc.Qt.AlignCenter)

        p2_temp = qw.QLabel("ANOTHER PATHWAY OPTION")
        p2_temp.setFrameShape(qw.QFrame.Box)
        pathwayOptionsLayout.addWidget(p2_temp,1,0,1,1,qc.Qt.AlignCenter)

        pathwayWidget.setLayout(pathwayOptionsLayout)
        topWidgetsLayout.addWidget(pathwayWidget,0,0,1,1)

        # umap/tsne map
        umap_temp = FigureCanvasQTAgg(self.umap_figure)
        topWidgetsLayout.addWidget(umap_temp,0,1,1,1,qc.Qt.AlignCenter) # add to top widget

        # gene-gene graph
        ggg_temp = qw.QLabel("GENE-GENE GRAPH HERE")
        ggg_temp.setFrameShape(qw.QFrame.Box)
        topWidgetsLayout.addWidget(ggg_temp,0,2,1,1,qc.Qt.AlignCenter) # add to top widget


        # add top and bottom widgets to top-most level layout
        outputWidgets.setLayout(topWidgetsLayout)
        saveWidget.setLayout(saveLayout)

        layout.addWidget(outputWidgets,0,0,1,1)
        layout.addWidget(saveWidget, 1,0,1,1)

        self.wid.setLayout(layout)

if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = OutputPage()
    ex.show()
    sys.exit(app.exec_())
