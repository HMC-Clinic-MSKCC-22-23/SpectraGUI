# import our libraries
import sys, os, yaml
import scanpy as sc
import PyQt5.QtWidgets as qw
from PyQt5.QtCore import QThread

import InputPage as ip
import OutputPage as op

# find spectra installation
with open('config.yaml', 'r') as file:
    config = yaml.safe_load(file)

sys.path.append(os.path.join(sys.path[0], config["spectra_path"]))

if config["gpu_enable"]:
    from spectra import spectra_gpu as spc
else:
    from spectra import spectra as spc

class ModelBuilder(QThread):
    
    def __init__(self, gui_class):
        QThread.__init__(self)
        self.gui_class = gui_class

    def run(self):

        print("The model building has begun. Standby for a progress bar.")

        self.model = spc.est_spectra(adata = self.gui_class.adata,  gene_set_dictionary = self.gui_class.gene_dict, cell_type_key = self.gui_class.cell_type_key, 
                                     lam = self.gui_class.lambda_val, use_highly_variable = self.gui_class.highly_var, rho = self.gui_class.rho_val, 
                                     delta = self.gui_class.delta_val, kappa = self.gui_class.kappa_val, use_weights = self.gui_class.use_weights, 
                                     n_top_vals = self.gui_class.top_genes, num_epochs = self.gui_class.epochs)

class SpectraGUI:

    def __init__(self, width, height) -> None:
        self.screen_width = width
        self.screen_height = height
        self.input = ip.InputPage(self.screen_width, self.screen_height)
        self.input.run_button.clicked.connect(self.run_and_quit)


    def run_and_quit(self):
        if(self.input.run_button_press()):
            anndata_path = self.input.annd_box.text()
            self.cell_type_key = self.input.ctype_box.text()


            self.gene_dict = self.input.path_ann.genes_dict

            try:
                self.lambda_val = float(self.input.lam_box.text())
            except:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("Lambda value isn't a number - try again")
                newWindow.exec()
                return

            try:
                self.epochs = int(self.input.epochs_box.text())
            except:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("'Number of epochs' must be an integer - try again")
                newWindow.exec()
                return

            self.highly_var = self.input.adv_op.hv_label.isChecked()

            try:
                self.rho_val = float(self.input.adv_op.rho_box.text())
            except:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("Rho value isn't a number - try again")
                newWindow.exec()
                return

            try:
                self.delta_val = float(self.input.adv_op.delta_box.text())
            except:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("Delta value isn't a number - try again")
                newWindow.exec()
                return

            try:
                self.kappa_val = float(self.input.adv_op.kappa_box.text())
            except:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("Kappa value isn't a number - try again")
                newWindow.exec()
                return

            self.use_weights = self.input.adv_op.weights_label.isChecked()
            
            try:
                self.top_genes = int(self.input.adv_op.genes_box.text())
            except:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("'Top genes per factor' isn't a number - try again")
                newWindow.exec()
                return


            self.adata = sc.read_h5ad(anndata_path)

            if self.cell_type_key not in self.adata.obs_keys():
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("Cell type key not found in adata.obs - try again")
                newWindow.exec()
                return
            
            self.check_gene_dict()

            if len(self.gene_dict["global"]) == 0:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("The provided 'global' gene set cannot be empty - try again")
                newWindow.exec()
                return

            self.workerThread = ModelBuilder(self)
            self.workerThread.finished.connect(self.end_thread)
            self.workerThread.start()

            newWindow = qw.QMessageBox(self.input.wid)
            newWindow.setText("Running in progress. Check the console for progress.")
            newWindow.exec()

            self.input.hide()
            


    def end_thread(self):

        self.output = op.OutputPage(self.screen_width, self.screen_height, anndata = self.adata, model = self.workerThread.model, 
                                    cell_type_key = self.cell_type_key)
        self.output.reRunButton.clicked.connect(self.run_spectra_again)
        self.output.MainWindow.show()

        self.workerThread.deleteLater()


    def run_spectra_again(self):
        reply = qw.QMessageBox.question(self.output.MainWindow, "Rerun", "Are you sure you want to rerun the training process? This will discard the trained model and updated Annotated Data.",
                                        qw.QMessageBox.Yes | qw.QMessageBox.No, qw.QMessageBox.No)
        if reply == qw.QMessageBox.Yes:
            self.output.MainWindow.close()
            self.input.show()
    
    def check_gene_dict(self):
        # must contain "global" key in addition to every unique cell type under .obs.<cell_type_key>
        gene_set_list = list(self.gene_dict.keys())
        for cell_type in self.adata.obs[self.cell_type_key].unique():
            if cell_type not in gene_set_list:
                self.gene_dict[cell_type] = {}

        if "global" not in gene_set_list:
            self.gene_dict["global"] = {}


if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    app.setStyle("Windows")
    ex = SpectraGUI(app.primaryScreen().size().width(), app.primaryScreen().size().height())
    sys.exit(app.exec_())
