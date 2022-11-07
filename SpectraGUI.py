 # import our libraries
import sys
import PyQt5.QtWidgets as qw
import InputPage as ip
import OutputPage2 as op
import scanpy as sc
from spectra import spectra as spc

class SpectraGUI:

    def __init__(self, height, width) -> None:
        self.screen_height = height
        self.screen_width = width
        self.input = ip.InputPage(self.screen_height, self.screen_width)
        self.input.run_button.clicked.connect(self.run_and_quit)


    def run_and_quit(self):
        if(self.input.run_button_press()):
            anndata_path = self.input.annd_box.text()
            gene_dict = self.input.path_ann.genes_dict
            cell_type_key = self.input.ctype_box.text()

            try:
                lambda_val = float(self.input.lam_box.text())
            except:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("Lambda value isn't a number - try again")
                newWindow.exec()

            highly_var = self.input.adv_op.hv_label.isChecked()
            try:
                rho_val = float(self.input.adv_op.rho_box.text())
            except:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("Rho value isn't a number - try again")
                newWindow.exec()
            try:
                delta_val = float(self.input.adv_op.delta_box.text())
            except:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("Delta value isn't a number - try again")
                newWindow.exec()
            try:
                kappa_val = float(self.input.adv_op.kappa_box.text())
            except:
                newWindow = qw.QMessageBox(self.input.wid)
                newWindow.setText("Kappa value isn't a number - try again")
                newWindow.exec()
            use_weights = self.input.adv_op.weights_label.isChecked()
            top_genes = int(self.input.adv_op.genes_box.text())

            self.input.hide()

            adata = sc.read_h5ad(anndata_path)
            model = spc.est_spectra(adata = self.adata,  gene_set_dictionary = gene_dict, cell_type_key = cell_type_key, lam = lambda_val, use_highly_variable = highly_var, 
                                    rho = rho_val, delta = delta_val, kappa = kappa_val, use_weights = use_weights, n_top_vals = top_genes)

            self.output = op.OutputPage2(adata, model, self.screen_height, self.screen_width)
            self.output.reRunButton.clicked.connect(self.run_spectra_again)
            self.output.MainWindow.show()
  
    def run_spectra_again(self):
        self.output.MainWindow.hide()
        self.input.show()
    
if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = SpectraGUI(app.primaryScreen().size().height(), app.primaryScreen().size().width())
    sys.exit(app.exec_())


