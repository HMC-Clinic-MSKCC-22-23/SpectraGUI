# import our libraries
import sys
import PyQt5.QtWidgets as qw
import InputPage as ip
import OutputPage as op
import scanpy as sc
from spectra import spectra as spc

class SpectraGUI:

    def __init__(self) -> None:
        self.input = ip.InputPage()
        self.input.run_button.clicked.connect(self.run_and_quit)

        self.output = op.OutputPage()


    def run_and_quit(self):
        if(self.input.run_button_press()):
            anndata = sc.read_h5ad(self.input.annd_box.text())
            gene_dict = self.input.path_ann.genes_dict

            lambda_val = int(self.input.lam_box.text())

            highly_var = self.input.adv_op.hv_label.isChecked()
            rho_val = int(self.input.adv_op.rho_box.text())
            delta_val = int(self.input.adv_op.delta_box.text())
            kappa_val = int(self.input.adv_op.kappa_box.text())
            use_weights = self.input.adv_op.weights_label.isChecked()
            top_genes = int(self.input.adv_op.genes_box.text())

            model = spc.est_spectra(adata = anndata,  gene_set_dictionary = gene_dict, cell_type_key = "cell_type", use_highly_variable = highly_var, lam = lambda_val,
                                    rho = rho_val, delta = delta_val, kappa = kappa_val, use_weights = use_weights, n_top_vals = top_genes)

            self.input.hide()
            self.output.show()
    
if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = SpectraGUI()
    sys.exit(app.exec_())


