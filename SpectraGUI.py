 # import our libraries
import sys
import PyQt5.QtWidgets as qw
import InputPage as ip
import OutputPage2 as op
import scanpy as sc
from spectra import spectra as spc

class SpectraGUI:

    def __init__(self) -> None:
        self.input = ip.InputPage()
        self.input.run_button.clicked.connect(self.run_and_quit)


    def run_and_quit(self):
        if(self.input.run_button_press()):
            anndata_path = self.input.annd_box.text()
            gene_dict = self.input.path_ann.genes_dict
            cell_type_key = self.input.ctype_box.text()

            lambda_val = float(self.input.lam_box.text())

            highly_var = self.input.adv_op.hv_label.isChecked()
            rho_val = float(self.input.adv_op.rho_box.text())
            delta_val = float(self.input.adv_op.delta_box.text())
            kappa_val = float(self.input.adv_op.kappa_box.text())
            use_weights = self.input.adv_op.weights_label.isChecked()
            top_genes = int(self.input.adv_op.genes_box.text())

            self.input.hide()

            self.output = op.OutputPage2(anndata_path, gene_dict, cell_type_key, lambda_val, highly_var, rho_val, delta_val, kappa_val, use_weights, top_genes)

            self.output.MainWindow.show()
    
if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = SpectraGUI()
    sys.exit(app.exec_())


