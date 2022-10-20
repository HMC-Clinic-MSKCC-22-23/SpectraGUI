# import our libraries
import sys
import PyQt5.QtWidgets as qw
import InputPage as ip
import OutputPage as op
from spectra import spectra as spc

class SpectraGUI:

    def __init__(self) -> None:
        self.input = ip.InputPage()
        self.input.run_button.clicked.connect(self.run_and_quit)

        self.output = op.OutputPage()


    def run_and_quit(self):
        self.input.run_button_press()

        spc.est_spectra()

        self.input.hide()
        self.output.show()
    
if __name__ == '__main__':
    app = qw.QApplication(sys.argv)
    ex = SpectraGUI()
    sys.exit(app.exec_())


