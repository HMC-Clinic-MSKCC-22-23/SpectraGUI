import sys,os,csv,json
from PyQt5 import QtCore, QtGui, QtWidgets, QtPrintSupport
from PyQt5.QtWidgets import QMainWindow, QWidget, QPlainTextEdit,QApplication

class Example(QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):

        self.setGeometry(300, 300, 250, 150)
        self.setWindowTitle('CSVMarshall')

        self.genes = QPlainTextEdit(self)
        self.genes.move(300, 30)
        self.pushButtonLoad = QtWidgets.QPushButton(self)
        self.pushButtonLoad.setText("Load CSV")
        self.pushButtonLoad.clicked.connect(self.loadCsv)
        self.genes_dict = {}
        self.show()

    def loadCsv(self,fileName):
        self.genes_dict = {}
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open CSV",
                                                (QtCore.QDir.homePath()), "CSV (*.csv *.tsv)")
        ff = open(fileName, 'r')
        mytext = ff.read()

        ff.close()

        if fileName:
            f = open(fileName, 'r')
            with f:
                self.fname = os.path.splitext(str(fileName))[0].split("/")[-1]
                if mytext.count(';') <= mytext.count(','):  #tab?
                    reader = csv.reader(f, delimiter=',')
                    self.genes.clear()
                    for row in reader:
                        items = [field for field in row]

       
                        i =2
                        if items[0] not in self.genes_dict:
                            self.genes_dict[items[0]] = {}
                        if items[1] not in self.genes_dict[items[0]]:
                            self.genes_dict[items[0]][items[1]] = []
                        while len(items)>i:
                            if items[i] != "":
                                self.genes_dict[items[0]][items[1]] += [items[i]]
                            i += 1
                    self.genes.insertPlainText(str(self.genes_dict))

  
                else:
                    reader = csv.reader(f, delimiter='\t')
                    self.genes.clear()
                    for row in reader:
                        items = [field for field in row]
    
                        i = 2
                        if items[0] not in self.genes_dict:
                            self.genes_dict[items[0]] = {}
                        if items[1] not in self.genes_dict[items[0]]:
                            self.genes_dict[items[0]][items[1]] = []
                        while len(items) > i:
                            if items[i] != "":
                                self.genes_dict[items[0]][items[1]] += [items[i]]
                            i += 1
                    self.genes.insertPlainText(str(self.genes_dict))


def main():
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
