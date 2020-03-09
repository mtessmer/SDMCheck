from PyQt5 import QtWidgets, QtCore, QtGui
import pyqtgraph as pg
import sys, os
from functools import partial
from Bio import SeqIO, AlignIO, pairwise2
from Bio.Align.Applications import ClustalOmegaCommandline

format_dict = {'ab1': 'abi',
               'sta': 'fasta',
               'seq': 'fasta'}

channels = ('DATA9', 'DATA10', 'DATA11', 'DATA12')

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.setWindowTitle("PySeq")
        self.resize(1200, 800)

        self.centralwidget = QtWidgets.QWidget(self)

        self.listWidget = FileList(self.centralwidget)
        self.listWidget.setGeometry(QtCore.QRect(10, 10, 500, 190))

        self.list_button_widget = QtWidgets.QWidget(self.centralwidget)
        self.list_button_widget.setGeometry(QtCore.QRect(510, 10, 90, 190))

        self.verticalLayout = QtWidgets.QVBoxLayout(self.list_button_widget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")

        self.add_file = QtWidgets.QPushButton("Add File", self.list_button_widget)
        self.verticalLayout.addWidget(self.add_file)
        self.add_file.clicked.connect(self.listWidget.add)

        self.remove_file = QtWidgets.QPushButton("Remove File", self.list_button_widget)
        self.verticalLayout.addWidget(self.remove_file)
        self.remove_file.clicked.connect(self.listWidget.remove)

        self.markRef = QtWidgets.QPushButton("Mark Reference", self.list_button_widget)
        self.verticalLayout.addWidget(self.markRef)
        self.markRef.clicked.connect(self.listWidget.mark_ref)

        self.align = QtWidgets.QPushButton("Align", self.list_button_widget)
        self.verticalLayout.addWidget(self.align)
        self.align.clicked.connect(self.run_alignment)

        self.init_menu()

        self.text = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.text.setGeometry(QtCore.QRect(10, 210, 1000, 200))
        self.text.setFont(QtGui.QFont("Courier", 15))
        self.text.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        self.highlight = SnpHighlighter(self.text.document())

        self.text_button_widget = QtWidgets.QWidget(self.centralwidget)
        self.text_button_widget.setGeometry(QtCore.QRect(1010, 202, 90, 200))
        self.verticalLayout2 = QtWidgets.QVBoxLayout(self.text_button_widget)
        self.verticalLayout2.setAlignment(QtCore.Qt.AlignTop)
        self.verticalLayout2.setSpacing(0)

        self.graphWidget = pg.PlotWidget(self.centralwidget)
        self.graphWidget.setGeometry(QtCore.QRect(10, 420, 1000, 200))
        self.graphWidget.setBackground('w')
        self.graphWidget.setRange(xRange=[300, 800])
        self.graphWidget.setMouseEnabled(y=False)

        self.fileListWidget = FileList()

        self.setCentralWidget(self.centralwidget)

    def init_menu(self):
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu('File')

        # Create actions
        quit_action = QtWidgets.QAction('Quit', self)
        quit_action.triggered.connect(self.quit_app)

        # Add actions to menu
        file_menu.addAction(quit_action)

    def run_alignment(self):

        if not self.listWidget.ref:
            return None

        self.alignment, self.seq_objs = self.listWidget.run_alignment()

        self.show_alignment()

    def show_alignment(self):
        self.button_list = []
        self.text.insertPlainText(str(self.alignment[0].seq) + '\n')
        for i, elem in enumerate(self.alignment[1:]):
            self.text.insertPlainText(str(elem.seq) + '\n')
            self.button_list.append(QtWidgets.QPushButton("View Trace"))
            self.button_list[i].clicked.connect(partial(self.show_trace, seq_id=self.listWidget.item(i+1).text()))
            self.verticalLayout2.addWidget(self.button_list[i])

    def show_trace(self, seq_id):
        self.graphWidget.clear()
        self.graphWidget.addLegend()
        seq_obj = self.seq_objs[seq_id]

        for i, channel in enumerate(channels):
            self.graphWidget.plot(seq_obj.annotations['abif_raw'][channel], pen=(i, 4), name="GATC"[i])

    def quit_app(self):
        QtWidgets.qApp.quit()

    def fileDropped(self, list_of_files):
        print('dropped:', list_of_files)


class FileList(QtWidgets.QListWidget):
    dropped = QtCore.pyqtSignal(list)
    ref = None
    def __init__(self, parent=None):
        super(FileList, self).__init__(parent)
        self.setAcceptDrops(True)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        if event.mimeData().hasUrls:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
            files = []
            for url in event.mimeData().urls():
                files.append(str(url.toLocalFile()))
            self.addItems(files)
        else:
            event.ignore()

    def add(self):
        files, _ = QtWidgets.QFileDialog.getOpenFileNames()
        self.addItems(files)

    def remove(self):
        for item in self.selectedItems():
            self.takeItem(self.row(item))

    def mark_ref(self):

        if (not self.count()) or (self.currentRow() == -1):
            return None

        for i in range(self.count()):
            self.item(i).setForeground(QtCore.Qt.black)

        self.ref = self.currentItem().text()
        i = self.currentRow()
        self.item(i).setForeground(QtCore.Qt.red)

        # Make reference first item
        self.insertItem(0, self.takeItem(i))

    def run_alignment(self):

        if not self.ref:
            return None

        seq_objs = {}
        seq_objs[self.ref] = read_file(self.ref)

        for i in range(1, self.count()):
            ref_seq = str(seq_objs[self.ref].seq)
            test_seq = self.item(i).text()
            seq_objs[test_seq] = read_file(test_seq)

            # Check reverse compliment
            taln = pairwise2.align.localxx(ref_seq, str(seq_objs[test_seq].seq), score_only=True)
            if (taln / len(ref_seq)) < 0.7:

                # Reverse Compliment ABI records GATC -->CTAG
                tmp = seq_objs[test_seq].annotations['abif_raw']
                tmp[channels[0]], tmp[channels[1]], tmp[channels[2]], tmp[channels[3]] = \
                    tmp[channels[3]][::-1], tmp[channels[2]][::-1], tmp[channels[1]][::-1], tmp[channels[0]][::-1]

                # Reverse compliment sequence
                seq_objs[test_seq] = seq_objs[test_seq].reverse_complement(id=seq_objs[test_seq].id + '_rc')
                seq_objs[test_seq].annotations['abif_raw'] = tmp.copy()

                taln = pairwise2.align.localxx(ref_seq, str(seq_objs[test_seq].seq), score_only=True)
                if (taln / len(ref_seq)) < 0.7:
                    print(self.item(i).text(), 'has low sequence homology and will not be used in the alignments')
                    del seq_objs[test_seq]

        # Write FASTA file
        SeqIO.write(seq_objs.values(), 'seqs.fasta', "fasta")

        # Perform alignment
        clustal_omega_exe = os.path.join(os.path.dirname(__file__ ), "extra\clustal-omega-1.2.2-win64\clustalo.exe")
        if os.path.exists('seqs.aln'):
            os.remove('seqs.aln')

        cline = ClustalOmegaCommandline(clustal_omega_exe, infile="seqs.fasta", outfile="seqs.aln")
        stdout, stderr = cline()

        alignment = AlignIO.read('seqs.aln', 'fasta')

        return alignment, seq_objs

class MyPlotWidget(pg.PlotWidget):

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.MidButton:  # or Qt.MiddleButton
            self.__prevMousePos = event.pos()
        else:
            super(MyPlotWidget, self).mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if event.buttons() == QtCore.Qt.MidButton:  # or Qt.MiddleButton
            offset = self.__prevMousePos - event.pos()
            self.__prevMousePos = event.pos()

            self.verticalScrollBar().setValue(self.verticalScrollBar().value() + offset.y())
            self.horizontalScrollBar().setValue(self.horizontalScrollBar().value() + offset.x())
        else:
            super(MyPlotWidget, self).mouseMoveEvent(event)

class SnpHighlighter(QtGui.QSyntaxHighlighter):

    def __init__(self, document):
        QtGui.QSyntaxHighlighter.__init__(self, document)
        self.document = document

    def highlightBlock(self, text):
        ref = self.document.toPlainText().split('/n')[0]
        if len(ref) < len(text):
            diff = len(text) - len(ref)
            ref += diff * '-'

        idx = [i for i in range(len(text)) if text[i] != ref[i]]

        for i in idx:
            self.setFormat(i, 1, QtCore.Qt.red)

        pass



def read_file(file_name):
    ext = file_name[-3:]
    if ext == 'ab1':
        return read_abi(file_name)
    elif ext == 'sta':
        return read_fasta(file_name)
    else:
        return None


def read_fasta(file_name):
    file_format = format_dict[file_name[-3:]]
    records = SeqIO.read(file_name, file_format)
    return records


def read_abi(file_name):
    file_format = format_dict[file_name[-3:]]
    records = SeqIO.read(file_name, file_format)
    return records

def main():
    app = QtWidgets.QApplication(sys.argv)
    main_window = MainWindow()

    main_window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()