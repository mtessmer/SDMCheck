from PyQt5 import QtWidgets, QtCore, QtGui
import pyqtgraph as pg
import numpy as np
from scipy.interpolate import interp1d
import sys, os
from functools import partial
from Bio import SeqIO, AlignIO, pairwise2
from Bio.Align.Applications import ClustalOmegaCommandline

format_dict = {'ab1': 'abi',
               'sta': 'fasta',
               'seq': 'fasta'}

channels = ('DATA9', 'DATA10', 'DATA11', 'DATA12')


# TODO: Comments
#       Docstrings
#       Fix Window size
#       Center Widgets
#       Add Save funciton
#       Add labels
#       Hide Trace
#       Highlight sequence of active trace & Gray out selected trace
#       Number Amino Acids
#       Shifter to adjust Amino Acid nunber
#       Codon Table With Ecoli freq
#       Live update trace alignment to edited sequence

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.setWindowTitle("SDMCheck")
        self.resize(1200, 800)
        self.pad = 10

        self.centralwidget = QtWidgets.QWidget(self)

        self.list_width =500
        self.list_height = 200
        self.listWidget = FileList(self.centralwidget)
        self.listWidget.setGeometry(QtCore.QRect(self.pad, self.pad, self.list_width, self.list_height))

        self.button_width = 90
        self.list_button_widget = QtWidgets.QWidget(self.centralwidget)
        self.list_button_widget.setGeometry(QtCore.QRect(self.list_width + 2 * self.pad, self.pad,
                                                         self.button_width, self.list_height))

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

        self.moveUp = QtWidgets.QPushButton("Move Up", self.list_button_widget)
        self.verticalLayout.addWidget(self.moveUp)
        self.moveUp.clicked.connect(self.listWidget.move_up)

        self.moveDown = QtWidgets.QPushButton("Move Down", self.list_button_widget)
        self.verticalLayout.addWidget(self.moveDown)
        self.moveDown.clicked.connect(self.listWidget.move_down)

        self.align = QtWidgets.QPushButton("Align", self.list_button_widget)
        self.verticalLayout.addWidget(self.align)
        self.align.clicked.connect(self.run_alignment)

        self.init_menu()

        self.text_width = 1000
        self.text_height = 300
        self.graph_text_diff = 36
        self.graphWidget = pg.PlotWidget(self.centralwidget)
        self.graphWidget.setGeometry(QtCore.QRect(self.pad, 2 * self.pad + self.list_height,
                                                  self.text_width + self.graph_text_diff,
                                                  self.text_height - 30))
        self.graphWidget.setBackground('w')


        self.text = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.text.setGeometry(QtCore.QRect(self.pad + self.graph_text_diff, 2 * self.pad + self.list_height,
                                           self.text_width, self.text_height))
        self.text.setFont(QtGui.QFont("Courier", 15))
        self.base_per_window = (self.text_width) / QtGui.QFontMetrics(QtGui.QFont("Courier", 15)).averageCharWidth()

        self.text.viewport().setAutoFillBackground(False)
        self.text.setFrameStyle(QtWidgets.QFrame.NoFrame)
        self.text.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        self.highlight = SnpHighlighter(self.text.document())

        self.graphWidget.setRange(xRange=[0, self.base_per_window], padding=0)
        self.graphWidget.setMouseEnabled(x=False)

        self.scroll = self.text.horizontalScrollBar()
        self.scroll.valueChanged.connect(self.scroll_update)

        self.text_button_widget = QtWidgets.QWidget(self.centralwidget)
        self.text_button_widget.setGeometry(QtCore.QRect(self.pad * 2 + 25 + self.text_width,
                                                         self.pad * 2 + self.list_height + 35,
                                                         self.button_width, self.text_height))

        self.verticalLayout2 = QtWidgets.QVBoxLayout(self.text_button_widget)
        self.verticalLayout2.setAlignment(QtCore.Qt.AlignTop)
        self.verticalLayout2.setSpacing(0)

        self.fileListWidget = FileList()

        self.setCentralWidget(self.centralwidget)

    def scroll_update(self):
        r = self.scroll.value()
        step_size = self.scroll.pageStep()
        l1 = r * self.base_per_window / step_size
        l2 = l1 + self.base_per_window
        self.graphWidget.setRange(xRange=[l1, l2], padding=0)

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
        self.ref_translation = list(self.seq_objs.values())[0].translate()
        self.ref_translate_aln = ""
        PROT_counter = 0
        DNA_counter = 0
        for ltr in self.alignment[0]:
            if ltr != '-' and PROT_counter < len(self.ref_translation):
                if (DNA_counter % 3) == 0:
                    self.ref_translate_aln += self.ref_translation[PROT_counter]
                    PROT_counter += 1
                    DNA_counter += 1

                else:
                    self.ref_translate_aln += " "
                    DNA_counter += 1
            else:
                self.ref_translate_aln += " "

        self.text.insertPlainText(self.ref_translate_aln + '\n')
        self.text.insertPlainText(str(self.alignment[0].seq) + '\n')
        for i, elem in enumerate(self.alignment[1:]):
            self.text.insertPlainText(str(elem.seq) + '\n')
            self.button_list.append(QtWidgets.QPushButton("View Trace"))
            self.button_list[i].clicked.connect(partial(self.show_trace, seq_id=self.listWidget.item(i+1).text(),
                                                        seq_idx=i + 1))

            self.verticalLayout2.addWidget(self.button_list[i])

    def show_trace(self, seq_id, seq_idx):
        self.graphWidget.clear()
        self.graphWidget.addLegend()
        seq_obj = self.seq_objs[seq_id]

        xticks = seq_obj.annotations['abif_raw']['PLOC2']
        xtick_labels = list(self.alignment[seq_idx])
        trace_idx = [0] + list(xticks) + [len(seq_obj.annotations['abif_raw'][channels[0]])]
        aln_idx = [i + 0.9 for i, ltr in enumerate(self.alignment[seq_idx]) if ltr != '-']
        aln_idx = [aln_idx[0] - 1] + aln_idx + [aln_idx[-1] + 1]
        print(trace_idx[-2], aln_idx[-2], len(seq_obj.seq))

        f = interp1d(trace_idx, aln_idx)
        new_x = f(np.arange(len(seq_obj.annotations['abif_raw'][channels[0]])))

        for i, channel in enumerate(channels):
            self.graphWidget.plot(new_x, seq_obj.annotations['abif_raw'][channel], name="GATC"[i], pen=(i,4))

        ax =self.graphWidget.getAxis('bottom')
        ax.setTicks([list(zip(np.arange(len(xtick_labels)) + 0.9, xtick_labels))])
        print("world")

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

    def move_up(self):
        currentRow = self.currentRow()
        if (self.ref) and currentRow == 1:
            return None
        currentItem = self.takeItem(currentRow)
        self.insertItem(currentRow - 1, currentItem)
        self.setCurrentRow(currentRow - 1)

    def move_down(self):
        currentRow = self.currentRow()
        if currentRow == self.count() - 1:
            return None
        currentItem = self.takeItem(currentRow)
        self.insertItem(currentRow + 1, currentItem)
        self.setCurrentRow(currentRow + 1)


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
                x = np.arange(len(tmp[channels[0]]))[::-1]
                tmp['PLOC2'] = [x[i] for i in tmp['PLOC2'][::-1]]
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

        if len(text) > 0 and text[0] == " ":
            return None

        ref = self.document.toPlainText().split('\n')[1]
        if len(ref) < len(text):
            diff = len(text) - len(ref)
            ref += diff * '-'

        idx = [i for i in range(len(text)) if text[i] != ref[i]]

        for i in idx:
            self.setFormat(i, 1, QtCore.Qt.red)

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