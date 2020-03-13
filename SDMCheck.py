from PyQt5 import QtWidgets, QtCore, QtGui
import pyqtgraph as pg
import numpy as np
from scipy.interpolate import interp1d
import sys, os, shutil, platform
from functools import partial
from Bio import SeqIO, AlignIO, pairwise2
from Bio.Align.Applications import ClustalOmegaCommandline

format_dict = {'ab1': 'abi',
               'fasta': 'fasta',
               'seq': 'fasta'}

channels = ('DATA9', 'DATA10', 'DATA11', 'DATA12')

clustal_paths = {'Windows': os.path.normpath('extra\clustal-omega-1.2.2-win64\clustalo.exe'),
                 'Darwin': os.path.normpath('extra\clustal-omega-1.2.3-macosx'),
                 'Linux': os.path.normpath('extra\clustalo-1.2.4-linux-x86_64')}

CLUSTAL_PATH = clustal_paths[platform.system()]


# TODO: Add support for other filetypes
#       Implement Model|View design pattern
#       Shifter to adjust Amino Acid nunber
#       Codon Table With E. coli codon freq


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Set widget size parameters
        self.text_height, self.text_width = 300, 1000
        self.list_height, self.list_width = 200, 500
        self.codon_height, self.codon_width = 200, 200
        self.label_width = 300
        self.button_width = 90
        self.graph_x_diff, self.graph_y_diff = 36, 36
        self.pad = 10

        # Set character parameters
        self.base_per_window = self.text_width / QtGui.QFontMetrics(QtGui.QFont("Courier", 15)).averageCharWidth()
        self.letter_per_label = int(self.label_width / QtGui.QFontMetrics(QtGui.QFont("Courier", 15)).averageCharWidth())
        self.char_height = QtGui.QFontMetrics(QtGui.QFont("Courier", 15)).height()

        # Set window parameters
        self.window_width = 4 * self.pad + self.label_width + self.text_width + self.graph_x_diff + self.button_width
        self.window_height = 4 * self.pad + self.list_height + self.text_height + self.codon_height
        self.setWindowTitle("SDMCheck")
        self.resize(self.window_width, self.window_height)
        self.list_start = (self.window_width - (self.list_width + self.button_width + self.pad)) / 2

        # Initialize file menu
        self.init_menu()

        # Create central widget
        self.centralwidget = QtWidgets.QWidget(self)

        # Create file list widget
        self.listWidget = FileList(self.centralwidget)
        self.listWidget.setGeometry(QtCore.QRect(self.list_start, self.pad, self.list_width, self.list_height))

        # Create list button widget
        self.list_button_widget = QtWidgets.QWidget(self.centralwidget)
        self.file_button_layout = QtWidgets.QVBoxLayout(self.list_button_widget)
        self.file_button_layout.setContentsMargins(0, 0, 0, 0)
        self.list_button_widget.setGeometry(QtCore.QRect(self.list_start + self.list_width + self.pad,
                                                         self.pad, self.button_width, self.list_height))

        # Define buttons in list button widget
        self.add_file = QtWidgets.QPushButton("Add File", self.list_button_widget)
        self.file_button_layout.addWidget(self.add_file)
        self.add_file.clicked.connect(self.listWidget.add)

        self.remove_file = QtWidgets.QPushButton("Remove File", self.list_button_widget)
        self.file_button_layout.addWidget(self.remove_file)
        self.remove_file.clicked.connect(self.listWidget.remove)

        self.markRef = QtWidgets.QPushButton("Mark Reference", self.list_button_widget)
        self.file_button_layout.addWidget(self.markRef)
        self.markRef.clicked.connect(self.listWidget.mark_ref)

        self.moveUp = QtWidgets.QPushButton("Move Up", self.list_button_widget)
        self.file_button_layout.addWidget(self.moveUp)
        self.moveUp.clicked.connect(self.listWidget.move_up)

        self.moveDown = QtWidgets.QPushButton("Move Down", self.list_button_widget)
        self.file_button_layout.addWidget(self.moveDown)
        self.moveDown.clicked.connect(self.listWidget.move_down)

        self.align = QtWidgets.QPushButton("Align", self.list_button_widget)
        self.file_button_layout.addWidget(self.align)
        self.align.clicked.connect(self.run_alignment)

        self.reset = QtWidgets.QPushButton("Reset", self.list_button_widget)
        self.file_button_layout.addWidget(self.reset)
        self.reset.clicked.connect(self._reset)

        # Create label widget for sequence alignment
        self.seq_label_widget = QtWidgets.QWidget(self.centralwidget)
        self.seq_label_layout = QtWidgets.QVBoxLayout(self.seq_label_widget)
        self.seq_label_layout.setAlignment(QtCore.Qt.AlignTop)
        self.seq_label_layout.setContentsMargins(0,0.5,0,0)
        self.seq_label_layout.setSpacing(0)
        self.seq_label_widget.setGeometry(QtCore.QRect(self.pad, 2 * self.pad + self.list_height + self.graph_y_diff,
                                                       self.label_width, self.text_height))

        # Create Graph widget for plotting traces
        self.graphWidget = pg.PlotWidget(self.centralwidget)
        self.graphWidget.setBackground('w')
        self.graphWidget.showAxis('top')
        self.graphWidget.setRange(xRange=[0, self.base_per_window], padding=0)
        self.graphWidget.setMouseEnabled(x=False)
        self.graphWidget.setGeometry(QtCore.QRect(2 * self.pad + self.label_width, 2 * self.pad + self.list_height,
                                                  self.text_width + self.graph_x_diff,
                                                  self.text_height + 6))

        # Create text widget superimposed on graph widget for sequence alignments
        self.text = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.text.setFont(QtGui.QFont("Courier", 15))
        self.text.viewport().setAutoFillBackground(False)
        self.text.setFrameStyle(QtWidgets.QFrame.NoFrame)
        self.text.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        self.highlight = SnpHighlighter(self.text.document())
        self.text.setGeometry(QtCore.QRect(2 * self.pad  + self.label_width + self.graph_x_diff,
                                           2 * self.pad + self.list_height + self.graph_y_diff,
                                           self.text_width, self.text_height))

        # Link graph widget and text widget scrollbar
        self.scroll = self.text.horizontalScrollBar()
        self.scroll.valueChanged.connect(self.scroll_update)

        # Create button widget for "show trace" buttons created after sequence alignment
        self.text_button_widget = QtWidgets.QWidget(self.centralwidget)
        self.trace_button_layout = QtWidgets.QVBoxLayout(self.text_button_widget)
        self.trace_button_layout.setAlignment(QtCore.Qt.AlignTop)
        self.trace_button_layout.setSpacing(0)
        self.trace_button_layout.setContentsMargins(0,0,0,0)
        self.text_button_widget.setGeometry(
            QtCore.QRect(self.pad * 3 + self.graph_x_diff + self.text_width + self.label_width,
                         self.pad * 2 + self.graph_y_diff + self.list_height + 2 * self.char_height,
                         self.button_width, self.text_height))

        # Set central widget
        self.setCentralWidget(self.centralwidget)

        # Initialize variables used by other methods
        self.labels = []
        self.button_list = []
        self.alignment, self.seq_objs, self.ref_translation, self.ref_translate_aln, self.active_trace = [None] * 5

    def scroll_update(self):
        """
        Update trace graph to match sequences position of text box using sync'd scrollbar
        """
        r = self.scroll.value()
        step_size = self.scroll.pageStep()
        l1 = r * self.base_per_window / step_size
        l2 = l1 + self.base_per_window
        self.graphWidget.setRange(xRange=[l1, l2], padding=0)


    def run_alignment(self):
        """
        Perform sequence alignment
        """
        # Warn user if no sequences are loaded
        if not self.listWidget.count():
            warn_user("Warning: There are no sequences to align.")

        # Warn user if no reference sequence is selected
        if not self.listWidget.ref:
            warn_user("Warning: You cannot perform alignment without first marking a reference sequence.")

        # Remove current alignmnet if it exists
        if self.button_list:
            self.remove_alignment()

        self.alignment, self.seq_objs = self.listWidget.run_alignment()
        self.show_alignment()

    def add_sequence_label(self, name):
        """
        Write sequence file name next to sequence in text box
        :param name: str
            Name of sequence file
        """
        # Use only the last n letters of the file name that fit in the label widget layout
        name = os.path.basename(name)
        start_val = self.letter_per_label
        pre_text = '...'
        if len(name) < start_val + 4:
            start_val = len(name)
            pre_text = ''

        # Create label widget for filename and add to  label widget layout
        self.labels.append(pre_text + name[-start_val:-4])
        q_lab = QtWidgets.QLabel(self.centralwidget)
        q_lab.setText(self.labels[-1])
        q_lab.setFont(QtGui.QFont("Courier", 15))
        q_lab.setAlignment(QtCore.Qt.AlignRight)
        self.seq_label_layout.addWidget(q_lab)

    def translate_reference_sequence(self):
        """
        Translate reference sequence and align amino acid to first nucleotide of the corresponding codon
        """
        # Translate sequence with Biopython
        self.ref_translation = list(self.seq_objs.values())[0].translate()

        # Perform alignment of AA Seq to Nucleotide
        self.ref_translate_aln = ""
        self.top_axis_ticks = []
        prot_counter = 0
        DNA_counter = 0
        for i, ltr in enumerate(self.alignment[0]):
            if ltr != '-' and prot_counter < len(self.ref_translation):
                if (DNA_counter % 3) == 0:
                    self.top_axis_ticks.append(i+0.9)
                    self.ref_translate_aln += self.ref_translation[prot_counter]
                    prot_counter += 1
                    DNA_counter += 1
                else:
                    self.ref_translate_aln += " "
                    DNA_counter += 1
            else:
                self.ref_translate_aln += " "

        # Show amino acid numbers
        self.ax_top = self.graphWidget.getAxis('top')
        self.ax_top.setTicks([list(zip(self.top_axis_ticks, np.arange(len(self.top_axis_ticks)) + 1))])

    def show_alignment(self):
        """
         Print MSA to text box. Add all Labels and "Show Trace" buttons.
        """
        # Translate reference sequence
        self.translate_reference_sequence()

        # Add amino acid sequence to text box
        self.add_sequence_label('Amino Acid Reference Seq')
        self.text.insertPlainText(self.ref_translate_aln + '\n')

        # Add reference sequence to text box
        self.add_sequence_label('DNA Reference Seq')
        self.text.insertPlainText(str(self.alignment[0].seq) + '\n')

        # Add all other sequences in MSA. +1 offset to skip the reference sequence
        for i, elem in enumerate(self.alignment[1:]):
            # Add sequence to text box
            self.add_sequence_label(self.listWidget.item(i + 1).text())
            self.text.insertPlainText(str(elem.seq) + '\n')

            # Add buttons and connect to traces
            self.button_list.append(QtWidgets.QPushButton("View Trace"))
            self.trace_button_layout.addWidget(self.button_list[i])
            self.button_list[i].clicked.connect(partial(self.show_trace, seq_id=self.listWidget.item(i+1).text(),
                                                        seq_idx=i + 1))

        # Add button to hide traces
        self.hide_trace_button = QtWidgets.QPushButton("Hide Trace")
        self.trace_button_layout.addWidget(self.hide_trace_button)
        self.hide_trace_button.clicked.connect(self.hide_trace)

    def remove_alignment(self):
        """Clear textbox, graph widget, labels and "View Trace" buttons"""
        self.text.clear()
        self.hide_trace()

        if self.seq_label_layout.count():
            self.seq_label_layout.itemAt(0).widget().setParent(None)
            for i in reversed(range(self.trace_button_layout.count())):
                self.trace_button_layout.itemAt(i).widget().setParent(None)
                self.seq_label_layout.itemAt(i).widget().setParent(None)

        self.button_list = []

    def _reset(self):
        """Remove all items from list, all alignments, "Show trace" buttons and sequence lables"""
        self.remove_alignment()
        self.listWidget.clear()

    def hide_trace(self):
        """Remove trace from graph"""
        self.graphWidget.clear()
        if self.active_trace is not None:
            self.active_trace.setText("View Trace")
            self.active_trace.setStyleSheet("QPushButton { color: black;}")
            self.ax.setTicks([list(zip(np.arange(0, len(self.alignment[0]), 10),
                                       np.arange(0, len(self.alignment[0]), 10)))])


    def show_trace(self, seq_id, seq_idx):
        """
        Show trace of `seq_id` sequence corresponding to the `seq_idx` index

        :param seq_id: str
            file name of sequence for plotting
        :param seq_idx:
            index of sequence in alignment
        """

        # Hide traces if any are plotted and reset "Show Trace"button
        if self.active_trace is not None:
            self.hide_trace()
            self.active_trace.setText("View Trace")
            self.active_trace.setStyleSheet("QPushButton { color: black;}")

        # Set active trace
        self.active_trace = self.button_list[seq_idx - 1]

        # Align traces to the sequence as it appears in the text box
        doc = self.text.document()
        seq_obj = self.seq_objs[seq_id]
        x_ticks = seq_obj.annotations['abif_raw']['PLOC2']
        x_tick_labels = list(doc.findBlockByLineNumber(seq_idx + 1).text())
        trace_idx = [0] + list(x_ticks) + [len(seq_obj.annotations['abif_raw'][channels[0]])]
        aln_idx = [i + 0.9 for i, ltr in enumerate(doc.findBlockByLineNumber(seq_idx + 1).text()) if ltr != '-']
        aln_idx = [aln_idx[0] - 1] + aln_idx + [aln_idx[-1] + 1]
        f = interp1d(trace_idx, aln_idx)
        new_x = f(np.arange(len(seq_obj.annotations['abif_raw'][channels[0]])))

        # Plot traces corresponding to bases
        for i, channel in enumerate(channels):
            self.graphWidget.plot(new_x, seq_obj.annotations['abif_raw'][channel], name="GATC"[i], pen=(i,4))

        # Set x-axis labels to corresponding nucleotide
        self.ax = self.graphWidget.getAxis('bottom')
        self.ax.setTicks([list(zip(np.arange(len(x_tick_labels)) + 0.9, x_tick_labels))])

        # Change "Show Trace" button to "Update" and mark red to show it is the active trace
        self.active_trace.setText("Update")
        self.active_trace.setStyleSheet("QPushButton { color: red;}")

    def init_menu(self):
        """
        Initialize file menue
        """
        # Create menu bar
        menu_bar = self.menuBar()

        # Create actions and bind to functions
        file_menu = menu_bar.addMenu('File')
        save_action = QtWidgets.QAction('Save', self)
        save_action.triggered.connect(self.save)
        quit_action = QtWidgets.QAction('Quit', self)
        quit_action.triggered.connect(self.quit_app)

        # Add actions to menu
        file_menu.addAction(save_action)
        file_menu.addAction(quit_action)

    def set_save_dir(self, line_edit):
        self.save_dir = QtWidgets.QFileDialog().getExistingDirectory()
        line_edit.setText(self.save_dir)

    def save(self):
        """Protocol for saving sequencing results in a standardized strain database format"""
        self.save_window = QtWidgets.QWidget()
        self.save_window.setWindowTitle("Save Files")
        self.save_dict = {}
        popup_layout = QtWidgets.QGridLayout(self.save_window)

        directory_label = QtWidgets.QLabel()
        directory_label.setText("Directory:")
        directory_input = QtWidgets.QLineEdit()
        directory_button = QtWidgets.QPushButton()
        directory_button.setText("Set Directory")
        directory_button.clicked.connect(partial(self.set_save_dir, line_edit=directory_input))
        popup_layout.addWidget(directory_label, 0, 0)
        popup_layout.addWidget(directory_input, 0, 1, 1, 2)
        popup_layout.addWidget(directory_button, 0, 3)

        plasmid_label = QtWidgets.QLabel()
        plasmid_label.setText("Plasmid")
        plasmid_input = QtWidgets.QLineEdit()
        popup_layout.addWidget(plasmid_label, 1, 0)
        popup_layout.addWidget(plasmid_input, 1, 1, 1, 2)

        protein_label = QtWidgets.QLabel()
        protein_label.setText("Protein")
        protein_input = QtWidgets.QLineEdit()
        popup_layout.addWidget(protein_label, 2, 0)
        popup_layout.addWidget(protein_input, 2, 1, 1, 2)
        le = []
        for i, label in enumerate(self.labels[2:]):

            j = popup_layout.count() + 1
            q_lab = QtWidgets.QLabel()
            q_lab.setText(label)
            chkbox = QtWidgets.QCheckBox()
            le.append(QtWidgets.QLineEdit())
            le[i].setPlaceholderText('Primer')
            chkbox.stateChanged.connect(partial(self.check_label_event, checkbox_widget=chkbox,
                                                label_widget=q_lab, le=le[i]))

            popup_layout.addWidget(q_lab, j, 0)
            popup_layout.addWidget(chkbox, j, 1)
            popup_layout.addWidget(le[i])

        save_button = QtWidgets.QPushButton()
        save_button.setText("Save")
        save_button.clicked.connect(partial(self._save, plasmid=plasmid_input,
                                            protein=protein_input))

        popup_layout.addWidget(save_button, popup_layout.count() + 1, 3)
        self.save_window.show()

    def check_label_event(self, checkbox_widget, label_widget, le):
        if checkbox_widget.isChecked():
            source = self.listWidget.get_full_path(partial=label_widget.text())
            self.save_dict[label_widget.text()] = (source, le)
        elif label_widget.text() in self.save_dict:
            del self.save_dict[label_widget.text()]
        else:
            return None

    def _save(self, plasmid, protein):
        plasmid = plasmid.text()
        protein = protein.text()
        protein = "_".join(protein.split(" "))
        new_dir = os.path.join(self.save_dir, '_'.join([plasmid, protein]))
        os.makedirs(new_dir)

        for key in self.save_dict:
            source = self.save_dict[key][0]
            ext = source.split(".")[-1]
            dest = os.path.join(new_dir, '_'.join([plasmid, protein]) + '_' + '.'.join([self.save_dict[key][1].text(), ext]))
            shutil.copy(source, dest)

        self.save_window.close()


    def quit_app(self):
        """Exit the application"""
        QtWidgets.qApp.quit()


class FileList(QtWidgets.QListWidget):
    """Subclass QListWidget to add features to list"""

    def __init__(self, parent=None):
        super(FileList, self).__init__(parent)
        self.ref = None
        self.setAcceptDrops(True)

    def dragEnterEvent(self, event):
        """Allow file objects to be dragged into listWidget"""
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event):
        """Show copy icon and follow mouse if dragged object has a file path"""
        if event.mimeData().hasUrls:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        """Add file to list when dropped"""
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
        """Use file dialog to add files to list"""
        files, _ = QtWidgets.QFileDialog.getOpenFileNames()
        self.addItems(files)

    def remove(self):
        """Remove file from list"""
        for item in self.selectedItems():
            self.takeItem(self.row(item))

    def get_full_path(self, partial):
        full_path = [self.item(i).text() for i in range(self.count()) if partial in self.item(i).text()][0]
        return full_path

    def move_up(self):
        """Move file up one in list"""
        currentRow = self.currentRow()
        if (self.ref) and currentRow == 1:
            return None
        currentItem = self.takeItem(currentRow)
        self.insertItem(currentRow - 1, currentItem)
        self.setCurrentRow(currentRow - 1)

    def move_down(self):
        """Move file down one in list"""
        currentRow = self.currentRow()
        if currentRow == self.count() - 1:
            return None
        currentItem = self.takeItem(currentRow)
        self.insertItem(currentRow + 1, currentItem)
        self.setCurrentRow(currentRow + 1)

    def mark_ref(self):
        """Assign reference sequence file and move reference to top"""
        # Do nothing if: no files in list, no files highlighted or highlighted file is currently marked as ref
        if (not self.count()) or (self.currentRow() == -1) or (self.ref == self.currentItem().text()):
            return None

        # Remove reference marker if it is on a different file
        if not self.ref:
            self.item(0).setForeground(QtCore.Qt.black)

        # Store ref file name, color red and move to top
        self.ref = self.currentItem().text()
        self.item(self.currentRow()).setForeground(QtCore.Qt.red)
        self.insertItem(0, self.takeItem(self.currentRow()))

    def run_alignment(self):
        """
        Read in all sequence information, test and fix reverse compliment strands, remove sequences with poor homology
        :return alignment, seq_obj: Bio.Align.MultipleSequenceAlignment, dict[str] = Bio.SeqRecord
            Returns multiple sequence alignment (MSA) and sequence records for each qualifying sequence in the alignment
        """
        # Read reference sequence information into dictionary
        seq_objs = {self.ref: read_file(self.ref)}

        # Temporarily store sequence of reference as string
        ref_seq = str(seq_objs[self.ref].seq)

        # Read sequence information of all other files into `seq_objs` dict
        for i in range(1, self.count()):
            # Read sequence information
            test_seq = self.item(i).text()
            seq_objs[test_seq] = read_file(test_seq)

            # If alignment to ref is less than 70% identical test reverse compliment
            taln = pairwise2.align.localxx(ref_seq, str(seq_objs[test_seq].seq), score_only=True)
            if (taln / len(ref_seq)) < 0.7:

                # Store annotation information
                tmp = seq_objs[test_seq].annotations['abif_raw']

                # Reverse compliment sequence
                seq_objs[test_seq] = seq_objs[test_seq].reverse_complement(id=seq_objs[test_seq].id + '_rc')

                # If reverse compliment aligns with less than 70% identity remove the sequence
                taln = pairwise2.align.localxx(ref_seq, str(seq_objs[test_seq].seq), score_only=True)
                if (taln / len(ref_seq)) < 0.7:
                    warn_user(self.item(i).text() + ' has low sequence homology and will not be used in the alignments')
                    self.takeItem(i)
                    del seq_objs[test_seq]

                # Reverse compliment trace
                tmp[channels[0]], tmp[channels[1]], tmp[channels[2]], tmp[channels[3]] = \
                    tmp[channels[3]][::-1], tmp[channels[2]][::-1], tmp[channels[1]][::-1], tmp[channels[0]][::-1]
                x = np.arange(len(tmp[channels[0]]))[::-1]
                tmp['PLOC2'] = [x[i] for i in tmp['PLOC2'][::-1]]
                seq_objs[test_seq].annotations['abif_raw'] = tmp.copy()

        # Write FASTA file with all sequences
        SeqIO.write(seq_objs.values(), 'seqs.fasta', "fasta")

        # Remove old alignment file if it exists
        if os.path.exists('seqs.aln'):
            os.remove('seqs.aln')

        # Perform alignment
        clustal_omega_exe = os.path.join(os.path.dirname(__file__), CLUSTAL_PATH)
        cline = ClustalOmegaCommandline(clustal_omega_exe, infile="seqs.fasta", outfile="seqs.aln")
        stdout, stderr = cline()

        # Read alignment file
        alignment = AlignIO.read('seqs.aln', 'fasta')

        return alignment, seq_objs


class SnpHighlighter(QtGui.QSyntaxHighlighter):
    """Subclass QtGui.QSyntaxHighlighter to create highlighting rules for single nucleotide polymorphisms (SNPs)"""

    def __init__(self, document):
        QtGui.QSyntaxHighlighter.__init__(self, document)
        self.document = document

    def highlightBlock(self, text):
        """Highlight block of `text`
        :param text, str
            block of text to be highlighted
        """
        # If there are no lines do nothing
        if not self.document.toPlainText():
            return None

        # If line begins with a space then ignore it (No nucleotides)
        if len(text) > 0 and text[0] == " ":
            return None

        # Assign reference sequence to the second line of text in the editor
        ref = self.document.toPlainText().split('\n')[1]

        # Add dummy nucleotides to reference sequence if `text` is longer
        if len(ref) < len(text):
            diff = len(text) - len(ref)
            ref += diff * '-'

        # Get indices of `text` nucleotides that do not match reference sequence
        idx = [i for i in range(len(text)) if text[i] != ref[i]]

        # Highlight mismatches with reference sequence as red text
        for i in idx:
            self.setFormat(i, 1, QtCore.Qt.red)


def read_file(file_name):
    """Helper function to read different sequence file types
    :param file_name, str
        name of file to read
    """
    # get file extension and call appropriate Biopython subfunction
    ext = file_name.split('.')[-1]
    file_format = format_dict[ext]
    return SeqIO.read(file_name, file_format)


def warn_user(text):
    """
    Trigger popup with warning text
    :param text: str

    :return:
    """
    msg = QtWidgets.QMessageBox()
    msg.setWindowTitle("Warning")
    msg.setText(text)
    msg.setIcon(QtWidgets.QMessageBox.Warning)
    msg.exec_()



def main():
    """Main Qt Application"""
    app = QtWidgets.QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
