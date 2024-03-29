from PyQt5 import QtWidgets, QtCore, QtGui
import pyqtgraph as pg
import numpy as np
import sys, os, shutil, platform, csv, io, tempfile
from pathlib import Path
from functools import partial
from Bio import SeqIO, AlignIO, pairwise2
from Bio.Align.Applications import ClustalOmegaCommandline

format_dict = {'ab1': 'abi',
               'fasta': 'fasta',
               'seq': 'fasta'}

channels = ('DATA9', 'DATA10', 'DATA11', 'DATA12')

clustal_paths = {'Windows': Path('resources/base/clustal-omega-1.2.2-win64/clustalo.exe'),
                 'Darwin': Path('resources/base/clustal-omega-1.2.3-macosx'),
                 'Linux': Path('resources/base/clustalo-1.2.4-linux-x86_64')}

CLUSTAL_PATH = clustal_paths[platform.system()]
CODON_PATH = Path("resources/base/codon")

# TODO: Add support for other filetypes
#       Fix bug when updating/showing trace of sequences that have deleted nucleotides
#       Gray out sequences that are not related to the active trace


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Set widget size parameters
        self.text_height, self.text_width = 300, 1000
        self.list_height, self.list_width = 400, 400
        self.label_width = 300
        self.table_width = 526
        self.button_width = 90
        self.graph_x_diff, self.graph_y_diff = 36, 36
        self.pad = 10

        # Set character parameters
        self.base_per_window = self.text_width / QtGui.QFontMetrics(QtGui.QFont("Courier", 15)).averageCharWidth()
        self.letter_per_label = int(self.label_width / QtGui.QFontMetrics(QtGui.QFont("Courier", 15)).averageCharWidth())
        self.char_height = QtGui.QFontMetrics(QtGui.QFont("Courier", 15)).height()

        # Set window parameters
        self.window_width = 4 * self.pad + 2 * self.label_width + self.text_width + self.graph_x_diff
        self.window_height = 3 * self.pad + self.list_height + self.text_height + 100
        self.setWindowTitle("SDMCheck")
        self.resize(self.window_width, self.window_height)
        self.list_start = 2 * self.pad + self.label_width #(self.window_width - (self.list_width + self.button_width + 2 * self.pad + self.table_width)) / 2

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

        # Create Codon Widget
        self.codon_widget = QtWidgets.QWidget(self.centralwidget)
        self.codon_layout = QtWidgets.QVBoxLayout(self.codon_widget)
        self.codon_layout.setContentsMargins(0, 0, 0, 0)
        self.sp_select = QtWidgets.QComboBox()
        self.codon_layout.addWidget(self.sp_select)
        self.codon_widget.setGeometry(QtCore.QRect(self.list_start + self.list_width + self.pad * 2 + self.button_width,
                                                   self.pad, self.table_width, self.list_height))

        codon_files = [f for f in os.listdir(CODON_PATH)]
        for codon_file in codon_files:
            self.sp_select.addItem(codon_file[:-4])
            self.sp_select.activated[str].connect(self.set_codon_table)

        self.codon_model = QtGui.QStandardItemModel()
        self.codon_table = QtWidgets.QTableView()
        self.codon_table.horizontalHeader().hide()
        self.codon_table.verticalHeader().hide()
        self.codon_table.setFont(QtGui.QFont("Roboto", 9))

        self.codon_table.setModel(self.codon_model)
        self.codon_layout.addWidget(self.codon_table)

        # Create label widget for sequence alignment
        self.seq_label_widget = QtWidgets.QWidget(self.centralwidget)
        self.seq_label_layout = QtWidgets.QVBoxLayout(self.seq_label_widget)
        self.seq_label_layout.setAlignment(QtCore.Qt.AlignTop)
        self.seq_label_layout.setContentsMargins(0, 1, 0, 0)
        self.seq_label_layout.setSpacing(0)
        self.seq_label_widget.setGeometry(QtCore.QRect(self.pad, 2 * self.pad + self.list_height + self.graph_y_diff,
                                                       self.label_width, self.text_height))

        # Create Graph widget for plotting traces
        pg.setConfigOption('foreground', 'k')
        self.graphWidget = pg.PlotWidget(self.centralwidget)
        self.graphWidget.setBackground('darkGray')
        self.graphWidget.showAxis('top')
        self.bottom_ax = self.graphWidget.getAxis('bottom')
        self.bottom_ax.setLabel('Nucleotide')
        self.top_ax = self.graphWidget.getAxis('top')
        self.left_ax = self.graphWidget.getAxis('left')
        self.left_ax.setTicks([])
        self.left_ax.setWidth(35)

        self.graphWidget.setRange(xRange=[0, self.base_per_window], padding=0)
        self.graphWidget.setMouseEnabled(x=False)
        self.graphWidget.setGeometry(QtCore.QRect(2 * self.pad + self.label_width, 2 * self.pad + self.list_height,
                                                  self.text_width + self.graph_x_diff,
                                                  self.text_height + 6))

        self.legend = pg.LegendItem(offset=(1, 1), colCount=4)
        self.legend.setParentItem(self.bottom_ax)
        self.colors = 'gbyr'
        for i, n in enumerate("GATC"):
            item = pg.PlotDataItem(pen=(self.colors[i]))
            self.legend.addItem(item, n)

        # Create text widget superimposed on graph widget for sequence alignments
        self.text = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.text.setFont(QtGui.QFont("Courier", 15))
        color_palette = self.text.palette()
        color_palette.setColor(QtGui.QPalette.Text, QtCore.Qt.black)
        self.text.setPalette(color_palette)
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
        self.trace_button_layout.setSpacing(1)
        self.trace_button_layout.setContentsMargins(0, 0, 0, 0)
        self.text_button_widget.setGeometry(
            QtCore.QRect(self.pad * 3 + self.graph_x_diff + self.text_width + self.label_width,
                         self.pad * 2 + self.graph_y_diff + self.list_height + 2 * self.char_height + 5,
                         self.button_width, self.text_height))

        # Setup amino acid number increment widget
        self.aa_inc_widget = QtWidgets.QWidget(self.centralwidget)
        self.aa_inc_layout = QtWidgets.QHBoxLayout(self.aa_inc_widget)
        self.aa_inc_layout.setSpacing(0)
        self.aa_inc_layout.setContentsMargins(0, 0, 0, 0)
        self.aa_inc_widget.setGeometry(
            QtCore.QRect(self.pad * 3 + self.graph_x_diff + self.text_width + self.label_width,
                         self.pad * 2 + self.list_height,
                         self.button_width, 2 * self.char_height))

        # Set central widget
        self.setCentralWidget(self.centralwidget)

        # Initialize variables used by other methods
        self.button_list = []
        self.top_axis_ticks = []
        self.top_axis_values = []
        self.alignment, self.seq_objs, self.ref_translation, self.ref_translate_aln, self.active_trace = [None] * 5

    @property
    def labels(self):
        # Use only the last n letters of the file name that fit in the label widget layout
        start_val = self.letter_per_label

        labels = []
        for i in range(0, self.listWidget.count()):
            name = self.listWidget.item(i).text()
            labels.append(self.get_label(name))

        return labels

    def get_label(self, name):
        start_val = self.letter_per_label
        pre_text = '...'
        if len(name) < start_val + 4:
            start_val = len(name)
            pre_text = ''

        return pre_text + name[-start_val:-4]

    def set_codon_table(self, species):
        self.codon_model.clear()
        with open(CODON_PATH / (species + '.csv')) as codon_file:
            for row in csv.reader(codon_file):
                items = [ QtGui.QStandardItem(field) for field in row]
                self.codon_model.appendRow(items)

        self.codon_table.resizeColumnsToContents()
        self.codon_table.resizeRowsToContents()

    def scroll_update(self):
        """Update trace graph to match sequences position of text box using sync'd scrollbar"""
        r = self.scroll.value()

        bases = len(self.text.document().findBlockByLineNumber(0).text())
        n_windows = bases / self.base_per_window
        step = (r / self.scroll.maximum()) * (n_windows - 1)
        delta = self.base_per_window + 0.045
        l1 = step * delta
        l2 = l1 + delta
        self.graphWidget.setRange(xRange=[l1, l2], padding=0)

    def run_alignment(self):
        """Perform sequence alignment"""
        # Warn user if no sequences are loaded
        if not self.listWidget.count():
            warn_user("Warning: There are no sequences to align.")
            return None

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


        # Create label widget for filename and add to  label widget layout
        q_lab = QtWidgets.QLabel(self.centralwidget)
        q_lab.setText(self.get_label(name))
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
        self.top_axis_values = []
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
        self.top_axis_values = np.arange((len(self.top_axis_ticks))) + 1
        self.set_top_xticks()

    def adjust_xticks(self, value):
        self.top_axis_values += value
        self.set_top_xticks()

    def set_top_xticks(self):
        self.top_ax.setTicks([list(zip(self.top_axis_ticks, self.top_axis_values.astype(str)))])

    def show_alignment(self):
        """
         Print MSA to text box. Add all Labels and "Show Trace" buttons.
        """
        # Translate reference sequence
        self.translate_reference_sequence()
        self.add_aa_inc_buttons()

        # Add amino acid sequence to text box
        self.add_sequence_label('Amino Acid Reference Seq')
        self.text.insertPlainText(self.ref_translate_aln + '\n')

        # Add reference sequence to text box
        self.add_sequence_label('DNA Reference Seq')
        self.text.insertPlainText(str(self.alignment[0].seq))

        # Add all other sequences in MSA. +1 offset to skip the reference sequence
        for i, elem in enumerate(self.alignment[1:]):
            # Add sequence to text box
            self.add_sequence_label(self.listWidget.item(i + 1).text())
            self.text.insertPlainText('\n' + str(elem.seq))

            # Add buttons and connect to traces
            self.button_list.append(QtWidgets.QPushButton("View Trace"))
            self.trace_button_layout.addWidget(self.button_list[i])
            self.button_list[i].clicked.connect(partial(self.show_trace, seq_id=self.listWidget.item(i+1).text(),
                                                        seq_idx=i + 1))

        # Add button to hide traces
        self.hide_trace_button = QtWidgets.QPushButton("Hide Trace")
        self.trace_button_layout.addWidget(self.hide_trace_button)
        self.hide_trace_button.clicked.connect(self.hide_trace)

        # Pad remaining lines with whitespace
        for i in range(11 - self.text.document().blockCount()):
            self.text.insertPlainText("\n" + " " * len(elem.seq))

    def add_aa_inc_buttons(self):
        """add buttons to increment amino acid label number """
        minus_ten = QtWidgets.QPushButton("--", self.aa_inc_widget)
        minus_ten.clicked.connect(partial(self.adjust_xticks, value=-10))
        self.aa_inc_layout.addStretch()
        self.aa_inc_layout.addWidget(minus_ten)

        minus = QtWidgets.QPushButton("-", self.aa_inc_widget)
        minus.clicked.connect(partial(self.adjust_xticks, value=-1))
        self.aa_inc_layout.addStretch()
        self.aa_inc_layout.addWidget(minus)

        plus = QtWidgets.QPushButton("+", self.aa_inc_widget)
        plus.clicked.connect(partial(self.adjust_xticks, value=1))
        self.aa_inc_layout.addStretch()
        self.aa_inc_layout.addWidget(plus)

        plus_ten = QtWidgets.QPushButton("++", self.aa_inc_widget)
        plus_ten.clicked.connect(partial(self.adjust_xticks, value=10))
        self.aa_inc_layout.addStretch()
        self.aa_inc_layout.addWidget(plus_ten)


    def remove_alignment(self):
        """Clear textbox, graph widget, labels and "View Trace" buttons"""
        self.text.clear()
        self.hide_trace()


        if self.seq_label_layout.count():
            # Remove "View Trace" and "Hide Trace" buttons
            self.seq_label_layout.itemAt(0).widget().setParent(None)
            for i in reversed(range(self.trace_button_layout.count())):
                self.trace_button_layout.itemAt(i).widget().setParent(None)
                self.seq_label_layout.itemAt(i).widget().setParent(None)

            # Remove amino acid increment buttons
            for i in reversed(range(self.aa_inc_layout.count())):
                item = self.aa_inc_layout.takeAt(i)
                del item

        self.button_list = []

    def _reset(self):
        """Remove all items from list, all alignments, "Show trace" buttons and sequence lables"""
        self.remove_alignment()
        self.listWidget.clear()
        self.listWidget.ref = None

    def hide_trace(self):
        """Remove trace from graph"""
        self.graphWidget.clear()
        if self.active_trace is not None:
            self.active_trace.setText("View Trace")
            self.active_trace.setStyleSheet("QPushButton { color: black;}")
            self.bottom_ax.setTicks(None)

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

        xx = np.arange(len(seq_obj.annotations['abif_raw'][channels[0]]))
        new_x = np.interp(xx, trace_idx, aln_idx)

        # Plot traces corresponding to bases
        for i, channel in enumerate(channels):
            self.graphWidget.plot(new_x, seq_obj.annotations['abif_raw'][channel], name="GATC"[i], pen=self.colors[i])

        # Set x-axis labels to corresponding nucleotide
        self.bottom_ax.setTicks([list(zip(np.arange(len(x_tick_labels)) + 0.9, x_tick_labels))])

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

        # TODO: fix this to keep track of full path
        for i, label in enumerate(self.labels):

            j = popup_layout.count() + 1
            q_lab = QtWidgets.QLabel()
            q_lab.setText(label)
            q_lab.mypath = "Test"
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
        new_dir = Path(self.save_dir) / '_'.join([plasmid, protein])

        # Check for missing primer names
        for key in self.save_dict:
            if self.save_dict[key][1].text() == "":
                warn_user("You must provide a primer name for sequence " + self.save_dict[key][0])
                return None

        # Check for preexisting diriectories with the same name
        if Path(new_dir).exists():
            overwrite = prompt_user("{} already exists. Would you like to overwrite this record?".format(new_dir))
            if not overwrite:
                return None
            else:
                shutil.rmtree(new_dir)

        os.mkdir(new_dir)
        for key in self.save_dict:
            source = self.save_dict[key][0]
            ext = source.split(".")[-1]
            dest = new_dir / ('_'.join([plasmid, protein]) + '_' + '.'.join([self.save_dict[key][1].text(), ext]))
            shutil.copy(Path(source), dest)

        self.save_window.close()
        self.save_dict = {}

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

        # Remove elipsis if present
        if partial[:3] == '...':
            partial = partial[3:]

        # TODO: Fix this so that there can be only one
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
        if self.ref:
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
        if self.ref is None:

            # Warn user if no reference sequence is selected
            warn_user("Warning: No reference sequence has been selected. SDMCheck will try to guess the reference "
                      "sequence.")

            for i in range(1, self.count()):
                item = self.item(i)
                fname = item.text()

                if fname.endswith('.fasta'):
                    self.ref = self.item(i).text()
                    item.setForeground(QtCore.Qt.red)
                    self.insertItem(0, self.takeItem(i))
                    break
            else:
                warn_user('No candidate reference sequence found. Please use a .fasta file for reference sequences')
                return None



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

        # Write temporary FASTA file with all sequences
        temp_file = tempfile.NamedTemporaryFile('w+', delete=False)
        SeqIO.write(seq_objs.values(), temp_file, "fasta")
        temp_file.close()

        # Perform alignment
        cline = ClustalOmegaCommandline(str(CLUSTAL_PATH), infile=temp_file.name)
        stdout, stderr = cline()
        os.unlink(temp_file.name)

        # Read alignment using biopython
        alignment = io.StringIO(stdout)
        alignment = AlignIO.read(alignment, 'fasta')

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


def prompt_user(text, title="Overwrite record?"):
    msg = QtWidgets.QMessageBox()
    msg.setText(text)
    msg.setWindowTitle(title)
    msg.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
    returnValue = msg.exec()

    return returnValue == QtWidgets.QMessageBox.Yes


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
    QtWidgets.QApplication.setStyle('Fusion')
    app = QtWidgets.QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

