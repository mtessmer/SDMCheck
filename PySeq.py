from Bio import SeqIO
from PyQt5 import QtWidgets, QtCore, QtGui
import pyqtgraph as pg
import sys
from Bio import SeqIO


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.setWindowTitle("PySeq")
        self.resize(1200, 500)

        self.centralwidget = QtWidgets.QWidget(self)

        self.listWidget = FileList(self.centralwidget)
        self.listWidget.setGeometry(QtCore.QRect(10, 10, 500, 190))

        self.widget = QtWidgets.QWidget(self.centralwidget)
        self.widget.setGeometry(QtCore.QRect(510, 10, 90, 190))

        self.verticalLayout = QtWidgets.QVBoxLayout(self.widget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")

        self.add_file = QtWidgets.QPushButton("Add File", self.widget)
        self.verticalLayout.addWidget(self.add_file)
        self.add_file.clicked.connect(self.listWidget.add)

        self.remove_file = QtWidgets.QPushButton("Remove File", self.widget)
        self.verticalLayout.addWidget(self.remove_file)
        self.remove_file.clicked.connect(self.listWidget.remove)

        self.markRef = QtWidgets.QPushButton("Mark Reference", self.widget)
        self.verticalLayout.addWidget(self.markRef)
        self.markRef.clicked.connect(self.listWidget.mark_ref)

        self.align = QtWidgets.QPushButton("Align", self.widget)
        self.verticalLayout.addWidget(self.align)
        self.markRef.clicked.connect(self.listWidget.run_alignment)

        self.init_menu()

        self.graphWidget = MyPlotWidget(self.centralwidget)
        self.graphWidget.setGeometry(QtCore.QRect(10, 210, 1000, 200))
        self.graphWidget.setBackground('w')
        self.graphWidget.plotItem.vb.setLimits(yMin=0, yMax=1)
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

        for i in range(self.count()):
            self.item(i).setForeground(QtCore.Qt.black)

        self.ref = self.currentItem().text()
        i = self.currentRow()
        self.item(i).setForeground(QtCore.Qt.red)

    def run_alignment(self):
        pass


class MyPlotWidget(pg.PlotWidget):

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.MidButton:  # or Qt.MiddleButton
            self.__prevMousePos = event.pos()
        else:
            super(MyView, self).mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if event.buttons() == QtCore.Qt.MidButton:  # or Qt.MiddleButton
            offset = self.__prevMousePos - event.pos()
            self.__prevMousePos = event.pos()

            self.verticalScrollBar().setValue(self.verticalScrollBar().value() + offset.y())
            self.horizontalScrollBar().setValue(self.horizontalScrollBar().value() + offset.x())
        else:
            super(MyPlotWidget, self).mouseMoveEvent(event)


def main():
    app = QtWidgets.QApplication(sys.argv)
    main_window = MainWindow()

    main_window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()