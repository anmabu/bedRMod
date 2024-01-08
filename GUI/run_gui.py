import sys

from PySide6 import QtCore, QtWidgets, QtGui
from PySide6.QtWidgets import QLabel, QLineEdit, QFileDialog, QPushButton, QComboBox, QTextEdit, QFrame
from PySide6.QtCore import Qt


class MyWidget(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()

        input_label = QLabel("Select input file:")
        self.file_path = QLabel()
        self.file_path.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.file_path.setText("none selected")
        self.file_path.setStyleSheet("background-color: white")
        self.input_file = QPushButton("...")
        self.input_file.clicked.connect(self.select_input_file)

        config_label = QLabel("Select config file:")
        self.config_file_path = QLabel()
        self.config_file_path.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.config_file_path.setText("none selected")
        self.config_file_path.setStyleSheet("background-color: white")
        self.config_file = QPushButton("...")
        self.config_file.clicked.connect(self.select_config_file)

        delimiter_label = QLabel("Select file type / column delimiter")
        self.delimiter = QComboBox()
        self.delimiter.addItem("comma")
        self.delimiter.addItem("tab")
        self.delimiter.addItem("xlsx")

        info_text = QTextEdit("some info what to do here, lorem ipsum dolor et amit")
        info_text.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        font_metrics = info_text.fontMetrics()
        line_height = font_metrics.lineSpacing()

        # Set the height of the QTextEdit to the height of one line
        info_text.setFixedHeight(line_height * 1.8)
        info_text.isReadOnly()
        # info_text.setVerticalScrollBarPolicy(0)

        # ref_seg
        ref_seg_label = QLabel("Reference Segment / Chromosome")
        ref_seg_label.setToolTip("Select column containing reference segment information. "
                                 "One reference segment per row in the file.")
        self.ref_seg = QTextEdit()
        self.ref_seg.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.ref_seg.setText('"ref_seg"')
        self.ref_seg.setFixedHeight(line_height * 1.6)

        # pos
        pos_label = QLabel("Position")
        pos_label.setToolTip("Select column containing position of modification. "
                             "Only a single position per row in this column.")
        self.pos = QTextEdit()
        self.pos.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.pos.setText('"pos"')
        self.pos.setFixedHeight(line_height * 1.6)

        # modification type
        modi_label = QLabel("Modification type / column")
        modi_label.setToolTip("Select the column that contains the modifications or input the modomics shortname "
                              "for the modification type when only one is present in the file.")
        self.modi = QTextEdit()
        self.modi.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.modi.setText('"m1A"')
        self.modi.setFixedHeight(line_height * 1.6)

        # score
        score_label = QLabel("Score")
        score_label.setToolTip("Select column containing modification score in the range of [0; 1000]."
                               "If the score in this interval is not readily available, a function to convert the "
                               "given values can be passed."
                               "Also a single integer can be passed as a fixed score value for the whole file.")
        self.score = QTextEdit()
        self.score.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.score.setText('"pos"')
        self.score.setFixedHeight(line_height * 1.6)
        self.score_function = QTextEdit()
        self.score_function.setText("Score function")
        self.score_function.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.score_function.setFixedHeight(line_height * 1.6)

        # layout stuff
        layout = QtWidgets.QGridLayout()

        # input file
        layout.addWidget(input_label, 1, 0)
        layout.addWidget(self.file_path, 1, 1)
        layout.addWidget(self.input_file, 1, 2)

        # config file
        layout.addWidget(config_label, 2, 0)
        layout.addWidget(self.config_file_path, 2, 1)
        layout.addWidget(self.config_file, 2, 2)

        layout.addWidget(delimiter_label, 3, 0)
        layout.addWidget(self.delimiter, 3, 1)

        layout.addWidget(info_text, 4, 0, 1, 3)

        layout.addWidget(ref_seg_label, 5, 0, 1, 1)
        layout.addWidget(self.ref_seg, 5, 1, 1, 1)

        layout.addWidget(pos_label, 6, 0, 1, 1)
        layout.addWidget(self.pos, 6, 1, 1, 1)

        layout.addWidget(modi_label, 7, 0, 1, 1)
        layout.addWidget(self.modi, 7, 1, 1, 1)

        layout.addWidget(score_label, 8, 0, 1, 1)
        layout.addWidget(self.score, 8, 1, 1, 1)
        layout.addWidget(self.score_function, 8, 2, 1, 1)

        self.setLayout(layout)

    @QtCore.Slot()
    def select_input_file(self):
        pathFile, ok = QFileDialog.getOpenFileName(self,
                                                   "Open the input file",
                                                   "",
                                                   "All Files(*)")
        if pathFile:
            self.file_path.setText(pathFile)

    @QtCore.Slot()
    def select_config_file(self):
        pathFile, ok = QFileDialog.getOpenFileName(self,
                                                   "Open the config file",
                                                   "",
                                                   "All Files(*)")
        if pathFile:
            self.config_file_path.setText(pathFile)


if __name__ == "__main__":
    app = QtWidgets.QApplication([])

    widget = MyWidget()
    # widget.resize(800, 600)
    widget.show()

    sys.exit(app.exec())
