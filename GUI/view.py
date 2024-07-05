import os
import sys
import yaml

from PySide6 import QtCore, QtWidgets
from PySide6.QtWidgets import QLabel, QLineEdit, QFileDialog, QPushButton, QTextEdit, QFrame, QRadioButton, \
    QWidget, QVBoxLayout, QButtonGroup, QComboBox

from controller import Controller


class NewConfigWindow(QWidget):
    def __init__(self, file_path):
        super().__init__()
        NewConfigWindow.setWindowTitle(self, f"{file_path}")
        self.file_path = file_path
        self.text_edit = QTextEdit()
        config = yaml.safe_load(open("../test/test_config.yaml", "r"))
        formatted_yaml = yaml.dump(config, default_flow_style=False, sort_keys=False)
        self.text_edit.setText(formatted_yaml)
        self.initUI()

    def initUI(self):
        save_button = QPushButton('Save and Close', self)
        save_button.clicked.connect(self.save_and_close)

        layout = QVBoxLayout()
        layout.addWidget(self.text_edit)
        layout.addWidget(save_button)

        self.setLayout(layout)

    def save_and_close(self):
        content = self.text_edit.toPlainText()
        if not self.file_path.endswith(".yaml"):
            self.file_path += ".yaml"
        with open(self.file_path, 'w') as new_file:
            new_file.write(content)
        self.close()


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.controller = Controller(self)

        self.file_path = None
        self.input_file = None
        self.config_file_path = None
        self.config_file = None
        self.new_config_file = None
        self.outfile_path = None
        self.output_file = None
        self.delimiter = None
        self.xlsx_file = None
        self.sheet_info = None
        self.sheet_selector = None
        self.custom_file_type = None
        self.custom_file_delimiter = None
        self.ref_seg = None
        self.pos = None
        self.index_0_button = None
        self.index_1_button = None
        self.button_group = None
        self.modi = None
        self.modi_button = None
        self.score = None
        self.score_function = None
        self.strand = None
        self.coverage = None
        self.coverage_function = None
        self.frequency = None
        self.frequency_function = None
        self.convert = None
        self.info_text = None

        self.initUI()

    def initUI(self):
        MainWindow.setWindowTitle(self, "Convert to bedRMod")

        self.info_text = QTextEdit("some info what to do here, lorem ipsum dolor et amit")
        self.info_text.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        font_metrics = self.info_text.fontMetrics()
        line_height = font_metrics.lineSpacing()

        # Set the height of the QTextEdit to the height of one line
        self.info_text.setFixedHeight(line_height * 1.8)
        self.info_text.isReadOnly()

        # input file
        input_label = QLabel("Select input file:")
        self.file_path = QTextEdit()
        self.file_path.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.file_path.setText("none selected")
        self.file_path.setFixedHeight(line_height * 1.6)
        # self.file_path.setStyleSheet("background-color: white")
        self.input_file = QPushButton("...")
        self.input_file.clicked.connect(self.select_input_file)

        # config file
        config_label = QLabel("Select config file:")
        self.config_file_path = QTextEdit()
        self.config_file_path.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.config_file_path.setText("none selected")
        self.config_file_path.setFixedHeight(line_height * 1.6)
        # self.config_file_path.setStyleSheet("background-color: white")
        self.config_file = QPushButton("...")
        self.config_file.clicked.connect(self.select_config_file)
        self.new_config_file = QPushButton("New Config file")
        self.new_config_file.clicked.connect(self.create_new_file)

        # output file
        output_label = QLabel("Select output file path:")
        self.outfile_path = QTextEdit()
        self.outfile_path.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.outfile_path.setText("none selected")
        self.outfile_path.setFixedHeight(line_height * 1.6)
        # self.file_path.setStyleSheet("background-color: white")
        self.output_file = QPushButton("...")
        self.output_file.clicked.connect(self.select_output_file)

        # delimiter info
        delimiter_label = QLabel("Select file type / column delimiter")
        self.delimiter = QButtonGroup(self)
        self.xlsx_file = QRadioButton(".xlsx")
        self.delimiter.addButton(self.xlsx_file, 1)

        self.custom_file_type = QRadioButton("custom delimiter")
        self.delimiter.addButton(self.custom_file_type, 2)
        self.custom_file_delimiter = QLineEdit()
        self.custom_file_delimiter.setText("e.g ',', '\\t'")
        self.custom_file_delimiter.setEnabled(False)  # only enable if custom delimiter is selected

        self.custom_file_type.setChecked(True)
        self.xlsx_file.toggled.connect(self.on_delimiter_button_toggled)
        self.custom_file_type.toggled.connect(self.on_delimiter_button_toggled)

        self.sheet_selector = QComboBox()
        self.sheet_selector.currentIndexChanged.connect(self.controller.on_sheet_selection)

        # ref_seg
        ref_seg_label = QLabel("Reference Segment / Chromosome")
        ref_seg_label.setToolTip("Select column containing reference segment information. "
                                 "One reference segment per row in the file.")
        self.ref_seg = QTextEdit()
        self.ref_seg.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.ref_seg.setText('chrom')
        self.ref_seg.setFixedHeight(line_height * 1.6)

        # pos
        pos_label = QLabel("Position")
        pos_label.setToolTip("Select column containing position of modification. "
                             "Only a single position per row in this column.")
        self.pos = QTextEdit()
        self.pos.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.pos.setText('start_col')
        self.pos.setFixedHeight(line_height * 1.6)
        self.index_0_button = QRadioButton("0 indexed data")
        self.index_1_button = QRadioButton("1 indexed data")
        self.button_group = QButtonGroup()
        self.button_group.addButton(self.index_0_button)
        self.button_group.addButton(self.index_1_button)
        self.index_0_button.setChecked(True)
        # self.index_0_button.toggled.connect(self.onIndexButtonToggled)
        # self.index_1_button.toggled.connect(self.onIndexButtonToggled)

        # modification type
        modi_label = QLabel("Modification type / column")
        modi_label.setToolTip("Select the column that contains the modifications or input the modomics shortname "
                              "for the modification type when only one is present in the file.")
        self.modi = QTextEdit()
        self.modi.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.modi.setText('name')
        self.modi.setFixedHeight(line_height * 1.6)
        self.modi_button = QRadioButton("Column?")

        # score
        score_label = QLabel("Score")
        score_label.setToolTip("Select column containing modification score in the range of [0; 1000]."
                               "If the score in this interval is not readily available, a function to convert the "
                               "given values can be passed."
                               "Also a single integer can be passed as a fixed score value for the whole file.")
        self.score = QTextEdit()
        self.score.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.score.setText(self.controller.score)
        self.score.setFixedHeight(line_height * 1.6)
        self.score_function = QTextEdit()
        self.score_function.setText("Score function")
        self.score_function.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.score_function.setToolTip("When writing a function and refering to a column name in the calculation "
                                       "(e.g. FDR), please refer to this column name as row['FDR']. "
                                       "(Or do this calculation in a script and store the result in the same file)")
        self.score_function.setFixedHeight(line_height * 1.6)

        # strand
        strand_label = QLabel("Strandedness / strand column")
        strand_label.setToolTip("Select the column that contains the strand information. If strandedness is the same "
                                "for the whole file, '+' or '-' will work, too.")
        self.strand = QTextEdit()
        self.strand.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.strand.setText('strandedness')
        self.strand.setFixedHeight(line_height * 1.6)

        # coverage
        coverage_label = QLabel("Coverage")
        coverage_label.setToolTip("Select the column that contains the coverage information.")
        self.coverage = QTextEdit()
        self.coverage.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.coverage.setText('coverage_col')
        self.coverage.setFixedHeight(line_height * 1.6)
        self.coverage_function = QTextEdit()
        self.coverage_function.setText("Coverage function")
        self.coverage_function.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.coverage_function.setToolTip("When writing a function and refering to a column name in the calculation "
                                          "(e.g. cov), please refer to this column name as 'cov'. "
                                          "(Or do this calculation in a script and store the result in the same file)")
        self.coverage_function.setFixedHeight(line_height * 1.6)

        # frequency
        frequency_label = QLabel("Modification Frequency")
        frequency_label.setToolTip("Select the column that contains the modification frequency information."
                                   "If the modification frequency is not stored but can be calculated, pass"
                                   "the function to calculate it in the last field. ")
        self.frequency = QTextEdit()
        self.frequency.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.frequency.setText('frequency_col')
        self.frequency.setFixedHeight(line_height * 1.6)
        self.frequency_function = QTextEdit()
        self.frequency_function.setText("Frequency function")
        self.frequency_function.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.frequency_function.setToolTip("When writing a function and refering to a column name in the calculation "
                                           "(e.g. FDR), please refer to this column name as row['FDR']. "
                                           "(Or do this calculation in a script and store the result in the same file)")
        self.frequency_function.setFixedHeight(line_height * 1.6)

        # convert now!
        self.convert = QPushButton()
        self.convert.setText("Convert!")
        self.convert.clicked.connect(self.controller.convert2bedrmod)

        # layout stuff
        self.layout = QtWidgets.QGridLayout()
        self.layout.setColumnStretch(0, 1)
        self.layout.setColumnStretch(1, 3)
        self.layout.setColumnStretch(2, 1)
        self.layout.setColumnStretch(3, 1)

        # input file
        self.layout.addWidget(input_label, 1, 0)
        self.layout.addWidget(self.file_path, 1, 1)
        self.layout.addWidget(self.input_file, 1, 2, 1, 2)

        # config file
        self.layout.addWidget(config_label, 2, 0, 1, 1)
        self.layout.addWidget(self.config_file_path, 2, 1, 1, 1)
        self.layout.addWidget(self.config_file, 2, 2, 1, 1)
        self.layout.addWidget(self.new_config_file, 2, 3, 1, 1)

        # output file
        self.layout.addWidget(output_label, 3, 0)
        self.layout.addWidget(self.outfile_path, 3, 1)
        self.layout.addWidget(self.output_file, 3, 2, 1, 2)

        # delimiter stuff
        self.layout.addWidget(delimiter_label, 4, 0)
        # layout.addWidget(self.delimiter, 4, 1)
        self.layout.addWidget(self.xlsx_file, 4, 1, 1, 1)
        self.layout.addWidget(self.custom_file_type, 4, 2, 1, 1)
        self.layout.addWidget(self.custom_file_delimiter, 4, 3, 1, 1)

        # conversion info
        # self.layout.addWidget(self.info_text, 5, 0, 1, 4)

        # chrom column
        self.layout.addWidget(ref_seg_label, 6, 0, 1, 1)
        self.layout.addWidget(self.ref_seg, 6, 1, 1, 1)

        # start column
        self.layout.addWidget(pos_label, 7, 0, 1, 1)
        self.layout.addWidget(self.pos, 7, 1, 1, 1)
        self.layout.addWidget(self.index_0_button, 7, 2, 1, 1)
        self.layout.addWidget(self.index_1_button, 7, 3, 1, 1)

        # modification label
        self.layout.addWidget(modi_label, 8, 0, 1, 1)
        self.layout.addWidget(self.modi, 8, 1, 1, 1)
        self.layout.addWidget(self.modi_button, 8, 2, 1, 1)

        # score column
        self.layout.addWidget(score_label, 9, 0, 1, 1)
        self.layout.addWidget(self.score, 9, 1, 1, 1)
        self.layout.addWidget(self.score_function, 9, 2, 1, 2)

        # strand info
        self.layout.addWidget(strand_label, 10, 0, 1, 1)
        self.layout.addWidget(self.strand, 10, 1, 1, 1)

        # coverage info
        self.layout.addWidget(coverage_label, 11, 0, 1, 1)
        self.layout.addWidget(self.coverage, 11, 1, 1, 1)
        self.layout.addWidget(self.coverage_function, 11, 2, 1, 2)

        # frequency info
        self.layout.addWidget(frequency_label, 12, 0, 1, 1)
        self.layout.addWidget(self.frequency, 12, 1, 1, 1)
        self.layout.addWidget(self.frequency_function, 12, 2, 1, 2)

        # convert button
        self.layout.addWidget(self.convert, 13, 0, 1, 4)

        self.setLayout(self.layout)

    @QtCore.Slot()
    def select_input_file(self):
        pathFile, ok = QFileDialog.getOpenFileName(self,
                                                   "Open input file",
                                                   "",
                                                   "All Files(*)")
        if pathFile:
            self.file_path.setText(pathFile)
            file_type, file_delimiter = self.controller.detect_file_type_delimiter(pathFile)
            file_endings = (".odf", ".ods", ".odt", ".xlsx", ".xls", ".xlsb")

            if file_type in file_endings:
                self.xlsx_file.setChecked(True)
                self.custom_file_type.setChecked(False)
                self.custom_file_delimiter.setEnabled(False)
                if self.controller.sheetnames is not None:
                    self.sheet_selector.addItems(self.controller.sheetnames)
                    self.sheet_info = QLabel("Select sheet")
                    self.layout.addWidget(self.sheet_info, 5, 0, 1, 1)
                    self.layout.addWidget(self.sheet_selector, 5, 1, 1, 3)
            else:
                self.xlsx_file.setChecked(False)
                self.custom_file_type.setChecked(True)
                if file_delimiter == "\t":
                    file_delimiter = "\\t"
                self.custom_file_delimiter.setText(file_delimiter)
                self.custom_file_delimiter.setEnabled(True)

    @QtCore.Slot()
    def select_output_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        conf_file = QFileDialog()
        plain_file_path = self.file_path.toPlainText()
        print(plain_file_path)
        if plain_file_path:
            input_dir = os.path.dirname(plain_file_path)
            conf_file.setDirectory(input_dir)
            # print(input_dir)
            file_path, _ = conf_file.getSaveFileName(self, "New .bedrmod", input_dir,
                                                     "BedRMod Files (*.bedrmod);;All Files (*)",
                                                     options=options)
        else:
            file_path, _ = conf_file.getSaveFileName(self, "New .bedrmod", ".bedrmod",
                                                     "BedRMod Files (*.bedrmod);;All Files (*)",
                                                     options=options)
        if file_path:
            if not file_path.endswith(".bedrmod"):
                self.outfile_path.setText(file_path + ".bedrmod")
            else:
                self.outfile_path.setText(file_path)

    @QtCore.Slot()
    def select_config_file(self):
        pathFile, ok = QFileDialog.getOpenFileName(self,
                                                   "Open the config file",
                                                   "",
                                                   "All Files(*)")
        if pathFile:
            self.config_file_path.setText(pathFile)

    @QtCore.Slot()
    def create_new_file(self):
        # Ask the user to choose a location and name for the new file
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        conf_file = QFileDialog()
        file_path, _ = conf_file.getSaveFileName(self, "New config.yaml", ".yaml",
                                                 "Config Files (*.yaml);;All Files (*)",
                                                 options=options)

        # file_path = conf_file.exec_()
        # If the user selected a file, create a new file
        if file_path:
            self.editor = NewConfigWindow(file_path)
            self.editor.show()

    @QtCore.Slot()
    def onIndexButtonToggled(self):
        if self.index_0_button.isChecked():
            print(f"value index 0: {self.index_0_button.isChecked()}")
        elif self.index_1_button.isChecked():
            print(f"value index 1: {self.index_1_button.isChecked()}")
        pass

    @QtCore.Slot()
    def on_delimiter_button_toggled(self):
        if self.xlsx_file.isChecked():
            self.custom_file_delimiter.setEnabled(False)
        elif self.custom_file_type.isChecked():
            self.layout.removeWidget(self.sheet_info)
            self.layout.removeWidget(self.sheet_selector)
            self.sheet_selector.setParent(None)
            self.sheet_info.setParent(None)
            self.custom_file_delimiter.setEnabled(True)




def start_gui():
    app = QtWidgets.QApplication([])

    widget = MainWindow()
    widget.show()

    sys.exit(app.exec())
