import sys

from PySide6 import QtCore, QtWidgets, QtGui
from PySide6.QtWidgets import QLabel, QLineEdit, QFileDialog, QPushButton, QComboBox, QTextEdit, QFrame, QRadioButton, \
    QWidget, QVBoxLayout


class EditorWindow(QWidget):
    def __init__(self, file_path):
        super().__init__()
        EditorWindow.setWindowTitle(self, f"{file_path}")
        self.file_path = file_path
        self.text_edit = QTextEdit()
        self.text_edit.setText(
            'options:\n  modification_type: "RNA"\n  organism: 9606\n  assembly: "GRCh38"\n  annotation_source: null\n  '
            'annotation_version: null\n  sequencing_platform: "Illumina NovaSeq 6000"\n  basecalling: null\n  '
            'bioinformatics_workflow: "https://github.com/y9c/m6A-sacseq"\n  experiment: "2-50 ng of poly-A enriched or '
            'ribosome RNA-depleted RNAs were fragmented and ligated, then divided in a 2:1 ratio. 2/3 of the starting '
            'materials are labeled by MjDim1, while the remaining 1/3 serve as the untreated control. After reverse '
            'transcription with HIV reverse transcriptase (Worthington Biochemical Corp), the cyclic allyl m6A sites '
            'are converted to mismatches, while unconverted m6A sites in the control group are read as A."\n  '
            'external_source: "GEO;GSE198246"\n  methods: "m6A-SAC-seq"\n  references:\n    pubmed_id: "36434097"\n \n  '
            'modifications_file: "example_files/mod_indices.csv"')
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
        with open(self.file_path, 'w') as new_file:
            new_file.write(content)
        self.close()


class MyWidget(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        MyWidget.setWindowTitle(self, "Convert to bedRMod")

        info_text = QTextEdit("some info what to do here, lorem ipsum dolor et amit")
        info_text.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        font_metrics = info_text.fontMetrics()
        line_height = font_metrics.lineSpacing()

        # Set the height of the QTextEdit to the height of one line
        info_text.setFixedHeight(line_height * 1.8)
        info_text.isReadOnly()

        input_label = QLabel("Select input file:")
        self.file_path = QTextEdit()
        self.file_path.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.file_path.setText("none selected")
        self.file_path.setFixedHeight(line_height * 1.6)
        # self.file_path.setStyleSheet("background-color: white")
        self.input_file = QPushButton("...")
        self.input_file.clicked.connect(self.select_input_file)

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

        delimiter_label = QLabel("Select file type / column delimiter")
        self.delimiter = QComboBox()
        self.delimiter.addItem("comma")
        self.delimiter.addItem("tab")
        self.delimiter.addItem("xlsx")

        # ref_seg
        ref_seg_label = QLabel("Reference Segment / Chromosome")
        ref_seg_label.setToolTip("Select column containing reference segment information. "
                                 "One reference segment per row in the file.")
        self.ref_seg = QTextEdit()
        self.ref_seg.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.ref_seg.setText('ref_seg')
        self.ref_seg.setFixedHeight(line_height * 1.6)

        # pos
        pos_label = QLabel("Position")
        pos_label.setToolTip("Select column containing position of modification. "
                             "Only a single position per row in this column.")
        self.pos = QTextEdit()
        self.pos.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.pos.setText('pos')
        self.pos.setFixedHeight(line_height * 1.6)

        # modification type
        modi_label = QLabel("Modification type / column")
        modi_label.setToolTip("Select the column that contains the modifications or input the modomics shortname "
                              "for the modification type when only one is present in the file.")
        self.modi = QTextEdit()
        self.modi.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.modi.setText('modification')
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
        self.score.setText("score")
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
                                "for the whole file, '+' or '-' will work, too." )
        self.strand = QTextEdit()
        self.strand.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.strand.setText('strand')
        self.strand.setFixedHeight(line_height * 1.6)

        # coverage
        coverage_label = QLabel("Coverage")
        coverage_label.setToolTip("Select the column that contains the coverage information.")
        self.coverage = QTextEdit()
        self.coverage.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.coverage.setText('coverage')
        self.coverage.setFixedHeight(line_height * 1.6)

        # frequency
        frequency_label = QLabel("Modification Frequency")
        frequency_label.setToolTip("Select the column that contains the modification frequency information."
                                   "If the modification frequency is not stored but can be calculated, pass"
                                   "the function to calculate it in the last field. ")
        self.frequency = QTextEdit()
        self.frequency.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.frequency.setText('modification frequency')
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

        # layout stuff
        layout = QtWidgets.QGridLayout()

        # input file
        layout.addWidget(input_label, 1, 0)
        layout.addWidget(self.file_path, 1, 1)
        layout.addWidget(self.input_file, 1, 2)

        # config file
        layout.addWidget(config_label, 2, 0, 1, 1)
        layout.addWidget(self.config_file_path, 2, 1, 1, 1)
        layout.addWidget(self.config_file, 2, 2, 1, 1)
        layout.addWidget(self.new_config_file, 2, 3, 1, 1)

        layout.addWidget(delimiter_label, 3, 0)
        layout.addWidget(self.delimiter, 3, 1)

        layout.addWidget(info_text, 4, 0, 1, 3)

        layout.addWidget(ref_seg_label, 5, 0, 1, 1)
        layout.addWidget(self.ref_seg, 5, 1, 1, 1)

        layout.addWidget(pos_label, 6, 0, 1, 1)
        layout.addWidget(self.pos, 6, 1, 1, 1)

        layout.addWidget(modi_label, 7, 0, 1, 1)
        layout.addWidget(self.modi, 7, 1, 1, 1)
        layout.addWidget(self.modi_button, 7, 2, 1, 1)

        layout.addWidget(score_label, 8, 0, 1, 1)
        layout.addWidget(self.score, 8, 1, 1, 1)
        layout.addWidget(self.score_function, 8, 2, 1, 1)

        layout.addWidget(strand_label, 9, 0, 1, 1)
        layout.addWidget(self.strand, 9, 1, 1, 1)

        layout.addWidget(coverage_label, 10, 0, 1, 1)
        layout.addWidget(self.coverage, 10, 1, 1, 1)

        layout.addWidget(frequency_label, 11, 0, 1, 1)
        layout.addWidget(self.frequency, 11, 1, 1, 1)
        layout.addWidget(self.frequency_function, 11, 2, 1, 1)

        layout.addWidget(self.convert, 12, 0, 1, 4)

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

    @QtCore.Slot()
    def create_new_file(self):
        # Ask the user to choose a location and name for the new file
        conf_file = QFileDialog()
        conf_file.setWindowTitle("New config.yaml")
        conf_file.setDefaultSuffix("yaml")
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_path, _ = conf_file.getSaveFileName(self, "Save File", ".yaml", "Config Files (*.yaml);;All Files (*)",
                                                 options=options)
        # file_path = conf_file.exec_()
        # If the user selected a file, create a new file
        if file_path:
            self.editor = EditorWindow(file_path)
            self.editor.show()


if __name__ == "__main__":
    app = QtWidgets.QApplication([])

    widget = MyWidget()
    # widget.resize(800, 600)
    widget.show()

    sys.exit(app.exec())
