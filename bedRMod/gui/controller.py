import csv
import os
import pandas as pd

from PySide6 import QtCore
from PySide6.QtWidgets import QFileDialog, QLabel, QApplication, QMessageBox

from .view import bedRModWidget, NewConfigWindow, MainWindow

from bedRMod.convert2bedRMod import df2bedRMod
from bedRMod.helper import parse_excel_sheetnames, write_bioinformatics_keys, read_bioinformatics_keys, funcify


class Controller:
    def __init__(self):
        self.app = QApplication([])

        self.ui = bedRModWidget(self)
        self.window = MainWindow(self.ui)

        # set default values for Window
        self.columns = None
        self.sheetnames = None
        self.selected_sheet = None
        self.score = "score_column"

        self.delimiter = None

        self.strand = None
        self.start_func = None
        self.score_func = None
        self.coverage_func = None
        self.frequency_func = None

    def open_config_window(self):
        config_window = NewConfigWindow(self)
        config_window.show()

    @QtCore.Slot()
    def select_input_file(self):
        pathFile, ok = QFileDialog.getOpenFileName(self.ui,
                                                   "Open input file",
                                                   "",
                                                   "All Files(*)")
        if pathFile:
            self.ui.file_path.setText(pathFile)
            file_type, file_delimiter = self.detect_file_type_delimiter(pathFile)
            file_endings = (".odf", ".ods", ".odt", ".xlsx", ".xls", ".xlsb")

            if file_type in file_endings:
                self.ui.xlsx_file.setChecked(True)
                self.ui.custom_file_type.setChecked(False)
                self.ui.custom_file_delimiter.setEnabled(False)
                if self.ui.controller.sheetnames is not None:
                    self.ui.sheet_selector.addItems(self.sheetnames)
                    self.ui.sheet_info = QLabel("Select sheet")
                    self.ui.layout.addWidget(self.ui.sheet_info, 5, 0, 1, 1)
                    self.ui.layout.addWidget(self.ui.sheet_selector, 5, 1, 1, 3)
            else:
                self.ui.xlsx_file.setChecked(False)
                self.ui.custom_file_type.setChecked(True)
                if file_delimiter == "\t":
                    file_delimiter = "\\t"
                self.ui.custom_file_delimiter.setText(file_delimiter)
                self.ui.custom_file_delimiter.setEnabled(True)

            # set default path to output file
            input_file_split, ending = self.ui.file_path.toPlainText().rsplit('.', 1)
            output_file = input_file_split + '.bedrmod'
            self.ui.outfile_path.setText(output_file)

    @QtCore.Slot()
    def select_output_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        conf_file = QFileDialog()
        plain_file_path = self.ui.file_path.toPlainText()
        #print(plain_file_path)
        if plain_file_path:
            input_dir = os.path.dirname(plain_file_path)
            conf_file.setDirectory(input_dir)
            file_path, _ = conf_file.getSaveFileName(self.ui, "New .bedrmod", input_dir,
                                                     "BedRMod Files (*.bedrmod);;All Files (*)",
                                                     options=options)
        else:
            file_path, _ = conf_file.getSaveFileName(self.ui, "New .bedrmod", ".bedrmod",
                                                     "BedRMod Files (*.bedrmod);;All Files (*)",
                                                     options=options)
        if file_path:
            if not file_path.endswith(".bedrmod"):
                self.ui.outfile_path.setText(file_path + ".bedrmod")
            else:
                self.ui.outfile_path.setText(file_path)

    @QtCore.Slot()
    def select_config_file(self):
        pathFile, ok = QFileDialog.getOpenFileName(self.ui,
                                                   "Open the config file",
                                                   "",
                                                   "All Files(*)")
        if pathFile:
            self.ui.config_file_path.setText(pathFile)
            self.ui.controller.update_function_selection(pathFile)

    @QtCore.Slot()
    def create_new_file(self):
        # Ask the user to choose a location and name for the new file
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        conf_file = QFileDialog()
        file_path, _ = conf_file.getSaveFileName(self.ui, "New config.yaml", ".yaml",
                                                 "Config Files (*.yaml);;All Files (*)",
                                                 options=options)

        # If the user selected a file, create a new file
        if file_path:
            editor = NewConfigWindow(file_path)
            editor.show()
            # the editor shows the file contents but the file is not created yet.
            # this line trys to parse info from the file that does not exist!!
            self.update_function_selection(file_path)

    @QtCore.Slot()
    def on_index_button_toggled(self):
        if self.ui.index_0_button.isChecked():
            print(f"value index 0: {self.ui.index_0_button.isChecked()}")
        elif self.ui.index_1_button.isChecked():
            print(f"value index 1: {self.ui.index_1_button.isChecked()}")
        pass

    @QtCore.Slot()
    def on_delimiter_button_toggled(self):
        if self.ui.xlsx_file.isChecked():
            self.ui.custom_file_delimiter.setEnabled(False)
        elif self.ui.custom_file_type.isChecked():
            if self.ui.sheet_info is not None:
                self.ui.layout.removeWidget(self.ui.sheet_info)
            if self.ui.sheet_selector is not None:
                self.ui.layout.removeWidget(self.ui.sheet_selector)
            if self.ui.sheet_selector is not None:
                self.ui.sheet_selector.setParent(None)
            if self.ui.sheet_info is not None:
                self.ui.sheet_info.setParent(None)
            self.ui.custom_file_delimiter.setEnabled(True)

    @QtCore.Slot()
    def on_custom_modification_toggled(self):
        if self.ui.modi_button.isChecked():
            self.ui.layout.addWidget(self.ui.modi_custom, 8, 3, 1, 1)
            self.ui.modi_button.setChecked(True)
        else:
            self.ui.modi_button.setChecked(False)
            self.ui.layout.removeWidget(self.ui.modi_custom)
            self.ui.modi_custom.setParent(None)

    @QtCore.Slot()
    def on_custom_strand_toggled(self):
        if self.ui.strand_button.isChecked():
            self.ui.layout.addWidget(self.ui.strand_custom, 10, 3, 1, 1)
            self.ui.strand_button.setChecked(True)
        else:
            self.ui.strand_button.setChecked(False)
            self.ui.layout.removeWidget(self.ui.strand_custom)
            self.ui.strand_custom.setParent(None)

    def detect_file_type_delimiter(self, file):
        file_endings = (".odf", ".ods", ".odt", ".xlsx", ".xls", ".xlsb")
        if file.endswith(file_endings):
            self.sheetnames = parse_excel_sheetnames(file)
            self.columns = pd.read_excel(file, sheet_name=list(self.sheetnames)[0], nrows=0).columns.tolist()
            print(self.columns)

            return ".xlsx", None
        else:
            with open(file, 'r') as f:
                sample = f.read(1024)
                dialect = csv.Sniffer().sniff(sample)
                delimiter = dialect.delimiter
            self.delimiter = delimiter
            self.columns = pd.read_csv(file, nrows=0, delimiter=delimiter).columns.tolist()
            self.update_columns_selection()
            return "csv", delimiter

    def update_function_selection(self, file_path):
        """
        If a config file is selected and the functions are there an not none, they are displayed in their respective fields
        :param file_path:
        :return:
        """
        workflow, cov_f, func_f, score_f = read_bioinformatics_keys(file_path)
        if score_f:
            self.ui.score_function.setText(score_f)
        if cov_f:
            self.ui.coverage_function.setText(cov_f)
        if func_f:
            self.ui.frequency_function.setText(func_f)

    def update_columns_selection(self):
        self.ui.ref_seg.addItems(self.columns)
        self.ui.pos.addItems(self.columns)
        self.ui.modi.addItems(self.columns)
        self.ui.score.addItems(self.columns)
        self.ui.strand.addItems(self.columns)
        self.ui.coverage.addItems(self.columns)
        self.ui.frequency.addItems(self.columns)

    def on_sheet_selection(self):
        if self.ui.sheet_selector is not None:
            print(self.ui.sheet_selector.currentIndex())
            self.selected_sheet = self.ui.sheet_selector.currentIndex()
            self.columns = pd.read_excel(self.ui.file_path.toPlainText(), sheet_name=self.selected_sheet, nrows=0)\
                .columns.tolist()
            self.update_columns_selection()

    def done_message(self, success):
        msg = QMessageBox()

        if success == "Done":
            msg.setWindowTitle("Conversion Finished!")
            msg.setText(f"The conversion of {self.ui.file_path.toPlainText()} into {self.ui.outfile_path.toPlainText()} is finished!")
        elif success == "Error":
            msg.setWindowTitle("Error!")
            msg.setText(f"An Error occured during the the conversion of {self.ui.file_path.toPlainText()} into "
                        f"{self.ui.outfile_path.toPlainText()}. Check the input values and try again.")
        close_button = msg.addButton("Close", QMessageBox.ButtonRole.AcceptRole)

        close_button.clicked.connect(msg.accept)
        msg.exec()

    def convert2bedrmod(self):
        # as the input file can also be written directly into the field, check if it exists
        input_file = self.ui.file_path.toPlainText()
        if not os.path.exists(input_file):
            print(f"The file at {input_file} does not exist! "
                  f"Please make sure you selected a valid input file and try again.")
            return

        config_file = self.ui.config_file_path.toPlainText()
        if not os.path.exists(config_file):
            print(f"The file at {config_file} does not exist! "
                  f"Please make sure you selected a valid config file and try again.")
            return

        output_file = self.ui.outfile_path.toPlainText()

        # check delimiter of file.
        if self.delimiter is None:  # then its xlsx or other specified
            df = pd.read_excel(input_file, self.selected_sheet)
        else:
            df = pd.read_csv(input_file, delimiter=self.delimiter)
        # ref_seg
        ref_seg = self.ui.ref_seg.currentText()

        # pos
        pos = self.ui.pos.currentText()
        # check if 0-index or 1-indexed
        if self.ui.index_1_button.isChecked():
            self.start_func = funcify("lambda x: x - 1")

        # score
        score = self.ui.score.currentText()
        # if score != self.score:
        #     print(score)
        if self.ui.score_function.toPlainText() != "Score function":
            self.score_func = self.ui.score_function.toPlainText()
            write_bioinformatics_keys(config_file, score_function=self.ui.score_function.toPlainText())

        # strand
        strand = None
        if self.ui.strand_button.isChecked():
            strand = self.ui.strand_custom.toPlainText()
        else:
            strand = self.ui.strand.currentText()
        # modi
        modi = None
        modi_button_check = None  # modi button check is exactly the other way around because df2bedRMod() ist
        if self.ui.modi_button.isChecked():
            modi = self.ui.modi_custom.toPlainText()
            modi_button_check = False
        else:
            modi = self.ui.modi.currentText()
            modi_button_check = True

        # coverage
        cov = self.ui.coverage.currentText()
        if self.ui.coverage_function.toPlainText() != "Coverage function":
            self.coverage_func = self.ui.coverage_function.toPlainText()
            write_bioinformatics_keys(config_file, coverage_function=self.coverage_func)
        # frequency
        freq = self.ui.frequency.currentText()
        if self.ui.frequency_function.toPlainText() != "Frequency function":
            self.frequency_func = self.ui.frequency_function.toPlainText()
            write_bioinformatics_keys(config_file, frequency_function=self.frequency_func)

        if self.score_func:
            self.score_func = funcify("lambda score: " + self.score_func)
        if self.coverage_func:
            self.coverage_func = funcify("lambda coverage: " + self.coverage_func)
        if self.frequency_func:
            self.frequency_func = funcify("lambda frequency: " + self.frequency_func)

        print(f"input file path: {self.ui.file_path.toPlainText()}")
        print(f"config yaml path: {self.ui.config_file_path.toPlainText()}")
        print(f"output file path: {self.ui.outfile_path.toPlainText()}")
        print(f"chrom column: {self.ui.ref_seg.currentText()}")
        print(f"position column: {self.ui.pos.currentText()}")
        print(f"0 indexed? {self.ui.index_0_button.isChecked()}")
        print(f"1 indexed? {self.ui.index_1_button.isChecked()}")
        print(f"modification info: {self.ui.modi.currentText()}")
        print(f"custom modification? {self.ui.modi_button.isChecked()}")
        print(f"strand column: {self.ui.strand.currentText()}")
        print(f"score column: {self.ui.score.currentText()}")
        print(f"coverage column: {self.ui.coverage.currentText()}")
        print(f"frequency column: {self.ui.frequency.currentText()}")
        print(f"frequency function: {self.ui.frequency_function.toPlainText()}")
        print(f"score function: {self.ui.score_function.toPlainText()}")
        print(f"coverage function: {self.ui.coverage_function.toPlainText()}")
        print(f"delimiter: {self.delimiter}")

        function_success = df2bedRMod(df, config_file, output_file, ref_seg=ref_seg, start=pos,
                                      start_function=self.start_func, modi=modi,
                   modi_column=modi_button_check, score=score, score_function=self.score_func, strand=strand,
                   coverage=cov, coverage_function=self.coverage_func, frequency=freq,
                   frequency_function=self.frequency_func, from_gui=True)

        self.done_message(function_success)

    def run(self):
        self.window.show()
        self.app.exec()