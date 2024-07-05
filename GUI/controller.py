import csv
import os

import pandas as pd
import sympy
import sys

from tsv2bedRMod.tsv2bedRMod import csv2bedRMod
from tsv2bedRMod.helper import parse_excel_sheetnames


def funcify(expression):
    """
    Takes a string of an expression as an input and convert it into a python function.
    :return: function of passed expression string
    """
    x = sympy.symbols('x')
    expression = sympy.sympify(expression)
    func = sympy.lambdify(x, expression, "numpy")
    return func


class Controller:
    def __init__(self, ui):
        self.ui = ui

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

    def detect_file_type_delimiter(self, file):
        file_endings = (".odf", ".ods", ".odt", ".xlsx", ".xls", ".xlsb")
        if file.endswith(file_endings):
            self.sheetnames = parse_excel_sheetnames(file)
            self.columns = pd.read_excel(file, sheet_name=list(self.sheetnames)[0], nrows=0).columns.tolist()
            print(self.columns)

            return ".xlsx", None
        else:
            with open(file, 'r') as file:
                sample = file.read(1024)
                dialect = csv.Sniffer().sniff(sample)
                delimiter = dialect.delimiter
            self.delimiter = delimiter
            self.columns = pd.read_csv(file, nrows=0).columns.tolist()
            print(self.columns)
            return "csv", delimiter

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

    def convert2bedrmod(self):
        print(f"input file path: {self.ui.file_path.toPlainText()}")
        print(f"config yaml path: {self.ui.config_file_path.toPlainText()}")
        print(f"output file path: {self.ui.outfile_path.toPlainText()}")
        print(f"chrom column: {self.ui.ref_seg.toPlainText()}")
        print(f"position column: {self.ui.pos.toPlainText()}")
        print(f"0 indexed? {self.ui.index_0_button.isChecked()}")
        print(f"1 indexed? {self.ui.index_1_button.isChecked()}")
        print(f"modification info: {self.ui.modi.toPlainText()}")
        print(f"modification column? {self.ui.modi_button.isChecked()}")
        print(f"strand column {self.ui.strand.toPlainText()}")
        print(f"score column {self.ui.score.toPlainText()}")
        print(f"coverage column: {self.ui.coverage.toPlainText()}")
        print(f"frequency column: {self.ui.frequency.toPlainText()}")
        print(f"frequency function: {self.ui.frequency_function.toPlainText()}")
        print(f"score function: {self.ui.score_function.toPlainText()}")
        print(f"coverage function: {self.ui.coverage_function.toPlainText()}")
        print(f"delimiter: {self.ui.delimiter.checkedId()}")

        # as the input file can also be written directly into the field, check if it exists
        if not os.path.exists(self.ui.file_path.toPlainText()):
            print(f"The file at {self.ui.file_path.toPlainText()} does not exist! "
                  f"Please make sure you selected a valid file and try again.")
            return

        if not os.path.exists(self.ui.config_file_path.toPlainText()):
            print(f"The file at {self.ui.config_file_path.toPlainText()} does not exist! "
                  f"Please make sure you selected a valid file and try again.")
            return
        # how to handle outfile?

        # check delimiter of file.
        if self.delimiter is None:  # then its xlsx or other specified
            df = pd.read_excel(self.ui.file_path, self.selected_sheet)

        # pos
        # check if 0-index or 1-indexed
        if self.ui.index_1_button.isChecked():
            self.start_func = funcify("x - 1")

        # score
        score = self.ui.score.toPlainText()
        if score != self.score:
            print(score)
        if self.ui.score_function.toPlainText() != "Score function":
            self.score_func = funcify(self.ui.score_function.toPlainText())

        # strand
        strand = self.ui.strand.toPlainText()
        if strand == "+":
            self.strand = "+"
        elif strand == "-":
            self.strand = "-"
        elif strand != "strandedness":
            self.strand = strand
        # modi

        # coverage

        # frequency

        parsed_func = funcify("x * 2")
        print(parsed_func(3))

        csv2bedRMod(self.ui.file_path, self.ui.config_file_path, self.ui.outfile_path, self.ui.delimiter,
                    self.ui.ref_seg_column, self.ui.position_column, self.start_func, self.ui.modification_column,
                    self.ui.modi_button.isChecked(), self.ui.score_column, self.score_func, self.strand,
                    self.ui.coverage_column, self.coverage_func, self.ui.frequency_column, self.frequency_func)
