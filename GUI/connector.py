import os
import sympy
import sys

from tsv2bedRMod.tsv2bedRMod import csv2bedRMod


def funcify(expression):
    """
    Takes a string of an expression as an input and convert it into a python function.
    :return: function of passed expression string
    """
    x = sympy.symbols('x')
    expression = sympy.sympify(expression)
    func = sympy.lambdify(x, expression, "numpy")
    return func


class Connector:
    def __init__(self, ui):
        self.ui = ui

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
        print(f"delimiter: {self.ui.delimiter.currentData()}")
        self.ui.start_func = None
        self.ui.score_func = None
        self.ui.coverage_func = None
        self.ui.frequency_func = None

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

        parsed_func = funcify("x * 2")
        print(parsed_func(3))

        csv2bedRMod(self.ui.file_path, self.ui.config_file_path, self.ui.outfile_path, self.ui.delimiter,
                    self.ui.ref_seg_column, self.ui.position_column, self.ui.start_func, self.ui.modification_column,
                    self.ui.modi_button.isChecked(), self.ui.score_column, self.ui.score_func, self.ui.strand_column,
                    self.ui.coverage_column, self.ui.coverage_func, self.ui.frequency_column, self.ui.frequency_func)
