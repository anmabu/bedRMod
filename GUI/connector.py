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
        print(f"input file path: {self.ui.input_file_path}")
        print(f"config yaml path: {self.ui.config_yaml_path}")
        print(f"output file path: {self.ui.output_file_path}")
        print(f"chrom column: {self.ui.ref_seg_column}")
        print(f"position column: {self.ui.position_column}")
        print(f"0 indexed? {self.ui.index_0_button.isChecked()}")
        print(f"1 indexed? {self.ui.index_1_button.isChecked()}")
        print(f"modification info: {self.ui.modification_column}")
        print(f"modification column? {self.ui.modi_button.isChecked()}")
        print(f"strand column {self.ui.strand_column}")
        print(f"score column {self.ui.score_column}")
        print(f"coverage column: {self.ui.coverage_column}")
        print(f"frequency column: {self.ui.frequency_column}")
        print(f"frequency function: {self.ui.frequency_function.toPlainText()}")
        print(f"score function: {self.ui.score_function.toPlainText()}")
        print(f"coverage function: {self.ui.coverage_function.toPlainText()}")
        self.ui.start_func = None
        self.ui.score_func = None
        self.ui.coverage_func = None
        self.ui.frequency_func = None

        # use sympy to parse functions from strings into lambda functions

        parsed_func = funcify("x * 2")
        print(parsed_func(3))
