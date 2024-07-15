from PySide6.QtWidgets import QApplication
from controller import Controller


def start_gui():
    app = QApplication([])
    controller = Controller(app)
    app.exec()


if __name__ == "__main__":
    start_gui()
