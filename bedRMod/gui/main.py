from .controller import Controller


def start_gui():
    controller = Controller()
    controller.run()


if __name__ == "__main__":
    start_gui()
