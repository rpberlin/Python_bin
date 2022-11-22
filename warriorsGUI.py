#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QApplication, QWidget, QLineEdit, QLabel, QSpinBox, QSlider
import warriorsNameMaker as wN
import sys

def onClick():
    fullName = wN.getName()
    text_edit.setPlainText(fullName)
    #text_edit.setPlainText("Your Warrior Cat Description Goes Here")

def display():
    print(text_edit.toPlainText())

def main():

    app = QApplication(sys.argv)
    window = QWidget()
    window.setGeometry(400,400,300,300)
    window.setWindowTitle("Warrior Cats")

    text_edit = QtWidgets.QTextEdit(window)
    text_edit.setPlaceholderText("Your Character Will Be Described Here")



    button = QtWidgets.QPushButton(window)
    button.setText("Create Your Cat")
    button.clicked.connect(onClick)
    button.move(75,200)

    #widget.setPixmap(QPixmap('otje.jpg'))

    #slider = QtWidgets.QSlider(window)
    #slider.setText('This is a Slider')


    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
