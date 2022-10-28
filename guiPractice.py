#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QApplication, QWidget, QLineEdit, QLabel, QSpinBox, QSlider
import sys

def onClick():
    text_edit.setPlainText("Your Warrior Cat Description Goes Here")

def display():
    print(text_edit.toPlainText())

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

#slider = QtWidgets.QSlider(window)
#slider.setText('This is a Slider')


window.show()
sys.exit(app.exec())



if __name__ == '__main__':
    n = len(sys.argv)
    print("Total arguments passed:", n)
    for i in range(1, n):
        print(i," ",sys.argv[i])

    if n < 3:
        print("Correct Usage compare2GPXroutes.py filename1.gpx filename2.gpx")

    if n >= 3:
        inputfilename1 = sys.argv[1];
        inputfilename2 = sys.argv[2];
        compare2GPXroutes(inputfilename1,inputfilename2)
